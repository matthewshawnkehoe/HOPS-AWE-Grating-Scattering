"""3D HOPS/AWE + PINN: physics-informed summation (PI-sum) of the crossed-grating HOPS/AWE series.

The best combined method of the 2D studies (HOPS_PINN_Hybrid, HOPS_AWE_PINN), carried over to the
doubly periodic code hops3d / refl_map_3D.py:

  * per frequency window, HOPS/AWE (coupled 3D solver) gives the coefficient fields u_{n,m}(x, y, z'),
    w_{n,m}(x, y, z') of the joint (eps, delta) expansion on the Fourier x Fourier x Chebyshev TFE grid;
  * these fields are the hidden layer of a least-squares interface PINN:
        u_theta = sum_j c^u_j phi^u_j,   w_theta = sum_j c^w_j phi^w_j    (POD of {u_{n,m}: n, m <= order})
  * at each (eps, omega) the weights c minimise the PINN loss of the 3D problem AT THAT (eps, omega):
        Helmholtz  Delta u + 2i(alpha d_x + beta d_y) u + (k^2 - alpha^2 - beta^2) u = 0  in both layers
                   (physical coordinates, chain rule through the eps-dependent TFE map)
        interface  u - w = -exp(-i gamma g),   d_N u - tau^2 d_N w = (i gamma + i alpha g_x + i beta g_y) e^{-i gamma g}
        exact transparent (DtN) conditions at z = a and z = -b (2D FFT)
  * adaptive: the PINN loss of the full-order AWE Taylor sum (cheap) decides where the solve is needed.

New in 3D (needed for speed): the interior Helmholtz rows are multiplied by L^2 (L = a - g or b + g, the
TFE layer thickness), which makes them EXACTLY polynomial in eps (degree 2) and in s = 1 + delta
(degree 2, since k, alpha, beta all scale with omega):
        A_pde(eps, s) = sum_{i, j <= 2} eps^i s^j A_ij.
The column-stacked [A_00 ... A_22] is QR-factorised once per window, so the (Nx Ny (Nz-1))-row PDE block
of every map point collapses EXACTLY to the small triangular factor:  ||A_pde c||^2 = ||R Phi(eps, s) c||^2.
Only the 4 Nx Ny boundary/interface rows are assembled per point.

The interior rows are the TFE equation in exactly the discrete (conservative flux) form hops3d solves
(layer.py), so the HOPS/AWE Taylor sum has a residual at round-off + truncation level: on under-resolved
grids (cos 4x cos 4y on 32 x 32) the expanded chain-rule form would differ from it by aliasing errors.
"""
import os
import sys
import time
import warnings
from functools import partial

import numpy as np
import scipy.linalg as sla

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_ROOT = os.path.dirname(HERE)
_HP = os.environ.get('HOPS_PYTHON', os.path.join(_ROOT, 'HOPS_Python'))
if _HP not in sys.path:
    sys.path.insert(0, _HP)

import refl_map_3D as rm3                          # noqa: E402  (HOPS_Python)
import hops3d as h3                                # noqa: E402
from hops3d.grid import fft2, ifft2, outgoing_sqrt  # noqa: E402
from hops import cheb                              # noqa: E402

DEFAULT_TOL = 1e-14        # indicator threshold (PINN loss of the AWE Taylor sum); see error_analysis.py
DEFAULT_ORDER = 12         # least-squares basis: orders n, m <= 12
COMPRESS = 1e-13           # POD tolerance
LAM_BC = LAM_IF = 10.0     # boundary / interface row weights (as in 2D)


# ---------------------------------------------------------------------------------------------------
# spectral derivatives on the TFE grid
# ---------------------------------------------------------------------------------------------------
def _dxy(B, kx, ky, ox, oy):
    """d^ox/dx^ox d^oy/dy^oy along axes 0, 1 of B (Nx, Ny, ...)."""
    sh = (B.shape[0], B.shape[1]) + (1,) * (B.ndim - 2)
    mult = ((1j * kx[:, None]) ** ox * (1j * ky[None, :]) ** oy).reshape(sh)
    return ifft2(mult * fft2(B))


def derivs(B, kx, ky, Dzeta, cross=False):
    """dict of TFE derivatives of B (Nx, Ny, Nz+1, nc): f, x, y, z, xx, yy, zz  (z = d/dzeta, zeta in [0, 1])."""
    dz = lambda A: np.einsum('lk,ijkq->ijlq', Dzeta, A, optimize=True)
    Bh = fft2(B)
    sh = (B.shape[0], B.shape[1], 1, 1)
    ikx, iky = (1j * kx[:, None] + 0 * ky[None, :]).reshape(sh), (0 * kx[:, None] + 1j * ky[None, :]).reshape(sh)
    Bx, By = ifft2(ikx * Bh), ifft2(iky * Bh)
    Bz = dz(B)
    out = dict(f=B, x=Bx, y=By, z=Bz, xx=ifft2(ikx ** 2 * Bh), yy=ifft2(iky ** 2 * Bh), zz=dz(Bz))
    if cross:
        out.update(xz=dz(Bx), yz=dz(By))
    return out


# ---------------------------------------------------------------------------------------------------
class Window3D:
    """PI-sum of one refl_map_3D frequency window.

    vu, vw: volume fields (Nx, Ny, Nz+1, M+1, N+1) of the coupled solver (keep_volume=True)."""

    def __init__(self, vu, vw, P, n_u, n_w, omega_bar, alpha_bar, beta_bar, N, M,
                 basis_order=DEFAULT_ORDER, compress=COMPRESS, rcond=1e-13, form='discrete'):
        """form: 'discrete' (default) -- interior rows = the TFE equation in hops3d's discrete flux form;
                 'chain' -- the 2D hybrid's form (physical Helmholtz via the chain rule, x L^2), kept for comparison"""
        self.form = form
        t0 = time.time()
        self.N, self.M, self.rcond = N, M, rcond
        self.Nx, self.Ny = P.Nx, P.Ny
        self.Nz = vu.shape[2] - 1
        self.a, self.b, self.tau2 = P.a, P.b, complex(P.tau2)
        self.n_u, self.n_w, self.omega_bar = complex(n_u), complex(n_w), float(omega_bar)
        self.ab, self.bb = float(alpha_bar), float(beta_bar)
        self.kx, self.ky = P.kx, P.ky
        self.f = np.asarray(P.f, float)
        self.fx, self.fy = np.asarray(P.f_x, float), np.asarray(P.f_y, float)
        self.lapf = np.real(_dxy(self.f.astype(complex), self.kx, self.ky, 2, 0) +
                            _dxy(self.f.astype(complex), self.kx, self.ky, 0, 2))
        D, t = cheb(self.Nz)
        self.zeta = 0.5 * (t + 1)                 # ell = 0: zeta = 1;  ell = Nz: zeta = 0
        self.Dzeta = 2 * D
        # gamma_bar^2 = k_bar^2 - alpha_bar^2 - beta_bar^2 per layer
        self.kbar2 = {'u': (self.n_u * omega_bar) ** 2, 'w': (self.n_w * omega_bar) ** 2}
        self.g2bar = {l: self.kbar2[l] - self.ab ** 2 - self.bb ** 2 for l in 'uw'}
        # ---- bases ------------------------------------------------------------------------------
        self.nm = [(n, m) for m in range(M + 1) for n in range(N + 1)]
        K = len(self.nm)
        nb = min(basis_order, N), min(basis_order, M)
        self.sub = np.array([k for k, (n, m) in enumerate(self.nm) if n <= nb[0] and m <= nb[1]])
        self.raw, self.feat, self.nfeat = {}, {}, {}
        for layer, V in (('u', vu), ('w', vw)):
            B = V.reshape(self.Nx, self.Ny, self.Nz + 1, K)
            bad = ~np.all(np.isfinite(B), axis=(0, 1, 2))
            if bad.any():
                B = np.where(bad[None, None, None, :], 0, B)
            self.raw[layer] = B
            Bs = B[..., self.sub]
            F = Bs.reshape(-1, Bs.shape[-1])
            cn = np.linalg.norm(F, axis=0)
            cn[cn == 0] = 1
            _, s, Vh = np.linalg.svd(F / cn, full_matrices=False)
            keep = s > compress * s[0]
            Bp = (F / cn) @ Vh[keep].conj().T
            self.feat[layer] = Bp.reshape(self.Nx, self.Ny, self.Nz + 1, -1)
            self.nfeat[layer] = int(keep.sum())
        self.t_pod = time.time() - t0
        # ---- PDE rows: exact polynomial blocks in (eps, delta), QR once --------------------------
        Su, Sw = P.setups()                              # hops3d layer data (cached: same as the solve)
        self.S = {'u': Su, 'w': Sw}
        self.Rpde, self.npow = {}, {}
        self.bnd = {}
        for layer in 'uw':
            Dd = derivs(self.feat[layer], self.kx, self.ky, self.Dzeta, cross=form == 'chain')
            blocks = self._pde_blocks(layer, Dd)                 # list of ((i, j), array rows x nc)
            self.npow[layer] = [ij for ij, _ in blocks]
            A = np.hstack([Ab for _, Ab in blocks])
            Rq = sla.qr(A, mode='r', overwrite_a=True, check_finite=False)[0]
            self.Rpde[layer] = np.ascontiguousarray(Rq[:min(Rq.shape)])   # economic triangular factor
            del Rq
            self.bnd[layer] = self._boundary_data(layer, Dd)
            del Dd, blocks, A
        nu = self.nfeat['u']
        # energy: Fourier coefficients of u at z = a (ell = 0, upper) and of w at z = -b (ell = Nz, lower)
        self.uhat_top = fft2(self.feat['u'][:, :, 0, :]) / (self.Nx * self.Ny)
        self.what_bot = fft2(self.feat['w'][:, :, -1, :]) / (self.Nx * self.Ny)
        self.n_unknowns = nu + self.nfeat['w']
        self.t_setup = time.time() - t0

    # ---- geometry of the TFE map (per layer) ------------------------------------------------------
    def _geom(self, layer, zeta):
        """c1(zeta), L0, Lf (L = L0 + eps Lf), sp:  u_x = F_x + c1 g_x F_z / L, L^2 zeta_xx = c1 L lap g + 2 sp c1 |grad g|^2"""
        if layer == 'u':
            return zeta - 1.0, self.a, -self.f, 1.0
        return -zeta, self.b, self.f, -1.0

    def _pde_blocks(self, layer, Dd, interior=True):
        if self.form == 'chain':
            return self._pde_blocks_chain(layer, Dd, interior)
        return self._pde_blocks_discrete(layer, Dd, interior)

    def _pde_blocks_chain(self, layer, Dd, interior=True):
        """L^2 x (physical Helmholtz residual via the chain rule) = sum eps^i delta^j A_ij (the 2D hybrid's
        rows).  Analytically equal to the discrete form, but differs by aliasing on under-resolved grids."""
        sl = slice(1, self.Nz) if interior else slice(None)
        D = {k: v[:, :, sl, :] for k, v in Dd.items()}
        c1, L0, Lf, sp = self._geom(layer, self.zeta[sl])
        c1 = c1[None, None, :, None]
        E = lambda v: np.asarray(v)[:, :, None, None]
        fx, fy, lapf, G2, Lf = E(self.fx), E(self.fy), E(self.lapf), E(self.fx ** 2 + self.fy ** 2), E(Lf)
        lapF = D['xx'] + D['yy']
        gfF = fx * D['xz'] + fy * D['yz']
        P0 = [L0 ** 2 * lapF + D['zz'],
              2 * L0 * Lf * lapF + 2 * c1 * L0 * gfF + c1 * L0 * lapf * D['z'],
              Lf ** 2 * lapF + 2 * c1 * Lf * gfF + c1 ** 2 * G2 * D['zz'] + (c1 * Lf * lapf + 2 * sp * c1 * G2) * D['z']]
        g2 = self.g2bar[layer]
        F = D['f']
        L2 = [L0 ** 2 * F, 2 * L0 * Lf * F, Lf ** 2 * F]       # x gamma_bar^2 s^2 = g2 (1 + 2 delta + delta^2)
        out = {}
        for i in range(3):
            out[(i, 0)] = P0[i] + g2 * L2[i]
            out[(i, 1)] = 2 * g2 * L2[i]
            out[(i, 2)] = g2 * L2[i]
        if self.ab != 0 or self.bb != 0:
            AB = self.ab * D['x'] + self.bb * D['y']
            afb = self.ab * fx + self.bb * fy
            P1 = [2j * L0 ** 2 * AB, 2j * (2 * L0 * Lf * AB + c1 * L0 * afb * D['z']),
                  2j * (Lf ** 2 * AB + c1 * Lf * afb * D['z'])]          # x s = 1 + delta
            for i in range(3):
                out[(i, 0)] = out[(i, 0)] + P1[i]
                out[(i, 1)] = out[(i, 1)] + P1[i]
        nrow = self.Nx * self.Ny * D['f'].shape[2]
        w = 1.0 / (1 + abs(self.kbar2[layer])) / np.sqrt(nrow)
        return [(ij, w * A.reshape(nrow, -1)) for ij, A in out.items()]

    def _pde_blocks_discrete(self, layer, Dd, interior=True):
        """The TFE Helmholtz equation EXACTLY as hops3d discretises it (conservative flux form, products
        on the grid, spectral / Chebyshev derivatives of the fluxes -- layer.py::_rhs_hat), summed over all
        orders:  sum_{k, j <= 2} eps^k delta^j E_kj[u] = 0  on the interior Chebyshev nodes.
        Returns [((k, j), E_kj as rows x nc)] (weighted)."""
        S = self.S[layer]
        C, g2, al, be, Lc = S['C'], S['g2'], S['alpha'], S['beta'], S['L']
        a0 = (al != 0) or (be != 0)
        DzL = S['DzL']
        u, ux, uy = Dd['f'], Dd['x'], Dd['y']
        uz, uzz = Dd['z'] / Lc, Dd['zz'] / Lc ** 2               # d/dz' = (1/L) d/dzeta
        kx, ky = self.kx, self.ky
        dx = lambda A: _dxy(A, kx, ky, 1, 0)
        dy = lambda A: _dxy(A, kx, ky, 0, 1)
        dz = lambda A: np.einsum('lk,ijkq->ijlq', DzL, A, optimize=True)
        out = {(0, 0): uzz + Dd['xx'] + Dd['yy'] + 2j * (al * ux + be * uy) + g2 * u}
        for k in (1, 2):
            Axx, Axz, Ayz = C[f'A{k}_xx'], C[f'A{k}_xz'], C[f'A{k}_yz']
            Zf = Axz * ux + Ayz * uy
            P = C[f'B{k}_x'] * ux + C[f'B{k}_y'] * uy + g2 * C[f'S{k}'] * u
            if k == 2:
                Zf = Zf + C['A2_zz'] * uz
                P = P + C['B2_z'] * uz
            if a0:
                P = P + 2j * C[f'S{k}'] * (al * ux + be * uy) + 2j * C[f'Aab{k}'] * uz
            out[(k, 0)] = dx(Axx * ux + Axz * uz) + dy(Axx * uy + Ayz * uz) + dz(Zf) + P
        for k in (0, 1, 2):
            Sk = 1.0 if k == 0 else C[f'S{k}']
            E1 = 2 * g2 * Sk * u
            if a0:
                E1 = E1 + 2j * Sk * (al * ux + be * uy)
                if k >= 1:
                    E1 = E1 + 2j * C[f'Aab{k}'] * uz
            out[(k, 1)] = E1
            out[(k, 2)] = g2 * Sk * u
        sl = slice(1, self.Nz) if interior else slice(None)
        nrow = self.Nx * self.Ny * len(range(self.Nz + 1)[sl])
        w = 1.0 / (1 + abs(self.kbar2[layer])) / np.sqrt(nrow)
        return [(kj, w * np.asarray(A)[:, :, sl, :].reshape(nrow, -1)) for kj, A in out.items()]

    def _boundary_data(self, layer, Dd):
        """values / derivatives at the interface node and the transparent-boundary node"""
        i_if, i_tb = (self.Nz, 0) if layer == 'u' else (0, self.Nz)
        pick = lambda l: {k: Dd[k][:, :, l, :] for k in ('f', 'x', 'y', 'z')}
        return dict(iface=pick(i_if), tbc=pick(i_tb), tbc_hat=fft2(Dd['f'][:, :, i_tb, :]))

    # ---- per-point pieces -------------------------------------------------------------------------
    def _phys_if(self, d, layer, eps):
        """physical u, u_x, u_y, u_z at the interface (zeta = 0 upper / 1 lower: c1 = -1 for both)"""
        _, L0, Lf, _ = self._geom(layer, np.array([0.0]))
        L = (L0 + eps * Lf)[:, :, None]
        gx, gy = (eps * self.fx)[:, :, None], (eps * self.fy)[:, :, None]
        return d['f'], d['x'] - gx * d['z'] / L, d['y'] - gy * d['z'] / L, d['z'] / L

    def _gamma_pq(self, layer, s):
        k2 = self.kbar2[layer] * s ** 2
        ap = self.ab * s + self.kx[:, None]
        bq = self.bb * s + self.ky[None, :]
        return outgoing_sqrt(k2 - ap ** 2 - bq ** 2), np.real(ap ** 2 + bq ** 2) < np.real(k2)

    def boundary_rows(self, eps, s, data=None):
        """rows (A_u, A_w, rhs) of the interface and transparent conditions at (eps, s = omega/omega_bar)"""
        data = data or self.bnd
        Nxy = self.Nx * self.Ny
        al, be = self.ab * s, self.bb * s
        gu = complex(outgoing_sqrt(self.g2bar['u'])) * s
        g = eps * self.f
        gx, gy = eps * self.fx, eps * self.fy
        inc = np.exp(-1j * gu * g)
        rD = (-inc).ravel()
        rN = ((1j * gu + 1j * al * gx + 1j * be * gy) * inc).ravel()
        fl = lambda A: A.reshape(Nxy, -1)
        U = self._phys_if(data['u']['iface'], 'u', eps)
        W = self._phys_if(data['w']['iface'], 'w', eps)
        nrm = lambda Q: Q[3] - gx[:, :, None] * Q[1] - gy[:, :, None] * Q[2] - 1j * (al * gx + be * gy)[:, :, None] * Q[0]
        wD = np.sqrt(LAM_BC / Nxy)
        wN = np.sqrt(LAM_IF / Nxy) / (1 + abs(gu))
        Au = [wD * fl(U[0]), wN * fl(nrm(U))]
        Aw = [-wD * fl(W[0]), -wN * self.tau2 * fl(nrm(W))]
        rhs = [wD * rD, wN * rN]
        for layer, sgn in (('u', 1.0), ('w', -1.0)):
            d = data[layer]
            _, L0, Lf, _ = self._geom(layer, np.array([0.0]))
            L = (L0 + eps * Lf)[:, :, None]
            gp, _ = self._gamma_pq(layer, s)
            T = ifft2(sgn * 1j * gp[:, :, None] * d['tbc_hat'])
            A = fl(d['tbc']['z'] / L - T) * (np.sqrt(LAM_BC / Nxy) / (1 + abs(self.kbar2[layer]) ** 0.5 * s))
            if layer == 'u':
                Au.append(A); Aw.append(None)
            else:
                Au.append(None); Aw.append(A)
            rhs.append(np.zeros(Nxy, complex))
        return Au, Aw, rhs

    def pde_reduced(self, layer, eps, s):
        R = self.Rpde[layer]
        nc = self.nfeat[layer]
        out = np.zeros((R.shape[0], nc), complex)
        d = s - 1.0
        for b, (i, j) in enumerate(self.npow[layer]):
            out += (eps ** i * d ** j) * R[:, b * nc:(b + 1) * nc]
        return out

    def energy(self, s, cu, cw):
        gu, pu = self._gamma_pq('u', s)
        gw, pw = self._gamma_pq('w', s)
        uh = (self.uhat_top @ cu)
        wh = (self.what_bot @ cw)
        g0 = gu[0, 0]
        R = float(np.real(np.sum(np.where(pu, (gu / g0) * np.abs(uh) ** 2, 0))))
        T = float(np.real(np.sum(np.where(pw, self.tau2 * (gw / g0) * np.abs(wh) ** 2, 0))))
        return R, T, 1.0 - R - T

    # ---- the AWE Taylor sum (all orders) and its PINN loss: the indicator -----------------------------
    def taylor_weights(self, eps, delta):
        return np.array([eps ** n * delta ** m for (n, m) in self.nm], dtype=complex)

    def indicator(self, eps, omega):
        s = omega / self.omega_bar
        tw = self.taylor_weights(eps, s - 1)
        loss = 0.0
        data = {}
        for layer in 'uw':
            Fld = (self.raw[layer] @ tw)[..., None]
            Dd = derivs(Fld, self.kx, self.ky, self.Dzeta, cross=self.form == 'chain')
            r = sum((eps ** i * (s - 1) ** j) * A for (i, j), A in self._pde_blocks(layer, Dd))
            loss += float(np.sum(np.abs(r) ** 2))
            data[layer] = self._boundary_data(layer, Dd)
        Au, Aw, rhs = self.boundary_rows(eps, s, data)
        for au, aw, r in zip(Au, Aw, rhs):
            res = -r.copy()
            if au is not None:
                res = res + au[:, 0]
            if aw is not None:
                res = res + aw[:, 0]
            loss += float(np.sum(np.abs(res) ** 2))
        uh = fft2(data['u']['tbc']['f'][:, :, 0]) / (self.Nx * self.Ny)
        wh = fft2(data['w']['tbc']['f'][:, :, 0]) / (self.Nx * self.Ny)
        gu, pu = self._gamma_pq('u', s)
        gw, pw = self._gamma_pq('w', s)
        R = float(np.real(np.sum(np.where(pu, (gu / gu[0, 0]) * np.abs(uh) ** 2, 0))))
        T = float(np.real(np.sum(np.where(pw, self.tau2 * (gw / gu[0, 0]) * np.abs(wh) ** 2, 0))))
        return dict(loss_awe_taylor=loss, R_taylor=R, T_taylor=T, D_taylor=1 - R - T)

    # ---- least-squares weights --------------------------------------------------------------------
    def solve(self, eps, omega):
        s = omega / self.omega_bar
        nu, nw = self.nfeat['u'], self.nfeat['w']
        Ru, Rw = self.pde_reduced('u', eps, s), self.pde_reduced('w', eps, s)
        Au, Aw, rhs = self.boundary_rows(eps, s)
        rows = [np.hstack([Ru, np.zeros((Ru.shape[0], nw))]), np.hstack([np.zeros((Rw.shape[0], nu)), Rw])]
        for au, aw in zip(Au, Aw):
            m = (au if au is not None else aw).shape[0]
            rows.append(np.hstack([au if au is not None else np.zeros((m, nu)),
                                   aw if aw is not None else np.zeros((m, nw))]))
        A = np.vstack(rows)
        r = np.concatenate([np.zeros(Ru.shape[0] + Rw.shape[0], complex)] + rhs)
        cn = np.linalg.norm(A, axis=0)
        cn[cn == 0] = 1
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            c = sla.lstsq(A / cn, r, cond=self.rcond, lapack_driver='gelsd', check_finite=False)[0] / cn
        R, T, D = self.energy(s, c[:nu], c[nu:])
        return dict(R=R, T=T, D=D, loss=float(np.sum(np.abs(A @ c - r) ** 2)), c=c)

    def point(self, eps, omega, tol=DEFAULT_TOL):
        """adaptive: AWE (full-order) Taylor sum where its PINN loss <= tol, else least squares"""
        s = self.indicator(eps, omega)
        if s['loss_awe_taylor'] <= tol:
            s.update(R=s['R_taylor'], T=s['T_taylor'], D=s['D_taylor'], solved=False)
            return s
        s.update(self.solve(eps, omega), solved=True)
        return s

    def fields_xz(self, c, iy=0):
        """x-z' slices (Nx, Nz+1) of u_theta, w_theta at y index iy for weights c (pictures)"""
        nu = self.nfeat['u']
        return self.feat['u'][:, iy] @ c[:nu], self.feat['w'][:, iy] @ c[nu:]


# ---------------------------------------------------------------------------------------------------
# the refl_map_3D window, extended: AWE exactly as refl_map_3D._window + the PI-sum map
# ---------------------------------------------------------------------------------------------------
def window_problem(win, c):
    """(P, n_u, n_w, alpha, beta) of one window exactly as refl_map_3D._window builds them"""
    key, omega_bar, dmax, edges = win
    n_u, n_w = c['n_u_w'][key], c['n_w_w'][key]
    if c['theta'] is not None:
        th, ph = np.deg2rad(c['theta']), np.deg2rad(c['phi'])
        k_u = np.real(n_u) * omega_bar
        alpha, beta = k_u * np.sin(th) * np.cos(ph), k_u * np.sin(th) * np.sin(ph)
    else:
        alpha, beta = c['alpha'], c['beta']
    P = h3.make_problem(c['Nx'], c['Ny'], c['Nz'], c['N'], c['M'], n_u, n_w, omega_bar, alpha, beta, f=c['f'],
                        f_x=c['f_x'], f_y=c['f_y'], a=c['a'], b=c['b'], Mode=c['Mode'])
    return P, n_u, n_w, alpha, beta


def hybrid_window(win, c, opts):
    """Solve one window (AWE, as refl_map_3D._window) and evaluate the PI-sum on its map grid.
    opts: tol, basis_order, keep_window (return the Window3D object), points (extra (eps, omega) list)."""
    key, omega_bar, dmax, edges = win
    N, M, Eps = c['N'], c['M'], c['Eps']
    delta = np.array([0.0]) if c['N_delta'] == 1 else np.linspace(-dmax, dmax, c['N_delta'])
    omega = omega_bar * (1 + delta)
    t0 = time.time()
    P, n_u, n_w, alpha, beta = window_problem(win, c)
    zeta, psi = h3.setup_zeta_psi_n_m_3d(alpha, beta, P.gamma_u_bar, P.f, P.f_x, P.f_y, N, M)
    U, W, ubar, wbar, vu, vw = h3.two_layer_solve_3d_coupled(P, zeta, psi, keep_volume=True)
    t1 = time.time()
    kw = dict(sum_domain=c['sum_domain'])
    args = (P.tau2, ubar, wbar, P.kx, P.ky, alpha, beta, P.gamma_u_bar, P.gamma_w_bar, Eps, delta)
    ee_f, ru_f, rl_f = h3.energy_defect_3d(*args, 0, 0, 1, **kw)
    ee, ru, rl = h3.energy_defect_3d(*args, N, M, 1 if c['Taylor'] else 2,
                                     taylor_full_order=c.get('taylor_full_order', False), **kw)
    t_awe = t1 - t0 + (time.time() - t1)
    out = dict(key=key, omega_bar=omega_bar, delta=delta, omega=omega, lam=2 * np.pi / omega, Eps=Eps,
               ee_awe=ee, ru_awe=ru, rl_awe=rl, ru_flat=ru_f, n_u=n_u, n_w=n_w, alpha=alpha, beta=beta,
               edges=edges, dmax=dmax, t_awe=t_awe)
    t2 = time.time()
    Wn = Window3D(vu, vw, P, n_u, n_w, omega_bar, alpha, beta, N, M, basis_order=opts.get('basis_order', DEFAULT_ORDER),
                  form=opts.get('form', 'discrete'))
    del vu, vw
    tol = opts.get('tol', DEFAULT_TOL)
    if opts.get('map', True):
        ne, nd = len(Eps), len(delta)
        R, T, D, ind, sol = (np.zeros((ne, nd)) for _ in range(5))
        for i, e in enumerate(Eps):
            for j, om in enumerate(omega):
                s = Wn.point(float(e), float(om), tol)
                R[i, j], T[i, j], D[i, j], ind[i, j], sol[i, j] = s['R'], s['T'], s['D'], s['loss_awe_taylor'], s['solved']
        out.update(ru=R + 0j, rl=T + 0j, ee=D + 0j, indicator=ind, solved=sol.astype(bool))
    pts = opts.get('points')
    if pts:
        out['points'] = [dict(eps=e, omega=om, **{k: v for k, v in Wn.point(e, om, -1.0).items() if k != 'c'})
                         for e, om in pts]
    out.update(t_hybrid=time.time() - t2, t_setup=Wn.t_setup, n_unknowns=Wn.n_unknowns,
               n_pde_rows=int(sum(Wn.Rpde[l].shape[0] for l in 'uw')))
    if opts.get('keep_window'):
        if opts['keep_window'] == 'light':          # rendering only needs the least-squares basis
            Wn.raw = None
        out.update(W=Wn, ubar_n_m=ubar, wbar_n_m=wbar, W_n_m=W, U_n_m=U, gamma_u_bar=P.gamma_u_bar,
                   gamma_w_bar=P.gamma_w_bar, tau2=P.tau2, kx=P.kx, ky=P.ky, iy=0)
    if c['verbose']:
        print(f"  {key}: omega_bar {omega_bar:.3f}: AWE {t_awe:.1f} s, hybrid {out['t_hybrid']:.0f} s "
              f"({Wn.n_unknowns} unknowns, setup {Wn.t_setup:.0f} s)", flush=True)
    return out


def run(scenario, qq=None, N_Eps=15, N_delta=15, tol=DEFAULT_TOL, basis_order=DEFAULT_ORDER, workers=1,
        verbose=False, keep_window=False, map=True, form='discrete', points=None, **over):
    """refl_map_3D.run with the PI-sum: results in refl_map_3D format (ru, rl, ee, RR = hybrid) plus the
    AWE copies ru_awe, rl_awe, ee_awe, RR_awe, the indicator and the 'solved' mask."""
    opts = dict(tol=tol, basis_order=basis_order, keep_window=keep_window, map=map, form=form, points=points)
    orig = rm3._window
    rm3._window = partial(hybrid_window, opts=opts)
    t0 = time.time()
    try:
        res, info = rm3.run(scenario, qq=qq, N_Eps=N_Eps, N_delta=N_delta, verbose=verbose, workers=workers, **over)
    finally:
        rm3._window = orig
    for r in res:
        r['RR_awe'] = r['ru_awe'] / r['ru_flat'] if info['relative'] else r['ru_awe']
        if map:
            r['RR'] = r['ru'] / r['ru_flat'] if info['relative'] else r['ru']
    info.update(t_total=time.time() - t0, tol=tol, basis_order=basis_order,
                t_awe=sum(r['t_awe'] for r in res), t_hybrid=sum(r['t_hybrid'] for r in res))
    return res, info


def awe_view(res):
    """the AWE results in plain refl_map_3D format (for rm3.plot)"""
    return [dict(r, ru=r['ru_awe'], rl=r['rl_awe'], ee=r['ee_awe'], RR=r['RR_awe']) for r in res]
