"""HOPS/AWE + least-squares PINN: physics-informed summation of the joint (eps, delta) expansion.

HOPS/AWE (Kehoe & Nicholls 2024) computes, once per frequency band, the coefficient fields of the joint
expansion
    u(x, z'; eps, delta) = sum_{n<=N} sum_{m<=M} u_{n,m}(x, z') eps^n delta^m        (same for w)
on the TFE grid (x_j, z'_l), and then SUMS the series (Taylor or Pade) at every (eps, omega = omega_bar
(1 + delta)) of the reflectivity map.  The summation is where AWE loses accuracy: the series is truncated
at (N, M) (the MATLAB code even sums only to half order), it converges only inside a polydisc limited by
the Rayleigh/Wood singularities, and Pade may produce spurious poles.

The hybrid keeps HOPS' coefficient fields but replaces the fixed weights eps^n delta^m by the weights that
make the field satisfy the governing equations (6a)-(6f) best AT THAT (eps, omega):  the coefficient fields
are the hidden-layer features of a least-squares interface PINN (see PINN_HOPS/pinn_hops/lsq_pinn.py),
    u_theta = sum_j c^u_j phi^u_j,   w_theta = sum_j c^w_j phi^w_j,
    {phi^u_j} = {u_{n,m}} (HOPS/AWE basis)  [+ optional random Fourier x tanh/sin features in (x, zeta)],
and c = argmin (PINN loss) is one complex linear least-squares solve.  The residual rows are the exact
equations at (eps, omega): Helmholtz in the physical coordinates (chain rule through the TFE map, which
depends on eps), the interface conditions (6c, 6d) with the exact incident data, and the exact DtN
transparent conditions (6e, 6f) by FFT.  The Taylor weights are one admissible choice of c, so the
physics-informed sum always has a PDE loss no larger than the AWE Taylor sum.

By-product: the PINN loss of the AWE Taylor sum itself is an a-posteriori error indicator for the
AWE reflectivity map (computable without any reference solution).
"""
import os
import sys
import time
import warnings

import numpy as np
import scipy.linalg as sla

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HOPS_PY = os.environ.get('HOPS_PYTHON', os.path.join(os.path.dirname(HERE), 'HOPS_Python'))
PINN_PY = os.environ.get('PINN_HOPS', os.path.join(os.path.dirname(HERE), 'PINN_HOPS'))
for _p in (HOPS_PY, PINN_PY):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import refl_map as rm                                   # noqa: E402  (HOPS_Python)
from hops import cheb                                   # noqa: E402
from pinn_hops.problem import Grating2D, outgoing_sqrt  # noqa: E402  (PINN_HOPS)
from pinn_hops.lsq_pinn import Features                 # noqa: E402

TWO_PI = 2 * np.pi


def spectral_dx(f, order=1):
    p = np.fft.fftfreq(f.shape[0], 1.0 / f.shape[0])
    sh = (-1,) + (1,) * (f.ndim - 1)
    return np.fft.ifft((1j * p.reshape(sh)) ** order * np.fft.fft(f, axis=0), axis=0)


class AWEBand:
    """One refl_map.py frequency window, with the HOPS/AWE coefficient fields as a basis on the TFE grid.

    scenario/q/overrides are passed to refl_map.run (so the AWE numbers are EXACTLY those of refl_map.py).
    """

    def __init__(self, scenario, q=1, N_Eps=21, N_delta=21, **over):
        t0 = time.time()
        res, info = rm.run(scenario, qq=(q,), N_Eps=N_Eps, N_delta=N_delta, verbose=False, workers=1,
                           keep_fields=True, **over)
        self.t_hops = time.time() - t0
        if len(res) != 1:
            raise ValueError('use a scenario with one window per band (windows="paper", no max_delta)')
        r = res[0]
        self.r, self.info = r, info
        self.Nx, self.Nz, self.N, self.M = info['Nx'], info['Nz'], info['N'], info['M']
        self.a, self.b = info['a'], info['b']
        self.n_u, self.n_w = complex(r['n_u']), complex(r['n_w'])
        self.mode, self.alpha_bar, self.omega_bar = info['mode'], info['alpha'], r['omega_bar']
        self.tau2 = 1.0 + 0j if self.mode == 'TE' else (self.n_u / self.n_w) ** 2
        self.x = TWO_PI * np.arange(self.Nx) / self.Nx
        self.f = np.asarray(info['f'], float)
        self.fx = np.real(spectral_dx(self.f.astype(complex)))
        self.fxx = np.real(spectral_dx(self.f.astype(complex), 2))
        D, t = cheb(self.Nz)
        self.zeta = 0.5 * (t + 1)                   # zeta in [0, 1]: upper 1 = top, lower 1 = interface
        self.Dzeta = 2 * D
        # AWE maps of refl_map.py
        self.Eps, self.delta, self.omega = r['Eps'], r['delta'], r['omega']
        self.R_awe, self.T_awe, self.D_awe = np.real(r['ru']), np.real(r['rl']), np.real(r['ee'])
        self.R_flat = np.real(r['ru_flat'])
        # coefficient fields -> basis (Nx, Nz+1, K), K = (M+1)(N+1); column k = (m, n)
        K = (self.M + 1) * (self.N + 1)
        self.nm = [(n, m) for m in range(self.M + 1) for n in range(self.N + 1)]
        self.basis = {}
        for layer, key in (('u', 'u_n_m'), ('w', 'w_n_m')):
            B = r[key].reshape(self.Nx, self.Nz + 1, K)
            bad = ~np.all(np.isfinite(B), axis=(0, 1))
            if bad.any():                   # an order that overflowed (e.g. a Wood anomaly at omega_bar)
                warnings.warn(f'{bad.sum()} non-finite AWE coefficient fields set to zero')
                B = np.where(bad[None, None, :], 0, B)
            self.basis[layer] = self._derivs(B)

    def _derivs(self, B):
        Bx = spectral_dx(B)
        Dz = self.Dzeta
        return dict(f=B, x=Bx, z=np.einsum('lk,jkq->jlq', Dz, B), xx=spectral_dx(B, 2),
                    zz=np.einsum('lk,jkq->jlq', Dz @ Dz, B), xz=np.einsum('lk,jkq->jlq', Dz, Bx))

    def taylor_weights(self, eps, delta):
        return np.array([eps ** n * delta ** m for (n, m) in self.nm], dtype=complex)

    def grating(self, eps, omega):
        return Grating2D(n_u=self.n_u, n_w=self.n_w, omega=omega, eps=eps,
                         alpha=self.alpha_bar * omega / self.omega_bar, a=self.a, b=self.b, mode=self.mode)


class RandomFeatures:
    """LSQ-PINN enrichment: separable Fourier x random-activation features in (x, zeta), on the TFE grid."""

    def __init__(self, band, K=10, nz=8, activation='sin', rz=3.0, seed=0):
        rng = np.random.default_rng(seed)
        X, Z = np.meshgrid(band.x, band.zeta, indexing='ij')
        self.basis = {}
        for layer in ('u', 'w'):
            F = Features(nz * (2 * K + 1), K, activation, 'separable', (0.0, 1.0), rz=rz, rng=rng)
            d = F(X.ravel(), Z.ravel(), 2)
            self.basis[layer] = {k: v.reshape(X.shape + (F.n,)) for k, v in d.items()}


class PISum:
    """Physics-informed summation of an AWEBand (optionally enriched with random features)."""

    def __init__(self, band, enrich=None, lam_bc=10.0, lam_if=10.0, rcond=1e-13, use_awe=True, compress=None,
                 ridge=None):
        """compress: drop directions of the (column-normalised) AWE basis whose singular values (of the
        field values on the grid) are below compress * s_max -- a POD of the HOPS basis, fewer unknowns."""
        self.B, self.enrich, self.ridge = band, enrich, ridge      # ridge: Tikhonov instead of SVD truncation
        self.lam_bc, self.lam_if, self.rcond = lam_bc, lam_if, rcond
        blocks = []
        if use_awe and compress:
            blocks.append(self._pod(band.basis, compress))
        elif use_awe:
            blocks.append(band.basis)
        if enrich is not None:
            blocks.append(enrich.basis)
        self.feat = {layer: {k: np.concatenate([b[layer][k] for b in blocks], axis=-1)
                             for k in ('f', 'x', 'z', 'xx', 'zz', 'xz')} for layer in ('u', 'w')}
        self.nfeat = {layer: self.feat[layer]['f'].shape[-1] for layer in ('u', 'w')}
        self.n_awe = band.basis['u']['f'].shape[-1] if (use_awe and not compress) else 0
        self.has_awe = use_awe

    @staticmethod
    def _pod(basis, tol):
        out = {}
        for layer, d in basis.items():
            F = d['f'].reshape(-1, d['f'].shape[-1])
            cn = np.linalg.norm(F, axis=0)
            cn[cn == 0] = 1
            _, s, Vh = np.linalg.svd(F / cn, full_matrices=False)
            V = (Vh[s > tol * s[0]].conj().T) / cn[:, None]
            out[layer] = {k: v @ V for k, v in d.items()}
        return out

    # ---- physical derivatives via the TFE map ---------------------------------------------------
    def _phys(self, layer, eps, feat=None):
        B = self.B
        g, gx, gxx = eps * B.f, eps * B.fx, eps * B.fxx
        g, gx, gxx = g[:, None], gx[:, None], gxx[:, None]
        z = B.zeta[None, :]
        if layer == 'u':
            L = B.a - g
            zx = gx * (z - 1) / L
            zxx = (gxx * (z - 1) + 2 * gx ** 2 * (z - 1) / L) / L
        else:
            L = g + B.b
            zx = -z * gx / L
            zxx = -z * (gxx - 2 * gx ** 2 / L) / L
        F = (feat or self.feat)[layer]
        e = lambda v: v[:, :, None]
        return dict(f=F['f'], x=F['x'] + F['z'] * e(zx), z=F['z'] / e(L),
                    xx=F['xx'] + 2 * F['xz'] * e(zx) + F['zz'] * e(zx ** 2) + F['z'] * e(zxx),
                    zz=F['zz'] / e(L ** 2))

    def assemble(self, eps, omega, feat=None):
        """rows A c = r of the PINN loss at (eps, omega) for the features (default: all of this PISum)"""
        B = self.B
        G = B.grating(eps, omega)
        al = G.alpha
        k2 = {'u': complex(G.k_u) ** 2, 'w': complex(G.k_w) ** 2}
        gu, t2 = complex(G.gamma_u), complex(G.tau2)
        Nx, Nz = B.Nx, B.Nz
        feat = feat or self.feat
        nu, nw = feat['u']['f'].shape[-1], feat['w']['f'].shape[-1]
        P = {layer: self._phys(layer, eps, feat) for layer in ('u', 'w')}
        rows, rhs = [], []

        def add(Au, Aw, r, w):
            m = len(r)
            Au = np.zeros((m, nu)) if Au is None else Au
            Aw = np.zeros((m, nw)) if Aw is None else Aw
            rows.append(w * np.hstack([Au, Aw]))
            rhs.append(w * r)

        for layer in ('u', 'w'):                                     # (6a), (6b): interior nodes
            d = P[layer]
            A = d['xx'] + d['zz'] + 2j * al * d['x'] + (k2[layer] - al ** 2) * d['f']
            A = A[:, 1:Nz, :].reshape(-1, A.shape[-1])
            add(A if layer == 'u' else None, A if layer == 'w' else None, np.zeros(len(A), complex),
                1 / (1 + abs(k2[layer])) / np.sqrt(len(A)))
        # interface: upper l = Nz (zeta = 0), lower l = 0 (zeta = 1)
        iu, iw = Nz, 0
        g, gx = eps * B.f, eps * B.fx
        inc = np.exp(-1j * gu * g)
        zeta_d, psi = -inc, (1j * gu + 1j * al * gx) * inc
        du = {k: v[:, iu, :] for k, v in P['u'].items()}
        dw = {k: v[:, iw, :] for k, v in P['w'].items()}
        add(du['f'], -dw['f'], zeta_d, np.sqrt(self.lam_bc / Nx))
        Nu = du['z'] - gx[:, None] * du['x'] - 1j * al * gx[:, None] * du['f']
        Nw = dw['z'] - gx[:, None] * dw['x'] - 1j * al * gx[:, None] * dw['f']
        add(Nu, -t2 * Nw, psi, np.sqrt(self.lam_if / Nx) / (1 + abs(gu)))
        # transparent conditions: upper l = 0 (z = a), lower l = Nz (z = -b)
        for layer, l, sgn in (('u', 0, 1), ('w', Nz, -1)):
            d = P[layer]
            gp = G.gamma_p(Nx, layer)
            T = np.fft.ifft(sgn * 1j * gp[:, None] * np.fft.fft(d['f'][:, l, :], axis=0), axis=0)
            A = d['z'][:, l, :] - T
            add(A if layer == 'u' else None, A if layer == 'w' else None, np.zeros(Nx, complex),
                np.sqrt(self.lam_bc / Nx) / (1 + abs(k2[layer]) ** 0.5))
        return np.vstack(rows), np.concatenate(rhs), G

    def solve(self, eps, omega):
        A, r, G = self.assemble(eps, omega)
        cn = np.linalg.norm(A, axis=0)
        cn[cn == 0] = 1
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if self.ridge:          # smooth in the data (needed for finite-difference Jacobians)
                n = A.shape[1]
                c, _, rank, sv = sla.lstsq(np.vstack([A / cn, self.ridge * np.linalg.norm(r) * np.eye(n)]),
                                           np.concatenate([r, np.zeros(n)]), lapack_driver='gelsd')
            else:
                c, _, rank, sv = sla.lstsq(A / cn, r, cond=self.rcond, lapack_driver='gelsd')
        c = c / cn
        nu = self.nfeat['u']
        out = self._energy(G, c[:nu], c[nu:])
        out.update(loss=float(np.sum(np.abs(A @ c - r) ** 2)), rank=int(rank), c=c)
        if self.has_awe:                   # a-posteriori indicator: PINN loss of the AWE Taylor sum
            out.update(self.indicator(eps, omega))
        return out

    def taylor_fields(self, eps, omega):
        """the AWE Taylor sum (all orders n <= N, m <= M) as a one-column feature set"""
        tw = self.B.taylor_weights(eps, omega / self.B.omega_bar - 1)
        return {layer: {k: (v @ tw)[..., None] for k, v in self.B.basis[layer].items()} for layer in ('u', 'w')}

    def indicator(self, eps, omega):
        """PINN loss of the AWE Taylor sum and its R, T, D -- cheap (no least-squares solve)"""
        F = self.taylor_fields(eps, omega)
        A, r, G = self.assemble(eps, omega, F)
        c = np.ones(2, complex)
        out = dict(loss_awe_taylor=float(np.sum(np.abs(A @ c - r) ** 2)))
        R, T, D = G.energy(F['u']['f'][:, 0, 0], F['w']['f'][:, -1, 0])
        out.update(R_taylor=R, T_taylor=T, D_taylor=D)
        return out

    def _energy(self, G, cu, cw):
        ut = self.feat['u']['f'][:, 0, :] @ cu                      # u(x, a)
        wb = self.feat['w']['f'][:, -1, :] @ cw                     # w(x, -b)
        R, T, D = G.energy(ut, wb)
        return dict(R=R, T=T, D=D)

    def map(self, Eps=None, delta=None, verbose=False, adaptive_tol=None):
        """physics-informed map on the band's (eps, delta) grid.
        adaptive_tol: solve the least-squares problem only where the PINN loss of the AWE Taylor sum
        exceeds adaptive_tol; elsewhere keep the (full-order) Taylor sum."""
        B = self.B
        Eps = B.Eps if Eps is None else Eps
        delta = B.delta if delta is None else delta
        keys = ('R', 'T', 'D', 'loss', 'loss_awe_taylor', 'R_taylor', 'D_taylor')
        out = {k: np.full((len(Eps), len(delta)), np.nan) for k in keys}
        t0 = time.time()
        for i, e in enumerate(Eps):
            for j, dl in enumerate(delta):
                om = float(B.omega_bar * (1 + dl))
                s = None
                if adaptive_tol is not None:
                    s = self.indicator(float(e), om)
                    if s['loss_awe_taylor'] <= adaptive_tol:
                        s.update(R=s['R_taylor'], T=s['T_taylor'], D=s['D_taylor'], loss=s['loss_awe_taylor'])
                        out.setdefault('n_solved', 0)
                    else:
                        s = None
                if s is None:
                    s = self.solve(float(e), om)
                    out['n_solved'] = out.get('n_solved', 0) + 1
                for k in keys:
                    if k in s:
                        out[k][i, j] = s[k]
            if verbose:
                print(f'  eps {e:.3f} done ({time.time() - t0:.0f} s)', flush=True)
        out['time'] = time.time() - t0
        return out
