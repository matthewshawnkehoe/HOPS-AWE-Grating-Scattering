"""TFE field recursion and DNOs for one layer in 3D (analogue of field_tfe_helmholtz_m_and_n[_lf].m,
dno_tfe_helmholtz_m_and_n[_lf].m and hops/operators.py::_layer_run).

Transformed Field Expansions, doubly periodic interface z = g = eps f(x, y)
------------------------------------------------------------------------------
Upper layer  {g < z < a}:  z' = a (z - g)/(a - g).   Lower layer {-b < z < g}:  z' = b (z - g)/(b + g).
Multiplying the transformed Helmholtz equation by (a - g)^2/a^2 (resp. (b + g)^2/b^2) gives,
for u~ = e^{-i(alpha x + beta y)} u,

  Delta' u + 2 i (alpha d_x' + beta d_y') u + gamma^2 u = F(eps, delta)

with (h = a - z' upper, h = b + z' lower; upper / lower coefficient)

  A1_xx = A1_yy = -2f/a | 2f/b        A2_xx = A2_yy = f^2/a^2 | f^2/b^2
  A1_xz = -h f_x/a | -h f_x/b         A2_xz = h f f_x/a^2 | -h f f_x/b^2
  A1_yz = -h f_y/a | -h f_y/b         A2_yz = h f f_y/a^2 | -h f f_y/b^2
  A2_zz = h^2 |grad f|^2/a^2 | h^2 |grad f|^2/b^2
  B1_x  = f_x/a | -f_x/b,  B1_y = f_y/a | -f_y/b
  B2_x  = -f f_x/a^2 | -f f_x/b^2,  B2_y = -f f_y/a^2 | -f f_y/b^2,  B2_z = -h|grad f|^2/a^2 | h|grad f|^2/b^2
  S1 = -2f/a | 2f/b,   S2 = f^2/a^2 | f^2/b^2

(the 2D coefficients of fields.py with f_x -> (f_x, f_y); there is no A_xy term because the
change of variables only couples x and y to z').  The right-hand side of order (n, m) is

  F_{n,m} = - div'( A1 grad' u_{n-1,m} ) - B1.grad' u_{n-1,m} - gamma^2 S1 u_{n-1,m}
            - 2i S1 (alpha u_x + beta u_y)_{n-1,m} - 2i (alpha A1_xz + beta A1_yz) d_z' u_{n-1,m}
            - (same with A2, B2, S2 at n-2)
            - [2 gamma^2 + 2i(alpha d_x + beta d_y)] (u_{n,m-1} + S1 u_{n-1,m-1} + S2 u_{n-2,m-1})
            - 2i (alpha A1_xz + beta A1_yz) d_z' u_{n-1,m-1} - 2i (alpha A2_xz + beta A2_yz) d_z' u_{n-2,m-1}
            - gamma^2 (u_{n,m-2} + S1 u_{n-1,m-2} + S2 u_{n-2,m-2}),

where div'(A grad' u) = d_x'(A_xx u_x + A_xz u_z) + d_y'(A_yy u_y + A_yz u_z)
                        + d_z'(A_xz u_x + A_yz u_y + A_zz u_z).
The 2 i alpha A_xz d_z' terms are the 3D form of config.ALPHA_FIX item 4 (always on in 3D).
Boundary conditions per order (identical to 2D, with the 2D Fourier multipliers):
  interface      u_{n,m}(z'=0) = xi_{n,m}
  artificial     d_z' u - T_0 u = sum_{r<m} T_{m-r} u_{n,r} - (f/a) sum_{r<=m} T_{m-r} u_{n-1,r}   (upper)
DNOs (N = (-g_x, -g_y, 1)):
  G_{n,m} = -u_z + (f_x u_x + f_y u_y)_{n-1} + (f/a) G_{n-1} - (f/a)(f_x u_x + f_y u_y)_{n-2} - |grad f|^2 u_z,{n-2}
  J_{n,m} =  w_z - (f_x w_x + f_y w_y)_{n-1} - (f/b) J_{n-1} - (f/b)(f_x w_x + f_y w_y)_{n-2} + |grad f|^2 w_z,{n-2}
"""
import hashlib
import numpy as np

from .grid import fft2, ifft2
from .expansions import T_dno_3d
from hops.bvp import bvp_matrix

_SETUP_CACHE = {}


def layer_setup(upper, f, f_x, f_y, kx, ky, alphap, betap, gammap, alpha, beta, gamma, Dz, L, Nz, M):
    """Data-independent pieces of one layer (coefficients, inverted collocation matrices,
    DNO symbol expansion).  Cached on all inputs."""
    key = hashlib.sha1(b''.join(np.ascontiguousarray(np.asarray(v, dtype=complex)).tobytes() for v in
                                (f, f_x, f_y, kx, ky, alphap, betap, gammap,
                                 [alpha, beta, gamma, L, Nz, M, upper], Dz))).hexdigest()
    if key in _SETUP_CACHE:
        return _SETUP_CACHE[key]
    Nx, Ny = gammap.shape
    k2 = alphap[0] ** 2 + betap[0] ** 2 + gammap[0, 0] ** 2
    ll = np.arange(Nz + 1)
    tz = np.cos(np.pi * ll / Nz)
    z_min, z_max = (0.0, L) if upper else (-L, 0.0)
    z = ((z_max - z_min) / 2.0) * (tz - 1.0) + z_max
    D = (2.0 / L) * Dz
    D2 = D @ D
    F = f[:, :, None, None]
    FX = f_x[:, :, None, None]
    FY = f_y[:, :, None, None]
    G2 = FX ** 2 + FY ** 2
    one = np.ones((1, 1, Nz + 1, 1))
    gp2 = (k2 - alphap[:, None] ** 2 - betap[None, :] ** 2).reshape(-1)
    gpf = gammap.reshape(-1)
    if upper:
        h = (L - z)[None, None, :, None]
        s = +1.0
        A = bvp_matrix(1.0, 0.0, gp2, Nx * Ny, np.eye(Nz + 1), D, D2, D[0, :], D[-1, :],
                       top_n=1.0, top_d=-1j * gpf, bot_n=0.0, bot_d=1.0)
        ell_bc, ell_if, sgn_T = 0, Nz, +1.0
        phase = np.exp(1j * gammap[:, :, None] * z[None, None, :])
    else:
        h = (L + z)[None, None, :, None]
        s = -1.0
        A = bvp_matrix(1.0, 0.0, gp2, Nx * Ny, np.eye(Nz + 1), D, D2, D[0, :], D[-1, :],
                       top_n=0.0, top_d=1.0, bot_n=1.0, bot_d=1j * gpf)
        ell_bc, ell_if, sgn_T = Nz, 0, -1.0
        phase = np.exp(-1j * gammap[:, :, None] * z[None, None, :])
    # s = +1 upper, -1 lower (see the table in the module docstring)
    C = dict(
        A1_xx=-s * (2.0 / L) * F * one, A1_xz=-(1.0 / L) * h * FX, A1_yz=-(1.0 / L) * h * FY,
        A2_xx=(1.0 / L ** 2) * F ** 2 * one, A2_xz=s * (1.0 / L ** 2) * h * F * FX,
        A2_yz=s * (1.0 / L ** 2) * h * F * FY, A2_zz=(1.0 / L ** 2) * h ** 2 * G2,
        B1_x=s * (1.0 / L) * FX * one, B1_y=s * (1.0 / L) * FY * one,
        B2_x=-(1.0 / L ** 2) * F * FX * one, B2_y=-(1.0 / L ** 2) * F * FY * one,
        B2_z=-s * (1.0 / L ** 2) * h * G2,
        S1=-s * (2.0 / L) * F * one, S2=(1.0 / L ** 2) * F ** 2 * one)
    for k in (1, 2):                                    # 2 i (alpha A_xz + beta A_yz) d_z' (item 4)
        C[f'Aab{k}'] = alpha * C[f'A{k}_xz'] + beta * C[f'A{k}_yz']
    Ainv = np.linalg.inv(A)                             # (Nx*Ny, Nz+1, Nz+1)
    T = T_dno_3d(alpha, beta, alphap, betap, gamma, gammap, k2, M)
    S = dict(upper=upper, f=f, f_x=f_x, f_y=f_y, G2=G2[:, :, 0, 0], Nx=Nx, Ny=Ny, Nz=Nz, L=L,
             alpha=alpha, beta=beta, g2=gamma ** 2, phase=phase, C=C, Ainv=Ainv, T=T,
             ell_bc=ell_bc, ell_if=ell_if, sgn_T=sgn_T, DzL=D, M=M, z=z,
             ikx=(1j * kx)[:, None, None, None], iky=(1j * ky)[None, :, None, None])
    if len(_SETUP_CACHE) >= 4:
        _SETUP_CACHE.clear()
    _SETUP_CACHE[key] = S
    return S


class LayerTFE:
    """Order-by-order TFE recursion for one layer, split so that the Dirichlet data of order (n, m)
    may be supplied AFTER the rest of that order is known (needed by the coupled two-layer solver):

        u_{n,m} = u^p_{n,m} + e^{+-i gamma_pq z'} xi^_{n,m}

    where the particular part u^p (zero Dirichlet data) depends only on lower orders and the
    homogeneous part is the exact solution of the order-(0,0) problem (as at order (0,0) in 2D).

      partial(n, m)          -> computes u^p_{n,m}; returns the DNO of order (n,m) WITHOUT the
                                contribution D_{0,0}[xi_{n,m}] (i.e. sum over lower orders)
      complete(n, m, xi_nm)  -> adds the homogeneous part, stores fields/traces, returns
                                (full DNO_{n,m}, trace_{n,m} at the artificial boundary)
    """

    def __init__(self, S, N, M, B=1, keep_volume=False):
        assert M <= S['M']
        self.S, self.N, self.M, self.B, self.keep = S, N, M, B, keep_volume
        Nx, Ny, Nz = S['Nx'], S['Ny'], S['Nz']
        self.shape = (Nx, Ny, Nz + 1, B)
        nlev = N + 1 if keep_volume else 3
        self.Vu = np.zeros((nlev, M + 1) + self.shape, dtype=complex)
        self.Vx = np.zeros_like(self.Vu); self.Vy = np.zeros_like(self.Vu); self.Vz = np.zeros_like(self.Vu)
        self.trh_bc = np.zeros((N + 1, M + 1, Nx, Ny, B), dtype=complex)   # fft2 of trace at artificial bdry
        self.tr_if = np.zeros((N + 1, M + 1, Nx, Ny, B), dtype=complex)
        self.tx_if = np.zeros_like(self.tr_if); self.ty_if = np.zeros_like(self.tr_if)
        self.tz_if = np.zeros_like(self.tr_if)
        self.D = np.zeros((N + 1, M + 1, Nx, Ny, B), dtype=complex)
        self.Tm = S['T'][:, :, :, None]
        self.phase = S['phase'][..., None]
        self.a0 = (S['alpha'] != 0) or (S['beta'] != 0)
        self.s = +1.0 if S['upper'] else -1.0
        self._uhp = None
        self._cleared = -1

    def lev(self, k):
        return k if self.keep else k % 3

    def _rhs_hat(self, n, m):
        S, C, g2, alpha, beta = self.S, self.S['C'], self.S['g2'], self.S['alpha'], self.S['beta']
        Vu, Vx, Vy, Vz, lev, shape = self.Vu, self.Vx, self.Vy, self.Vz, self.lev, self.shape
        P = np.zeros(shape, dtype=complex)          # pointwise terms
        Xx = np.zeros(shape, dtype=complex)         # fluxes differentiated in x
        Xy = np.zeros(shape, dtype=complex)         # fluxes differentiated in y
        Zf = np.zeros(shape, dtype=complex)         # fluxes differentiated in z'
        for k in (1, 2):                            # geometric terms at order (n-k, m)
            if n < k:
                continue
            L_ = lev(n - k)
            u, ux, uy, uz = Vu[L_, m], Vx[L_, m], Vy[L_, m], Vz[L_, m]
            Axx, Axz, Ayz = C[f'A{k}_xx'], C[f'A{k}_xz'], C[f'A{k}_yz']
            Xx += Axx * ux + Axz * uz
            Xy += Axx * uy + Ayz * uz
            Zf += Axz * ux + Ayz * uy
            P += C[f'B{k}_x'] * ux + C[f'B{k}_y'] * uy + g2 * C[f'S{k}'] * u
            if k == 2:
                Zf += C['A2_zz'] * uz
                P += C['B2_z'] * uz
            if self.a0:
                P += 2j * C[f'S{k}'] * (alpha * ux + beta * uy) + 2j * C[f'Aab{k}'] * uz
        for k in (0, 1, 2):                         # frequency terms at (n-k, m-1), (n-k, m-2)
            if n < k:
                continue
            L_ = lev(n - k)
            Sk = 1.0 if k == 0 else C[f'S{k}']
            if m >= 1:
                u, ux, uy, uz = Vu[L_, m - 1], Vx[L_, m - 1], Vy[L_, m - 1], Vz[L_, m - 1]
                P += 2 * g2 * Sk * u
                if self.a0:
                    P += 2j * Sk * (alpha * ux + beta * uy)
                    if k >= 1:
                        P += 2j * C[f'Aab{k}'] * uz
            if m >= 2:
                P += g2 * Sk * Vu[L_, m - 2]
        P += np.matmul(S['DzL'], Zf)
        return -(fft2(P) + S['ikx'] * fft2(Xx) + S['iky'] * fft2(Xy))

    def partial(self, n, m):
        S = self.S
        Nx, Ny, Nz, L, B = S['Nx'], S['Ny'], S['Nz'], S['L'], self.B
        if not self.keep and self._cleared != n:
            lv = n % 3
            self.Vu[lv] = 0.0; self.Vx[lv] = 0.0; self.Vy[lv] = 0.0; self.Vz[lv] = 0.0
            self._cleared = n
        if n == 0 and m == 0:
            self._uhp = None
            return np.zeros((Nx, Ny, B), dtype=complex)
        Fh = self._rhs_hat(n, m)
        Jh = np.zeros((Nx, Ny, B), dtype=complex)            # artificial-boundary data
        Tm, trh = self.Tm, self.trh_bc
        for r in range(m):
            Jh += S['sgn_T'] * Tm[:, :, m - r] * trh[n, r]
        if n >= 1:
            acc = np.zeros((Nx, Ny, B), dtype=complex)
            for r in range(m + 1):
                acc += Tm[:, :, m - r] * trh[n - 1, r]
            Jh -= (1.0 / L) * fft2(S['f'][:, :, None] * ifft2(acc))
        Fh[:, :, S['ell_bc'], :] = Jh
        Fh[:, :, S['ell_if'], :] = 0.0
        self._uhp = np.matmul(S['Ainv'], Fh.reshape(Nx * Ny, Nz + 1, B)).reshape(self.shape)
        # DNO of the particular part (+ lower-order terms)
        uz_if = ifft2(np.matmul(S['DzL'][S['ell_if']], self._uhp))     # (Nx, Ny, B)
        return self._dno(n, m, uz_if)

    def _dno(self, n, m, uz_if):
        S, s, L = self.S, self.s, self.S['L']
        fx, fy, fc, G2 = S['f_x'][..., None], S['f_y'][..., None], S['f'][..., None], S['G2'][..., None]
        G = -s * uz_if
        if n >= 1:
            G = G + s * (fx * self.tx_if[n - 1, m] + fy * self.ty_if[n - 1, m]) + s * (1.0 / L) * fc * self.D[n - 1, m]
        if n >= 2:
            G = G - (1.0 / L) * fc * (fx * self.tx_if[n - 2, m] + fy * self.ty_if[n - 2, m]) \
                - s * G2 * self.tz_if[n - 2, m]
        return G

    def complete(self, n, m, xi_nm=None):
        S = self.S
        ell_if, ell_bc = S['ell_if'], S['ell_bc']
        uh = self._uhp if self._uhp is not None else 0.0
        if xi_nm is not None:
            xi_nm = np.asarray(xi_nm, dtype=complex)
            if xi_nm.ndim == 2:
                xi_nm = xi_nm[:, :, None]
            uh = uh + self.phase * fft2(xi_nm)[:, :, None, :]
        if np.isscalar(uh):
            uh = np.zeros(self.shape, dtype=complex)
        u = ifft2(uh)
        L_ = self.lev(n)
        self.Vu[L_, m] = u
        self.Vx[L_, m] = ifft2(S['ikx'] * uh)
        self.Vy[L_, m] = ifft2(S['iky'] * uh)
        self.Vz[L_, m] = np.matmul(S['DzL'], u)
        self.trh_bc[n, m] = uh[:, :, ell_bc, :]
        self.tr_if[n, m] = u[:, :, ell_if, :]
        self.tx_if[n, m] = self.Vx[L_, m][:, :, ell_if, :]
        self.ty_if[n, m] = self.Vy[L_, m][:, :, ell_if, :]
        self.tz_if[n, m] = self.Vz[L_, m][:, :, ell_if, :]
        self.D[n, m] = self._dno(n, m, self.tz_if[n, m])
        self._uhp = None
        return self.D[n, m], u[:, :, ell_bc, :]

    def trace(self):
        return ifft2(self.trh_bc.transpose(2, 3, 0, 1, 4)).transpose(2, 3, 0, 1, 4)


def layer_run(S, xi, N, M, keep_volume=False):
    """TFE recursion + DNO for one layer with Dirichlet data xi.

    xi : (Nx, Ny) or (Nx, Ny, B)  -> data placed at order (0, 0) only (lean / operator use), or
         (Nx, Ny, M+1, N+1)       -> full interface data xi_{n,m} (field_tfe_helmholtz_* use).
    Returns dict(DNO (N+1, M+1, Nx, Ny, B), trace (N+1, M+1, Nx, Ny, B) at the artificial boundary,
                 vol (N+1, M+1, Nx, Ny, Nz+1, B) if keep_volume).
    """
    xi = np.asarray(xi, dtype=complex)
    full = xi.ndim == 4
    if not full and xi.ndim == 2:
        xi = xi[:, :, None]
    B = 1 if full else xi.shape[2]
    lay = LayerTFE(S, N, M, B, keep_volume)
    for n in range(N + 1):
        for m in range(M + 1):
            lay.partial(n, m)
            if full:
                lay.complete(n, m, xi[:, :, m, n])
            else:
                lay.complete(n, m, xi if (n == 0 and m == 0) else None)
    out = dict(DNO=lay.D, trace=lay.trace())
    if keep_volume:
        out['vol'] = lay.Vu                              # (N+1, M+1, Nx, Ny, Nz+1, B)
    return out


# ----------------------------------------------------------------------------
# MATLAB-style entry points (full interface data -> volume field -> DNO)
# ----------------------------------------------------------------------------
def _to_std(D):
    """(N+1, M+1, Nx, Ny, 1) -> interface layout (Nx, Ny, M+1, N+1)."""
    return np.transpose(D[..., 0], (2, 3, 1, 0))


def field_tfe_helmholtz_3d(xi_n_m, S, N, M):
    """Volume field u_{n,m} (Nx, Ny, Nz+1, M+1, N+1) of a layer with Dirichlet data xi_{n,m}."""
    out = layer_run(S, xi_n_m, N, M, keep_volume=True)
    return np.transpose(out['vol'][..., 0], (2, 3, 4, 1, 0)), out


def dno_tfe_helmholtz_3d(xi_n_m, S, N, M):
    """(DNO_{n,m}, trace_{n,m}) in the interface layout (Nx, Ny, M+1, N+1)."""
    out = layer_run(S, xi_n_m, N, M)
    return _to_std(out['DNO']), _to_std(out['trace'])
