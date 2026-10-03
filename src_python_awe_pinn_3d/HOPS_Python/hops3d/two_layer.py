"""Two-layer HOPS/AWE solver in 3D (analogue of AInverse.m, two_layer_solve_fast.m and
hops/operators.py).

Order by order in (eps, delta):
    U_{n,m} - W_{n,m}                         = zeta_{n,m}
    G_{0,0}[U_{n,m}] + tau2 J_{0,0}[W_{n,m}]  = -psi_{n,m} - i(alpha g_x + beta g_y)(U - tau2 W)|_{n,m}
                                               - sum_{(r,s)<(n,m)} (G_{n-r,m-s}[U_{r,s}] + tau2 J_{n-r,m-s}[W_{r,s}])
with G_{0,0} = -i gamma^u_pq, J_{0,0} = -i gamma^w_pq (Fourier multipliers).

Solvers
  two_layer_solve_3d          'lean': the contributions G_{.,.}[U_{r,s}] are accumulated as soon as
                              U_{r,s} is known (same arithmetic as two_layer_solve_fast.m, but memory
                              ~ one volume field); upper and lower layers run in two threads.
  two_layer_solve_3d_operator the linear maps U -> G_{p,r}[U], T_{p,r}[U] are built once on the
                              Fourier-mode basis (NxNy columns) -- faster for small grids and
                              large N, M (flop ratio ~ 4 Nx Ny / (N M)); memory ~ (N+1)(M+1)(NxNy)^2.
  two_layer_solve_3d_coupled  (default) one interleaved TFE recursion per layer: (N+1)(M+1)
                              order-solves per layer instead of ~(NM)^2/4 -- see its docstring.
  two_layer_solve_3d_auto     picks lean or operator.
"""
import time
import numpy as np

from .grid import fft2, ifft2
from .layer import layer_setup, layer_run, LayerTFE


def AInverse_3d(Q, R, gammap_u, gammap_w, tau2):
    """Invert the flat-interface 2x2 multiplier (per Fourier mode (p, q)); returns (U, W)."""
    Qh, Rh = fft2(Q), fft2(R)
    det = tau2 * gammap_w + gammap_u
    return ifft2(((tau2 * gammap_w) * Qh + 1j * Rh) / det), ifft2(((-gammap_u) * Qh + 1j * Rh) / det)


def phase_correction_3d(U, W, m, n, alpha_bar, beta_bar, f_x, f_y, tau2):
    """i (alpha(delta) g_x + beta(delta) g_y)(U - tau^2 W) at order (n, m)  (2D: config.ALPHA_FIX item 3)."""
    if n < 1 or (alpha_bar == 0 and beta_bar == 0):
        return 0.0
    X = U[:, :, m, n - 1] - tau2 * W[:, :, m, n - 1]
    if m >= 1:
        X = X + U[:, :, m - 1, n - 1] - tau2 * W[:, :, m - 1, n - 1]
    return 1j * (alpha_bar * f_x + beta_bar * f_y) * X


class Problem3D:
    """Everything that defines one frequency window (bundles the long 2D argument lists)."""

    def __init__(self, f, f_x, f_y, kx, ky, alphap, betap, gammap_u, gammap_w, alpha_bar, beta_bar,
                 gamma_u_bar, gamma_w_bar, Dz, a, b, Nz, N, M, tau2):
        self.__dict__.update(locals())
        del self.__dict__['self']
        self.Nx, self.Ny = gammap_u.shape

    def setups(self):
        Su = layer_setup(True, self.f, self.f_x, self.f_y, self.kx, self.ky, self.alphap, self.betap,
                         self.gammap_u, self.alpha_bar, self.beta_bar, self.gamma_u_bar, self.Dz, self.a,
                         self.Nz, self.M)
        Sw = layer_setup(False, self.f, self.f_x, self.f_y, self.kx, self.ky, self.alphap, self.betap,
                         self.gammap_w, self.alpha_bar, self.beta_bar, self.gamma_w_bar, self.Dz, self.b,
                         self.Nz, self.M)
        return Su, Sw


def two_layer_solve_3d(P, zeta, psi, verbose=False, threads=True):
    """Lean solver.  zeta, psi: (Nx, Ny, M+1, N+1).  Returns U, W, ubar, wbar (same layout)."""
    N, M, tau2 = P.N, P.M, P.tau2
    Su, Sw = P.setups()
    shp = zeta.shape[:2] + (M + 1, N + 1)
    U = np.zeros(shp, dtype=complex); W = np.zeros(shp, dtype=complex)
    acc = np.zeros(shp, dtype=complex)
    ubar = np.zeros(shp, dtype=complex); wbar = np.zeros(shp, dtype=complex)
    pool = None
    if threads:
        from concurrent.futures import ThreadPoolExecutor
        pool = ThreadPoolExecutor(1)
    t0 = time.time()
    for n in range(N + 1):
        for m in range(M + 1):
            R = (-psi[:, :, m, n] - acc[:, :, m, n]
                 - phase_correction_3d(U, W, m, n, P.alpha_bar, P.beta_bar, P.f_x, P.f_y, tau2))
            U[:, :, m, n], W[:, :, m, n] = AInverse_3d(zeta[:, :, m, n], R, P.gammap_u, P.gammap_w, tau2)
            if pool is not None:
                fu = pool.submit(layer_run, Su, U[:, :, m, n], N - n, M - m)
                ow = layer_run(Sw, W[:, :, m, n], N - n, M - m)
                ou = fu.result()
            else:
                ou = layer_run(Su, U[:, :, m, n], N - n, M - m)
                ow = layer_run(Sw, W[:, :, m, n], N - n, M - m)
            GJ = (ou['DNO'] + tau2 * ow['DNO'])[..., 0]          # (N-n+1, M-m+1, Nx, Ny)
            GJ[0, 0] = 0.0
            acc[:, :, m:, n:] += np.transpose(GJ, (2, 3, 1, 0))
            ubar[:, :, m:, n:] += np.transpose(ou['trace'][..., 0], (2, 3, 1, 0))
            wbar[:, :, m:, n:] += np.transpose(ow['trace'][..., 0], (2, 3, 1, 0))
        if verbose:
            print(f'  two_layer_solve_3d: n = {n}/{N} ({time.time() - t0:.1f} s)', flush=True)
    if pool is not None:
        pool.shutdown()
    return U, W, ubar, wbar


def two_layer_solve_3d_coupled(P, zeta, psi, verbose=False, keep_volume=False):
    """Coupled solver (new; default).  Runs ONE TFE recursion per layer on the full data
    U = sum U_{r,s} eps^r delta^s, interleaved with the interface solve:

        at order (n, m):  u^p_{n,m}, w^p_{n,m} from lower orders  ->  their DNOs equal
        sum_{(r,s)<(n,m)} G_{n-r,m-s}[U_{r,s}] (and J[W])  ->  AInverse gives U_{n,m}, W_{n,m}
        ->  u_{n,m} = u^p_{n,m} + e^{i gamma z} U^_{n,m}.

    By linearity and shift invariance of the recursion this is algebraically identical to
    two_layer_solve_3d / two_layer_solve_fast.m, but costs (N+1)(M+1) order-solves per layer
    instead of ~(N+1)^2 (M+1)^2 / 4  (~30x fewer for N = M = 10).
    keep_volume=True also returns the volume fields (upper, lower) as (Nx, Ny, Nz+1, M+1, N+1).
    """
    N, M, tau2 = P.N, P.M, P.tau2
    Su, Sw = P.setups()
    lu = LayerTFE(Su, N, M, 1, keep_volume)
    lw = LayerTFE(Sw, N, M, 1, keep_volume)
    shp = zeta.shape[:2] + (M + 1, N + 1)
    U = np.zeros(shp, dtype=complex); W = np.zeros(shp, dtype=complex)
    ubar = np.zeros(shp, dtype=complex); wbar = np.zeros(shp, dtype=complex)
    t0 = time.time()
    for n in range(N + 1):
        for m in range(M + 1):
            Gp = lu.partial(n, m)[..., 0]
            Jp = lw.partial(n, m)[..., 0]
            R = (-psi[:, :, m, n] - Gp - tau2 * Jp
                 - phase_correction_3d(U, W, m, n, P.alpha_bar, P.beta_bar, P.f_x, P.f_y, tau2))
            U[:, :, m, n], W[:, :, m, n] = AInverse_3d(zeta[:, :, m, n], R, P.gammap_u, P.gammap_w, tau2)
            ubar[:, :, m, n] = lu.complete(n, m, U[:, :, m, n])[1][..., 0]
            wbar[:, :, m, n] = lw.complete(n, m, W[:, :, m, n])[1][..., 0]
        if verbose:
            print(f'  two_layer_solve_3d_coupled: n = {n}/{N} ({time.time() - t0:.1f} s)', flush=True)
    if keep_volume:
        vol = lambda l: np.transpose(l.Vu[..., 0], (2, 3, 4, 1, 0))
        return U, W, ubar, wbar, vol(lu), vol(lw)
    return U, W, ubar, wbar


def layer_operators_3d(P):
    """(G, Tu, J, Tw) acting on Fourier coefficients c = fft2(U)/(Nx Ny); each (N+1, M+1, NxNy, NxNy)."""
    Su, Sw = P.setups()
    Nx, Ny = P.Nx, P.Ny
    I = np.eye(Nx * Ny).reshape(Nx, Ny, Nx * Ny)
    modes = ifft2(I) * (Nx * Ny)                                # column k = e^{i (p_k x + q_k y)}
    ops = []
    for S in (Su, Sw):
        o = layer_run(S, modes, P.N, P.M)
        ops += [o['DNO'].reshape(P.N + 1, P.M + 1, Nx * Ny, Nx * Ny),
                o['trace'].reshape(P.N + 1, P.M + 1, Nx * Ny, Nx * Ny)]
    return ops[0], ops[1], ops[2], ops[3]


def two_layer_solve_3d_operator(P, zeta, psi, verbose=False, ops=None):
    N, M, tau2, Nx, Ny = P.N, P.M, P.tau2, P.Nx, P.Ny
    G, Tu, J, Tw = ops if ops is not None else layer_operators_3d(P)
    K = Nx * Ny
    shp = (Nx, Ny, M + 1, N + 1)
    U = np.zeros(shp, dtype=complex); W = np.zeros(shp, dtype=complex)
    Uh = np.zeros((K, M + 1, N + 1), dtype=complex); Wh = np.zeros_like(Uh)
    for n in range(N + 1):
        for m in range(M + 1):
            R = (-psi[:, :, m, n] - phase_correction_3d(U, W, m, n, P.alpha_bar, P.beta_bar, P.f_x, P.f_y, tau2))
            R = R.astype(complex)
            if n > 0 or m > 0:
                rr, ss = np.meshgrid(np.arange(n + 1), np.arange(m + 1), indexing='ij')
                msk = (rr < n) | (ss < m)
                rr, ss = rr[msk], ss[msk]
                v = np.einsum('kij,jk->i', G[n - rr, m - ss], Uh[:, ss, rr])
                v += tau2 * np.einsum('kij,jk->i', J[n - rr, m - ss], Wh[:, ss, rr])
                R -= v.reshape(Nx, Ny)
            U[:, :, m, n], W[:, :, m, n] = AInverse_3d(zeta[:, :, m, n], R, P.gammap_u, P.gammap_w, tau2)
            Uh[:, m, n] = fft2(U[:, :, m, n]).ravel() / K
            Wh[:, m, n] = fft2(W[:, :, m, n]).ravel() / K
    ubar = np.zeros(shp, dtype=complex); wbar = np.zeros(shp, dtype=complex)
    for n in range(N + 1):
        for m in range(M + 1):
            rr, ss = np.meshgrid(np.arange(n + 1), np.arange(m + 1), indexing='ij')
            rr, ss = rr.ravel(), ss.ravel()
            ubar[:, :, m, n] = np.einsum('kij,jk->i', Tu[n - rr, m - ss], Uh[:, ss, rr]).reshape(Nx, Ny)
            wbar[:, :, m, n] = np.einsum('kij,jk->i', Tw[n - rr, m - ss], Wh[:, ss, rr]).reshape(Nx, Ny)
    return U, W, ubar, wbar


def two_layer_solve_3d_auto(P, zeta, psi, verbose=False, mem_limit=1.5e9):
    K = P.Nx * P.Ny
    mem = 4 * (P.N + 1) * (P.M + 1) * K * K * 16 + 3 * (P.M + 1) * K * (P.Nz + 1) * K * 16 * 4
    if 4 * K <= 2 * P.N * P.M and mem < mem_limit:
        return two_layer_solve_3d_operator(P, zeta, psi, verbose)
    return two_layer_solve_3d(P, zeta, psi, verbose)


SOLVERS_3D = {'coupled': two_layer_solve_3d_coupled, 'lean': two_layer_solve_3d, 'operator': two_layer_solve_3d_operator,
              'auto': two_layer_solve_3d_auto}
