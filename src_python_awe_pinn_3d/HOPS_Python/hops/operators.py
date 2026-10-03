"""Fast "operator" formulation of the HOPS/AWE two-layer solver (new, not in MATLAB).

Why it is faster
----------------
two_layer_solve_fast.m calls the field solver once for every order (q,s) with data
U_{q,s} placed at order (0,0), in total  sum_{q,s} (N-q+1)(M-s+1) ~ (NM)^2/4
order-solves per layer (18,496 for N=M=15).  But the map  U -> G_{p,r}[U]  is LINEAR
and the same for every (q,s).  So we compute, once per layer, the Nx x Nx matrices

    G_{p,r}  (upper DNO),  J_{p,r}  (lower DNO),
    T^u_{p,r}: U -> ubar_{p,r}  (trace at z=a),   T^w_{p,r}: W -> wbar_{p,r}  (trace at z=-b)

by running the TFE recursion ONCE with all Nx unit vectors as a batch of right-hand
sides ((N+1)(M+1) = 256 order-solves, each a batched matrix-matrix product).  The
two-layer recursion then only needs small matrix-vector products:

    R_{n,m} = -psi_{n,m} - sum_{(r,s)<(n,m)} ( G_{n-r,m-s} U_{r,s} + tau2 J_{n-r,m-s} W_{r,s} )
    ubar_{n,m} = sum_{r<=n, s<=m} T^u_{n-r,m-s} U_{r,s}

This is algebraically identical to two_layer_solve_fast (agreement ~1e-12, i.e.
round-off) but does ~2x fewer flops for Nx=32, N=M=15 and removes almost all Python
overhead (everything is batched).  Flop ratio operator/fast ~ 4*Nx/(N*M): for very
large Nx (e.g. Nx=256 with N=M=20) the classical routine can be cheaper; use
``two_layer_solve_auto`` which picks for you.
"""
import numpy as np
from .expansions import T_dno
from . import config
from .bvp import bvp_matrix, _inverse

_fft = lambda u: np.fft.fft(u, axis=0)
_ifft = lambda u: np.fft.ifft(u, axis=0)


# Dense Nx x Nx x-derivative matrices are ~25% faster for Nx = 32 but lose up to ~0.7
# digits in round-off-dominated regimes (e.g. paper Fig. 5d, n_w = 10.1, N = M = 16),
# so FFTs are used by default (as in MATLAB).  Set to e.g. 128 to trade accuracy for speed.
DENSE_X_MAX = 0


_SETUP_CACHE = {}


def _layer_setup(upper, f, p, gammap, alpha, gamma, Dz, L, Nx, Nz, M, alphap, dense_x=None):
    """Everything that does not depend on the data: coefficients, the inverted
    collocation matrices, the DNO symbol expansion T_dno, the x-derivative.
    Cached, so repeated calls (two_layer_solve_lean) pay for it once."""
    import hashlib
    key = hashlib.sha1(b''.join(np.ascontiguousarray(np.asarray(v, dtype=complex)).tobytes()
                                for v in (f, p, gammap, alphap, [alpha, gamma, L, Nx, Nz, M, upper, config.ALPHA_FIX,
                                                                  -1 if dense_x is None else dense_x],
                                          Dz))).hexdigest()
    if key in _SETUP_CACHE:
        return _SETUP_CACHE[key]
    f = np.asarray(f, dtype=float)
    k2 = alphap[0] ** 2 + gammap[0] ** 2
    fx_r = np.real(_ifft(1j * p * _fft(f)))       # f_x used in the field solvers (real)
    fx_c = _ifft(1j * p * _fft(f))                # f_x used in the DNOs (complex, as MATLAB)
    ll = np.arange(Nz + 1)
    tz = np.cos(np.pi * ll / Nz)
    if upper:
        z_min, z_max = 0.0, L
    else:
        z_min, z_max = -L, 0.0
    z = ((z_max - z_min) / 2.0) * (tz - 1.0) + z_max
    D = (2.0 / (z_max - z_min)) * Dz
    D2 = D @ D
    DzL = (2.0 / L) * Dz                          # the 'dz' operator of dz.m

    F = f[:, None]
    FX = fx_r[:, None]
    if upper:
        h = (L - z)[None, :]
        A1_xx = -(2.0 / L) * F + 0 * h
        A1_xz = -(1.0 / L) * h * FX
        A2_xx = (1.0 / L ** 2) * F ** 2 + 0 * h
        A2_xz = (1.0 / L ** 2) * h * (F * FX)
        A2_zz = (1.0 / L ** 2) * h ** 2 * FX ** 2
        B1_x = (1.0 / L) * FX + 0 * h
        B2_x = -(1.0 / L ** 2) * F * FX + 0 * h
        B2_z = -(1.0 / L ** 2) * h * FX ** 2
        S1 = -(2.0 / L) * F + 0 * h
        S2 = (1.0 / L ** 2) * F ** 2 + 0 * h
        ell_bc = 0              # row holding the DNO (Robin) condition / trace ubar
        ell_if = Nz             # interface row
        A = bvp_matrix(1.0, 0.0, k2 - alphap ** 2, Nx, np.eye(Nz + 1), D, D2, D[0, :], D[-1, :],
                       top_n=1.0, top_d=-1j * gammap, bot_n=0.0, bot_d=1.0)
        sgn_T = +1.0
        phase = np.exp(1j * np.outer(gammap, z))
    else:
        h = (L + z)[None, :]
        A1_xx = (2.0 / L) * F + 0 * h
        A1_xz = -(1.0 / L) * h * FX
        A2_xx = (1.0 / L ** 2) * F ** 2 + 0 * h
        A2_xz = -(1.0 / L ** 2) * h * (F * FX)
        A2_zz = (1.0 / L ** 2) * h ** 2 * FX ** 2
        B1_x = -(1.0 / L) * FX + 0 * h
        B2_x = -(1.0 / L ** 2) * F * FX + 0 * h
        B2_z = (1.0 / L ** 2) * h * FX ** 2
        S1 = (2.0 / L) * F + 0 * h
        S2 = (1.0 / L ** 2) * F ** 2 + 0 * h
        ell_bc = Nz
        ell_if = 0
        A = bvp_matrix(1.0, 0.0, k2 - alphap ** 2, Nx, np.eye(Nz + 1), D, D2, D[0, :], D[-1, :],
                       top_n=0.0, top_d=1.0, bot_n=1.0, bot_d=1j * gammap)
        sgn_T = -1.0
        phase = np.exp(-1j * np.outer(gammap, z))
    A1_zx = A1_xz
    A2_zx = A2_xz
    ex = lambda C: C[:, :, None]
    A1_xx, A1_xz, A1_zx, A2_xx, A2_xz, A2_zx, A2_zz, B1_x, B2_x, B2_z, S1, S2 = map(
        ex, (A1_xx, A1_xz, A1_zx, A2_xx, A2_xz, A2_zx, A2_zz, B1_x, B2_x, B2_z, S1, S2))
    Ainv = np.linalg.inv(A)
    del A
    T = T_dno(alpha, alphap if config.ALPHA_FIX else p, gamma, gammap, k2, Nx, M)
    ip = (1j * p)[:, None, None]
    ip2 = (1j * p)[:, None]
    bdz = lambda u: np.matmul(DzL, u)
    if dense_x is None:
        dense_x_ = Nx <= DENSE_X_MAX
    else:
        dense_x_ = dense_x
    if dense_x_:
        # For small Nx, x-derivatives and Fourier multipliers are applied as dense
        # Nx x Nx matrices (one BLAS call) instead of an FFT pair: same operator
        # (Dx = ifft(i p fft(I))), much less overhead.
        I = np.eye(Nx)
        Dx = _ifft(1j * p[:, None] * _fft(I))
        Tm = [_ifft(T[:, k, None] * _fft(I)) for k in range(T.shape[1])]
        mm = lambda A_, u: (A_ @ u.reshape(Nx, -1)).reshape(u.shape)
        bdx = lambda u: mm(Dx, u)
        tmul = lambda k, v: Tm[k] @ v
    else:
        bdx = lambda u: _ifft(ip * _fft(u))
        tmul = lambda k, v: _ifft(T[:, k, None] * _fft(v))
    g2 = gamma ** 2
    S = dict(upper=upper, f=f, fx_c=fx_c, Nx=Nx, Nz=Nz, L=L, alpha=alpha, g2=g2, phase=phase,
             A1_xx=A1_xx, A1_xz=A1_xz, A1_zx=A1_zx, A2_xx=A2_xx, A2_xz=A2_xz, A2_zx=A2_zx,
             A2_zz=A2_zz, B1_x=B1_x, B2_x=B2_x, B2_z=B2_z, S1=S1, S2=S2, Ainv=Ainv,
             ell_bc=ell_bc, ell_if=ell_if, sgn_T=sgn_T, DzL=DzL, ip2=ip2,
             bdx=bdx, bdz=bdz, tmul=tmul, M=M)
    if len(_SETUP_CACHE) >= 4:
        _SETUP_CACHE.clear()
    _SETUP_CACHE[key] = S
    return S


def _layer_run(S, xi, N, M):
    """Run the TFE recursion + DNO for data xi (Nx, B) placed at order (0,0).

    Returns (DNO (N+1, M+1, Nx, B), trace at the artificial boundary (N+1, M+1, Nx, B)).
    Memory: only 3 n-levels of the volume field are kept.
    """
    assert M <= S['M']
    (upper, f, fx_c, Nx, Nz, L, alpha, g2, phase, A1_xx, A1_xz, A1_zx, A2_xx, A2_xz, A2_zx, A2_zz,
     B1_x, B2_x, B2_z, S1, S2, Ainv, ell_bc, ell_if, sgn_T, DzL, ip2, bdx, bdz, tmul) = (
        S[k] for k in ('upper', 'f', 'fx_c', 'Nx', 'Nz', 'L', 'alpha', 'g2', 'phase', 'A1_xx', 'A1_xz',
                       'A1_zx', 'A2_xx', 'A2_xz', 'A2_zx', 'A2_zz', 'B1_x', 'B2_x', 'B2_z', 'S1', 'S2',
                       'Ainv', 'ell_bc', 'ell_if', 'sgn_T', 'DzL', 'ip2', 'bdx', 'bdz', 'tmul'))
    xi = np.asarray(xi, dtype=complex)
    if xi.ndim == 1:
        xi = xi[:, None]
    B = xi.shape[1]

    # storage: 3 rolling n-levels of volume fields; traces for all orders
    V = np.zeros((3, M + 1, Nx, Nz + 1, B), dtype=complex)
    tr_bc = np.zeros((N + 1, M + 1, Nx, B), dtype=complex)     # u at ell_bc (for T/J terms, ubar)
    tr_if = np.zeros((N + 1, M + 1, Nx, B), dtype=complex)     # u at interface
    tz_if = np.zeros((N + 1, M + 1, Nx, B), dtype=complex)     # u_z at interface
    xi_hat = _fft(xi)                                          # data at order (0,0)
    get = lambda m, n: V[n % 3, m]

    for n in range(N + 1):
        V[n % 3] = 0.0
        for m in range(M + 1):
            if n == 0 and m == 0:
                u = _ifft(phase[:, :, None] * xi_hat[:, None, :])
            else:
                Fn = np.zeros((Nx, Nz + 1, B), dtype=complex)
                a0 = alpha != 0          # skip the 2 i alpha d_x terms when alpha = 0
                afix = config.ALPHA_FIX  # + 2 i alpha A_xz d_z (see hops/config.py, item 4)
                if n >= 1:
                    u1 = get(m, n - 1)
                    ux = bdx(u1)
                    uz = bdz(u1)
                    Fn -= bdx(A1_xx * ux + A1_xz * uz) + bdz(A1_zx * ux) + B1_x * ux + g2 * S1 * u1
                    if a0:
                        Fn -= 2j * alpha * S1 * ux
                        if afix:
                            Fn -= 2j * alpha * A1_xz * uz
                if m >= 1:
                    u1 = get(m - 1, n)
                    Fn -= 2 * g2 * u1
                    if a0:
                        Fn -= 2j * alpha * bdx(u1)
                if n >= 1 and m >= 1:
                    u1 = get(m - 1, n - 1)
                    Fn -= 2 * g2 * S1 * u1
                    if a0:
                        Fn -= 2j * alpha * S1 * bdx(u1)
                        if afix:
                            Fn -= 2j * alpha * A1_xz * bdz(u1)
                if n >= 2:
                    u1 = get(m, n - 2)
                    ux = bdx(u1)
                    uz = bdz(u1)
                    Fn -= (bdx(A2_xx * ux + A2_xz * uz) + bdz(A2_zx * ux + A2_zz * uz)
                           + B2_x * ux + B2_z * uz + g2 * S2 * u1)
                    if a0:
                        Fn -= 2j * alpha * S2 * ux
                        if afix:
                            Fn -= 2j * alpha * A2_xz * uz
                if m >= 2:
                    Fn -= g2 * get(m - 2, n)
                if n >= 1 and m >= 2:
                    Fn -= g2 * S1 * get(m - 2, n - 1)
                if n >= 2 and m >= 1:
                    u1 = get(m - 1, n - 2)
                    Fn -= 2 * g2 * S2 * u1
                    if a0:
                        Fn -= 2j * alpha * S2 * bdx(u1)
                        if afix:
                            Fn -= 2j * alpha * A2_xz * bdz(u1)
                if n >= 2 and m >= 2:
                    Fn -= g2 * S2 * get(m - 2, n - 2)
                # boundary data at the artificial boundary (J_{n,m} upper / Q_{n,m} lower)
                Jn = np.zeros((Nx, B), dtype=complex)
                for r in range(m):
                    Jn += sgn_T * tmul(m - r, tr_bc[n, r])
                if n >= 1:
                    for r in range(m + 1):
                        Jn -= (1.0 / L) * f[:, None] * tmul(m - r, tr_bc[n - 1, r])
                Fh = _fft(Fn)
                Fh[:, ell_bc, :] = _fft(Jn)
                Fh[:, ell_if, :] = 0.0          # Dirichlet data xi_{n,m} = 0 for (n,m) != (0,0)
                u = _ifft(np.matmul(Ainv, Fh))
            V[n % 3, m] = u
            tr_bc[n, m] = u[:, ell_bc, :]
            tr_if[n, m] = u[:, ell_if, :]
            tz_if[n, m] = np.matmul(DzL[ell_if], u)          # (dz u) at interface row

    # DNO operators from the stored traces (dno_tfe_helmholtz_m_and_n[_lf].m)
    fxc = fx_c[:, None]
    fc = f[:, None]
    ddx = lambda v: _ifft(ip2 * _fft(v))
    Gop = np.zeros((N + 1, M + 1, Nx, B), dtype=complex)
    for n in range(N + 1):
        for m in range(M + 1):
            if upper:
                G = -tz_if[n, m]
                if n >= 1:
                    G = G + fxc * ddx(tr_if[n - 1, m]) + (1.0 / L) * (fc * Gop[n - 1, m])
                if n >= 2:
                    G = G - (1.0 / L) * (fc * (fxc * ddx(tr_if[n - 2, m]))) - fxc * (fxc * tz_if[n - 2, m])
            else:
                G = tz_if[n, m].copy()
                if n >= 1:
                    G = G - fxc * ddx(tr_if[n - 1, m]) - (1.0 / L) * (fc * Gop[n - 1, m])
                if n >= 2:
                    G = G - (1.0 / L) * (fc * (fxc * ddx(tr_if[n - 2, m]))) + fxc * (fxc * tz_if[n - 2, m])
            Gop[n, m] = G
    return Gop, tr_bc


def _layer_operators(upper, f, p, gammap, alpha, gamma, Dz, L, Nx, Nz, N, M, alphap, dense_x=None):
    """Return (DNO_ops, TRACE_ops), each (N+1, M+1, Nx, Nx) complex (data = all unit vectors).

    upper=True : field_tfe_helmholtz_m_and_n + dno_tfe_helmholtz_m_and_n  (z in [0,a])
    upper=False: field_tfe_helmholtz_m_and_n_lf + dno_tfe_helmholtz_m_and_n_lf (z in [-b,0])
    """
    S = _layer_setup(upper, f, p, gammap, alpha, gamma, Dz, L, Nx, Nz, M, alphap, dense_x)
    # Input basis = Fourier modes e^{i p_k x} (column k), NOT physical unit vectors:
    # round-off in each column is then relative to that mode's own response, so applying
    # the operator to smooth data (small high-mode coefficients) is as accurate as the
    # classic field solve.  (With a unit-vector basis the large responses of the highly
    # evanescent modes pollute every entry; for n_w = 10.1, N = M = 16 that costs ~1.5 digits.)
    modes = np.fft.ifft(np.eye(Nx) * Nx, axis=0)          # columns e^{i p_k x_j}
    return _layer_run(S, modes, N, M)


def layer_operators(f, pp, gamma_u_bar_p, gamma_w_bar_p, alpha_bar, gamma_u_bar, gamma_w_bar,
                    Dz, a, b, Nx, Nz, N, M, alpha_bar_p):
    """(G_ops, Tu_ops, J_ops, Tw_ops), each of shape (N+1, M+1, Nx, Nx)."""
    Gop, Tu = _layer_operators(True, f, pp, gamma_u_bar_p, alpha_bar, gamma_u_bar, Dz, a,
                               Nx, Nz, N, M, alpha_bar_p)
    Jop, Tw = _layer_operators(False, f, pp, gamma_w_bar_p, alpha_bar, gamma_w_bar, Dz, b,
                               Nx, Nz, N, M, alpha_bar_p)
    return Gop, Tu, Jop, Tw


def two_layer_solve_operator(tau2, zeta_n_m, psi_n_m, gamma_u_bar_p, gamma_w_bar_p, N, Nx, f, f_x,
                             pp, alpha_bar, gamma_u_bar, gamma_w_bar, Dz, a, b, Nz, M, identy,
                             alpha_bar_p, ops=None):
    """Drop-in replacement for two_layer_solve_fast (same arguments and outputs)."""
    from .two_layer import AInverse, phase_correction
    if ops is None:
        ops = layer_operators(f, pp, gamma_u_bar_p, gamma_w_bar_p, alpha_bar, gamma_u_bar,
                              gamma_w_bar, Dz, a, b, Nx, Nz, N, M, alpha_bar_p)
    Gop, Tu, Jop, Tw = ops              # act on Fourier coefficients  c = fft(U)/Nx
    U = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    W = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    Uh = np.zeros_like(U)
    Wh = np.zeros_like(W)
    for n in range(N + 1):
        for m in range(M + 1):
            R = -psi_n_m[:, m, n].astype(complex) - phase_correction(U, W, m, n, alpha_bar, f_x, tau2)
            if n > 0 or m > 0:
                # sum over (r,s) < (n,m): G_{n-r,m-s} U_{r,s} + tau2 J_{n-r,m-s} W_{r,s}
                rr, ss = np.meshgrid(np.arange(n + 1), np.arange(m + 1), indexing='ij')
                mask = (rr < n) | (ss < m)
                rr, ss = rr[mask], ss[mask]
                R -= np.einsum('kij,jk->i', Gop[n - rr, m - ss], Uh[:, ss, rr])
                R -= tau2 * np.einsum('kij,jk->i', Jop[n - rr, m - ss], Wh[:, ss, rr])
            U[:, m, n], W[:, m, n] = AInverse(zeta_n_m[:, m, n], R, gamma_u_bar_p, gamma_w_bar_p,
                                              Nx, tau2)
            Uh[:, m, n] = np.fft.fft(U[:, m, n]) / Nx
            Wh[:, m, n] = np.fft.fft(W[:, m, n]) / Nx
    ubar = np.zeros_like(U)
    wbar = np.zeros_like(W)
    for n in range(N + 1):
        for m in range(M + 1):
            rr, ss = np.meshgrid(np.arange(n + 1), np.arange(m + 1), indexing='ij')
            rr, ss = rr.ravel(), ss.ravel()
            ubar[:, m, n] = np.einsum('kij,jk->i', Tu[n - rr, m - ss], Uh[:, ss, rr])
            wbar[:, m, n] = np.einsum('kij,jk->i', Tw[n - rr, m - ss], Wh[:, ss, rr])
    return U, W, ubar, wbar


def two_layer_solve_lean(tau2, zeta_n_m, psi_n_m, gamma_u_bar_p, gamma_w_bar_p, N, Nx, f, f_x,
                         pp, alpha_bar, gamma_u_bar, gamma_w_bar, Dz, a, b, Nz, M, identy,
                         alpha_bar_p, verbose=False, threads=True):
    """Low-memory version of two_layer_solve_fast for large Nx (e.g. Nx = 1024, Nz = 128,
    N = M = 20, paper Figs. 6-8).  Same arithmetic as the classic routine, but

    * instead of storing G_{p,r}[U_{q,s}] in an (Nx, M+1, M+1, N+1, N+1) array (3.2 GB
      for Nx = 1024, N = M = 20) the contributions are added to an (Nx, M+1, N+1)
      accumulator as soon as U_{q,s} is known;
    * ubar/wbar are accumulated the same way (no extra full field solve);
    * the field recursion keeps only 3 n-levels of the volume field (130 MB instead of
      0.9 GB), and the collocation matrices are inverted once per layer.
    """
    from .two_layer import AInverse, phase_correction
    import time
    Su = _layer_setup(True, f, pp, gamma_u_bar_p, alpha_bar, gamma_u_bar, Dz, a, Nx, Nz, M, alpha_bar_p)
    Sw = _layer_setup(False, f, pp, gamma_w_bar_p, alpha_bar, gamma_w_bar, Dz, b, Nx, Nz, M, alpha_bar_p)
    U = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    W = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    acc = np.zeros((Nx, M + 1, N + 1), dtype=complex)     # sum_{(r,s)<(n,m)} G[U] + tau2 J[W]
    ubar = np.zeros_like(U)
    wbar = np.zeros_like(W)
    from concurrent.futures import ThreadPoolExecutor
    pool = ThreadPoolExecutor(1) if threads else None
    t0 = time.time()
    for n in range(N + 1):
        for m in range(M + 1):
            R = -psi_n_m[:, m, n] - acc[:, m, n] - phase_correction(U, W, m, n, alpha_bar, f_x, tau2)
            U[:, m, n], W[:, m, n] = AInverse(zeta_n_m[:, m, n], R, gamma_u_bar_p, gamma_w_bar_p,
                                              Nx, tau2)
            if pool is not None:      # upper and lower layers are independent -> 2 threads
                fu = pool.submit(_layer_run, Su, U[:, m, n], N - n, M - m)
                J, Tw = _layer_run(Sw, W[:, m, n], N - n, M - m)
                G, Tu = fu.result()
            else:
                G, Tu = _layer_run(Su, U[:, m, n], N - n, M - m)
                J, Tw = _layer_run(Sw, W[:, m, n], N - n, M - m)
            GJ = (G + tau2 * J)[..., 0]                   # (N-n+1, M-m+1, Nx)
            GJ[0, 0] = 0.0                                # (r,s) = (n,m) itself is excluded
            acc[:, m:, n:] += np.transpose(GJ, (2, 1, 0))
            ubar[:, m:, n:] += np.transpose(Tu[..., 0], (2, 1, 0))
            wbar[:, m:, n:] += np.transpose(Tw[..., 0], (2, 1, 0))
        if verbose:
            print(f'  two_layer_solve_lean: n = {n}/{N} done ({time.time() - t0:.0f} s)', flush=True)
    if pool is not None:
        pool.shutdown()
    return U, W, ubar, wbar


def two_layer_solve_auto(*args, **kw):
    """Pick the cheapest formulation: the coupled solver (hops/coupled.py) is fastest in every regime
    tested (Nx = 32 ... 1024) and needs no Nx x Nx operators."""
    from .coupled import two_layer_solve_coupled
    return two_layer_solve_coupled(*args)
    N, Nx, M, Nz = args[5], args[6], args[17], args[16]
    mem = 3 * (M + 1) * Nx * (Nz + 1) * Nx * 16 + 2 * (N + 1) * (M + 1) * Nx * Nx * 16 * 3
    if 4 * Nx <= 2 * N * M and mem < 2e9:
        return two_layer_solve_operator(*args, **kw)
    return two_layer_solve_lean(*args, **kw)
