"""Coupled HOPS/AWE two-layer solver for 2D (new; same algorithm as hops3d.two_layer_solve_3d_coupled).

Idea
----
two_layer_solve_fast.m needs, at order (n, m),

    S_{n,m} = sum_{(r,s) < (n,m)}  G_{n-r,m-s}[U_{r,s}]  (+ tau2 J[W]),

and gets it by running a separate TFE recursion for every U_{r,s} (~(N+1)^2 (M+1)^2 / 4 order-solves
per layer).  But the recursion is LINEAR and SHIFT-INVARIANT in (n, m): the field of the full data
U = sum U_{r,s} eps^r delta^s is  u = sum_{r,s} (field of U_{r,s} shifted by (r,s)).  So run ONE
recursion per layer on the full data, interleaved with the interface solve:

  at order (n, m):
    1. particular field u^p_{n,m}: RHS and artificial-boundary data from lower orders, zero Dirichlet data
    2. its DNO (with the lower-order DNO terms)  ==  S_{n,m}        (exactly the sum above)
    3. AInverse -> U_{n,m}, W_{n,m}
    4. u_{n,m} = u^p_{n,m} + e^{+-i gamma_p z} U^_{n,m}  (the exact homogeneous solution, which is what
       two_layer_solve_fast uses for data placed at order (0,0)); store traces, full DNO.

Cost: (N+1)(M+1) order-solves per layer instead of ~(N+1)^2(M+1)^2/4  (N = M = 15: 256 vs 18,496),
with no Nx x Nx operators (unlike two_layer_solve_operator), so it also scales to Nx = 1024.
Algebraically identical to two_layer_solve_fast / two_layer_solve_operator (agreement ~ round-off).
Uses the same data-independent setup as hops/operators.py (_layer_setup; honours config.ALPHA_FIX).
"""
import numpy as np
from . import config
from .operators import _layer_setup, _fft, _ifft


class _Layer:
    """Order-by-order TFE recursion for one layer, split into partial() / complete()."""

    def __init__(self, S, N, M):
        self.S, self.N, self.M = S, N, M
        Nx, Nz = S['Nx'], S['Nz']
        self.V = np.zeros((3, M + 1, Nx, Nz + 1, 1), dtype=complex)
        self.tr_bc = np.zeros((N + 1, M + 1, Nx, 1), dtype=complex)
        self.tr_if = np.zeros_like(self.tr_bc)
        self.tz_if = np.zeros_like(self.tr_bc)
        self.D = np.zeros_like(self.tr_bc)
        self.uhp = None
        self.cleared = -1
        ip2 = S['ip2']
        self.ddx = lambda v: _ifft(ip2 * _fft(v))

    def _rhs(self, n, m):
        S = self.S
        (L, alpha, g2, A1_xx, A1_xz, A1_zx, A2_xx, A2_xz, A2_zx, A2_zz, B1_x, B2_x, B2_z, S1, S2, bdx, bdz) = (
            S[k] for k in ('L', 'alpha', 'g2', 'A1_xx', 'A1_xz', 'A1_zx', 'A2_xx', 'A2_xz', 'A2_zx', 'A2_zz',
                           'B1_x', 'B2_x', 'B2_z', 'S1', 'S2', 'bdx', 'bdz'))
        get = lambda mm, nn: self.V[nn % 3, mm]
        Fn = np.zeros(self.V.shape[2:], dtype=complex)
        a0 = alpha != 0
        afix = config.ALPHA_FIX
        if n >= 1:
            u1 = get(m, n - 1); ux = bdx(u1); uz = bdz(u1)
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
            u1 = get(m, n - 2); ux = bdx(u1); uz = bdz(u1)
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
        return Fn

    def _dno(self, n, m, tz):
        S = self.S
        L, fxc, fc = S['L'], S['fx_c'][:, None], S['f'][:, None]
        tr_if, tz_if, D, ddx = self.tr_if, self.tz_if, self.D, self.ddx
        if S['upper']:
            G = -tz
            if n >= 1:
                G = G + fxc * ddx(tr_if[n - 1, m]) + (1.0 / L) * (fc * D[n - 1, m])
            if n >= 2:
                G = G - (1.0 / L) * (fc * (fxc * ddx(tr_if[n - 2, m]))) - fxc * (fxc * tz_if[n - 2, m])
        else:
            G = tz.copy()
            if n >= 1:
                G = G - fxc * ddx(tr_if[n - 1, m]) - (1.0 / L) * (fc * D[n - 1, m])
            if n >= 2:
                G = G - (1.0 / L) * (fc * (fxc * ddx(tr_if[n - 2, m]))) + fxc * (fxc * tz_if[n - 2, m])
        return G

    def partial(self, n, m):
        """Particular field of order (n, m); returns its DNO = sum over lower orders."""
        S = self.S
        if self.cleared != n:
            self.V[n % 3] = 0.0
            self.cleared = n
        if n == 0 and m == 0:
            self.uhp = None
            return np.zeros((S['Nx'], 1), dtype=complex)
        Fn = self._rhs(n, m)
        Jn = np.zeros((S['Nx'], 1), dtype=complex)
        for r in range(m):
            Jn += S['sgn_T'] * S['tmul'](m - r, self.tr_bc[n, r])
        if n >= 1:
            for r in range(m + 1):
                Jn -= (1.0 / S['L']) * S['f'][:, None] * S['tmul'](m - r, self.tr_bc[n - 1, r])
        Fh = _fft(Fn)
        Fh[:, S['ell_bc'], :] = _fft(Jn)
        Fh[:, S['ell_if'], :] = 0.0
        self.uhp = np.matmul(S['Ainv'], Fh)                    # Fourier coefficients (Nx, Nz+1, 1)
        tz = _ifft(np.matmul(S['DzL'][S['ell_if']], self.uhp))
        return self._dno(n, m, tz)

    def complete(self, n, m, xi):
        """Add the homogeneous part e^{+-i gamma z} xi^; store; return (DNO_{n,m}, trace at the artificial bdry)."""
        S = self.S
        uh = S['phase'][:, :, None] * _fft(xi)[:, None, None]
        if self.uhp is not None:
            uh = uh + self.uhp
        u = _ifft(uh)
        self.V[n % 3, m] = u
        self.tr_bc[n, m] = u[:, S['ell_bc'], :]
        self.tr_if[n, m] = u[:, S['ell_if'], :]
        self.tz_if[n, m] = np.matmul(S['DzL'][S['ell_if']], u)
        self.D[n, m] = self._dno(n, m, self.tz_if[n, m])
        self.uhp = None
        return self.D[n, m], self.tr_bc[n, m]


def two_layer_solve_coupled(tau2, zeta_n_m, psi_n_m, gamma_u_bar_p, gamma_w_bar_p, N, Nx, f, f_x, pp,
                            alpha_bar, gamma_u_bar, gamma_w_bar, Dz, a, b, Nz, M, identy, alpha_bar_p,
                            verbose=False):
    """Drop-in replacement for two_layer_solve_fast (same arguments, same outputs U, W, ubar, wbar)."""
    from .two_layer import AInverse, phase_correction
    Su = _layer_setup(True, f, pp, gamma_u_bar_p, alpha_bar, gamma_u_bar, Dz, a, Nx, Nz, M, alpha_bar_p)
    Sw = _layer_setup(False, f, pp, gamma_w_bar_p, alpha_bar, gamma_w_bar, Dz, b, Nx, Nz, M, alpha_bar_p)
    lu, lw = _Layer(Su, N, M), _Layer(Sw, N, M)
    U = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    W = np.zeros_like(U)
    ubar = np.zeros_like(U)
    wbar = np.zeros_like(U)
    for n in range(N + 1):
        for m in range(M + 1):
            Gp = lu.partial(n, m)[:, 0]
            Jp = lw.partial(n, m)[:, 0]
            R = -psi_n_m[:, m, n] - Gp - tau2 * Jp - phase_correction(U, W, m, n, alpha_bar, f_x, tau2)
            U[:, m, n], W[:, m, n] = AInverse(zeta_n_m[:, m, n], R, gamma_u_bar_p, gamma_w_bar_p, Nx, tau2)
            ubar[:, m, n] = lu.complete(n, m, U[:, m, n])[1][:, 0]
            wbar[:, m, n] = lw.complete(n, m, W[:, m, n])[1][:, 0]
        if verbose:
            print(f'  two_layer_solve_coupled: n = {n}/{N}', flush=True)
    return U, W, ubar, wbar
