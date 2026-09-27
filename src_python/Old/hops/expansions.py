"""Frequency/boundary expansions: gamma_exp.m, E_exp.m, E_exp_lf.m, A_exp.m, T_dno.m."""
import numpy as np


def gamma_exp(alpha_bar, alpha_bar_q, gamma_bar, gamma_bar_q, k_bar, M):
    """Taylor coefficients (in delta) of gamma_q(delta); returns length M+1."""
    g = np.zeros(M + 1, dtype=complex)
    g[0] = gamma_bar_q
    if M >= 1:
        g[1] = (k_bar ** 2 - alpha_bar * alpha_bar_q) / g[0]
    if M >= 2:
        g[2] = (gamma_bar ** 2 - g[1] ** 2) / (2 * g[0])
    for m in range(3, M + 1):
        num_sum = 0.0
        for r in range(1, m):
            num_sum = num_sum - g[m - r] * g[r]
        g[m] = num_sum / (2 * g[0])
    return g


def _E(gamma_m, f, N, M, sign):
    f = np.asarray(f)
    Nx = f.size
    E = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    E[:, 0, 0] = 1.0
    for n in range(1, N + 1):
        for m in range(M + 1):
            s = np.zeros(Nx, dtype=complex)
            for r in range(m + 1):
                s = s + E[:, r, n - 1] * (1j * gamma_m[m - r])
            E[:, m, n] = sign * f * s / n
    return E


def E_exp(gamma_m, f, N, M):
    """Expansion of exp(i*gamma*eps*f) -> E_nm, shape (Nx, M+1, N+1)."""
    return _E(gamma_m, f, N, M, +1.0)


def E_exp_lf(gamma_m, f, N, M):
    """Expansion of exp(-i*gamma*eps*f) -> E_nm, shape (Nx, M+1, N+1)."""
    return _E(gamma_m, f, N, M, -1.0)


def A_exp(alpha_m, f_x, N, M):
    """A_exp.m (note: MATLAB layout here is (Nx, N+1, M+1))."""
    f_x = np.asarray(f_x)
    Nx = f_x.size
    A = np.zeros((Nx, N + 1, M + 1), dtype=complex)
    A[:, 0, 0] = 1.0
    for n in range(1, N + 1):
        for m in range(M + 1):
            s = np.zeros(Nx, dtype=complex)
            for r in range(m + 1):
                s = s + A[:, n - 1, r] * (1j * alpha_m[m - r])
            A[:, n, m] = f_x * s / n
    return A


def T_dno(alpha, alphap, gamma, gammap, k2, Nx, M):
    """Frequency expansion of the DNO symbol i*gamma_p(delta); shape (Nx, M+1).

    Paper: HOPS/AWE, Sec. 4.1, eqs. (12)-(16).  Requires M >= 2 (as in MATLAB).
    """
    g = np.zeros((Nx, max(M, 2) + 1), dtype=complex)
    g[:, 0] = gammap
    g[:, 1] = 2 * (k2 - alpha * alphap) / (2 * gammap)
    g[:, 2] = (gamma ** 2 - g[:, 1] ** 2) / (2 * gammap)
    for m in range(3, M + 1):
        for r in range(1, m):
            g[:, m] = g[:, m] - (g[:, m - r] * g[:, r]) / (2 * gammap)
    return 1j * g
