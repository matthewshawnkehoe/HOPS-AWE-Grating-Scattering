"""Frequency / boundary expansions in 3D (analogues of gamma_exp.m, E_exp.m, E_exp_lf.m, T_dno.m).

Frequency perturbation (paper Sec. 3, extended to two lateral directions):
    omega = omega_bar (1 + delta),  k = k_bar (1 + delta),
    alpha = alpha_bar (1 + delta),  beta = beta_bar (1 + delta),
so that for every Fourier mode (p, q)
    gamma_pq(delta)^2 = k^2 - (alpha + p)^2 - (beta + q)^2
                      = gamma_pq_bar^2 + 2 delta (k_bar^2 - alpha_bar alpha_p - beta_bar beta_q)
                        + delta^2 gamma_bar^2,           gamma_bar^2 = k_bar^2 - alpha_bar^2 - beta_bar^2.
Writing gamma_pq(delta) = sum_m g_m delta^m and matching powers gives
    g_0 = gamma_pq_bar,  g_1 = (k_bar^2 - alpha_bar alpha_p - beta_bar beta_q) / g_0,
    g_2 = (gamma_bar^2 - g_1^2) / (2 g_0),  g_m = - sum_{r=1}^{m-1} g_{m-r} g_r / (2 g_0)  (m >= 3),
i.e. exactly the 2D recursion (paper eqs. (12)-(16)) with alpha alpha_p -> alpha alpha_p + beta beta_q.
"""
import numpy as np


def gamma_exp_3d(alpha_bar, beta_bar, alpha_r, beta_s, gamma_bar, gamma_rs, k_bar, M):
    """Taylor coefficients (in delta) of gamma_rs(delta) for ONE mode; length max(M,2)+1."""
    g = np.zeros(max(M, 2) + 1, dtype=complex)
    g[0] = gamma_rs
    g[1] = (k_bar ** 2 - alpha_bar * alpha_r - beta_bar * beta_s) / g[0]
    g[2] = (gamma_bar ** 2 - g[1] ** 2) / (2 * g[0])
    for m in range(3, M + 1):
        g[m] = -sum(g[m - r] * g[r] for r in range(1, m)) / (2 * g[0])
    return g[:M + 1]


def T_dno_3d(alpha_bar, beta_bar, alphap, betap, gamma_bar, gammap, k2, M):
    """Expansion of the artificial-boundary DNO symbol i gamma_pq(delta); shape (Nx, Ny, max(M,2)+1)."""
    Nx, Ny = gammap.shape
    g = np.zeros((Nx, Ny, max(M, 2) + 1), dtype=complex)
    g[..., 0] = gammap
    g[..., 1] = (k2 - alpha_bar * alphap[:, None] - beta_bar * betap[None, :]) / gammap
    g[..., 2] = (gamma_bar ** 2 - g[..., 1] ** 2) / (2 * gammap)
    for m in range(3, M + 1):
        for r in range(1, m):
            g[..., m] -= g[..., m - r] * g[..., r] / (2 * gammap)
    return 1j * g


def E_exp_3d(gamma_m, f, N, M, sign=+1.0):
    """Joint (eps, delta) expansion of exp(sign * i gamma(delta) eps f(x,y)).

    gamma_m: Taylor coefficients of gamma(delta) (length >= M+1);  f: (Nx, Ny).
    Returns E of shape (Nx, Ny, M+1, N+1)  with  E_{n,m} = sign f/n sum_r E_{r,n-1} i gamma_{m-r}.
    """
    f = np.asarray(f)
    gm = np.zeros(M + 1, dtype=complex)
    gm[:min(M + 1, len(gamma_m))] = np.asarray(gamma_m)[:M + 1]
    E = np.zeros(f.shape + (M + 1, N + 1), dtype=complex)
    E[..., 0, 0] = 1.0
    for n in range(1, N + 1):
        for m in range(M + 1):
            s = np.zeros(f.shape, dtype=complex)
            for r in range(m + 1):
                s = s + E[..., r, n - 1] * (1j * gm[m - r])
            E[..., m, n] = sign * f * s / n
    return E
