"""Two-layer HOPS/AWE solver: AInverse.m, two_layer_solve.m, two_layer_solve_fast.m.

Solves, order by order in (eps, delta),
    U_{n,m} - W_{n,m}                   = zeta_{n,m}
    G_{0,0}[U_{n,m}] + tau2 J_{0,0}[W_{n,m}] = -psi_{n,m} - sum_{(r,s)<(n,m)} (G_{n-r,m-s}[U_{r,s}] + tau2 J_{n-r,m-s}[W_{r,s}])
Note the sign convention: the lower DNO J is +d_N w, so psi = -nu_u - tau2*nu_w
where nu_w already carries the '+d_N' sign (see setup_xi_w_nu_w_n_m).
"""
import numpy as np
from . import config
from .fields import (field_tfe_helmholtz_m_and_n, field_tfe_helmholtz_m_and_n_lf,
                     dno_tfe_helmholtz_m_and_n, dno_tfe_helmholtz_m_and_n_lf)


def phase_correction(U, W, m, n, alpha_bar, f_x, tau2):
    """i alpha(delta) g_x (U - tau^2 W) at order (n, m), alpha(delta) = alpha_bar (1 + delta), g = eps f
    (see hops/config.py, item 3).  Zero when alpha_bar = 0 or config.ALPHA_FIX is False."""
    if not config.ALPHA_FIX or alpha_bar == 0 or n < 1:
        return 0.0
    X = U[:, m, n - 1] - tau2 * W[:, m, n - 1]
    if m >= 1:
        X = X + U[:, m - 1, n - 1] - tau2 * W[:, m - 1, n - 1]
    return 1j * alpha_bar * f_x * X


def AInverse(Q, R, gammap, gammapw, Nx, tau2):
    """Invert the flat-interface 2x2 Fourier multiplier; returns (U, W)."""
    Q_hat = np.fft.fft(Q)
    R_hat = np.fft.fft(R)
    det_p = tau2 * gammapw + gammap
    a = ((tau2 * gammapw) * Q_hat + 1j * R_hat) / det_p
    b = ((-gammap) * Q_hat + 1j * R_hat) / det_p
    return np.fft.ifft(a), np.fft.ifft(b)


def _traces(U_n_m, W_n_m, f, pp, gamma_u_bar_p, gamma_w_bar_p, alpha_bar, gamma_u_bar,
            gamma_w_bar, Dz, a, b, Nx, Nz, N, M, identy, alpha_bar_p):
    u_n_m = field_tfe_helmholtz_m_and_n(U_n_m, f, pp, gamma_u_bar_p, alpha_bar, gamma_u_bar,
                                        Dz, a, Nx, Nz, N, M, identy, alpha_bar_p)
    w_n_m = field_tfe_helmholtz_m_and_n_lf(W_n_m, f, pp, gamma_w_bar_p, alpha_bar, gamma_w_bar,
                                           Dz, b, Nx, Nz, N, M, identy, alpha_bar_p)
    ubar = u_n_m[:, 0, :, :].copy()      # ell_top
    wbar = w_n_m[:, Nz, :, :].copy()     # ell_bottom
    return ubar, wbar


def two_layer_solve(tau2, zeta_n_m, psi_n_m, gamma_u_bar_p, gamma_w_bar_p, N, Nx, f, f_x, pp,
                    alpha_bar, gamma_u_bar, gamma_w_bar, Dz, a, b, Nz, M, identy, alpha_bar_p):
    """Reference (slow) version: recomputes G_{n-r,m-s}[U_{r,s}] every time."""
    U = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    W = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    U[:, 0, 0], W[:, 0, 0] = AInverse(zeta_n_m[:, 0, 0], -psi_n_m[:, 0, 0],
                                      gamma_u_bar_p, gamma_w_bar_p, Nx, tau2)
    for n in range(N + 1):
        for m in range(M + 1):
            Q = zeta_n_m[:, m, n]
            R = -psi_n_m[:, m, n]
            # MSK 11/01/2021 - phase terms (MATLAB 1-based U_n_m(:,m,n-1) -> [:, m-1, n-2])
            if config.ALPHA_FIX:
                R = R - phase_correction(U, W, m, n, alpha_bar, f_x, tau2)
            elif n > 1 and m > 0:
                R = R - f_x * 1j * alpha_bar * U[:, m - 1, n - 2] + tau2 * f_x * 1j * alpha_bar * W[:, m - 1, n - 2]
                if m > 1:
                    R = R - f_x * 1j * alpha_bar * U[:, m - 2, n - 2] + tau2 * f_x * 1j * alpha_bar * W[:, m - 2, n - 2]
            for r in range(n + 1):
                for s in range(m + 1):
                    if r < n or s < m:
                        xi = np.zeros((Nx, M + 1, N + 1), dtype=complex)
                        xi[:, 0, 0] = U[:, s, r]
                        u = field_tfe_helmholtz_m_and_n(xi, f, pp, gamma_u_bar_p, alpha_bar, gamma_u_bar,
                                                        Dz, a, Nx, Nz, n - r, m - s, identy, alpha_bar_p)
                        G = dno_tfe_helmholtz_m_and_n(u, f, pp, Dz, a, Nx, Nz, n - r, m - s)
                        R = R - G[:, m - s, n - r]
                        xi = np.zeros((Nx, M + 1, N + 1), dtype=complex)
                        xi[:, 0, 0] = W[:, s, r]
                        w = field_tfe_helmholtz_m_and_n_lf(xi, f, pp, gamma_w_bar_p, alpha_bar, gamma_w_bar,
                                                           Dz, b, Nx, Nz, n - r, m - s, identy, alpha_bar_p)
                        J = dno_tfe_helmholtz_m_and_n_lf(w, f, pp, Dz, b, Nx, Nz, n - r, m - s)
                        R = R - tau2 * J[:, m - s, n - r]
            if n > 0 or m > 0:
                U[:, m, n], W[:, m, n] = AInverse(Q, R, gamma_u_bar_p, gamma_w_bar_p, Nx, tau2)
    ubar, wbar = _traces(U, W, f, pp, gamma_u_bar_p, gamma_w_bar_p, alpha_bar, gamma_u_bar,
                         gamma_w_bar, Dz, a, b, Nx, Nz, N, M, identy, alpha_bar_p)
    return U, W, ubar, wbar


def two_layer_solve_fast(tau2, zeta_n_m, psi_n_m, gamma_u_bar_p, gamma_w_bar_p, N, Nx, f, f_x, pp,
                         alpha_bar, gamma_u_bar, gamma_w_bar, Dz, a, b, Nz, M, identy, alpha_bar_p,
                         verbose=False):
    """Fast version: stores G_{p,r}[U_{q,s}] and J_{p,r}[W_{q,s}] once per (q,s).

    Returns U_n_m, W_n_m, ubar_n_m, wbar_n_m, each of shape (Nx, M+1, N+1).
    """
    U = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    W = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    # G_U[:, r, s, p, q] = G_{p,r}[U_{q,s}]  (MATLAB layout (Nx,M+1,M+1,N+1,N+1))
    G_U = np.zeros((Nx, M + 1, M + 1, N + 1, N + 1), dtype=complex)
    J_W = np.zeros((Nx, M + 1, M + 1, N + 1, N + 1), dtype=complex)

    def store(q, s):
        xi = np.zeros((Nx, M - s + 1, N - q + 1), dtype=complex)
        xi[:, 0, 0] = U[:, s, q]
        u = field_tfe_helmholtz_m_and_n(xi, f, pp, gamma_u_bar_p, alpha_bar, gamma_u_bar,
                                        Dz, a, Nx, Nz, N - q, M - s, identy, alpha_bar_p)
        G = dno_tfe_helmholtz_m_and_n(u, f, pp, Dz, a, Nx, Nz, N - q, M - s)
        G_U[:, :M - s + 1, s, :N - q + 1, q] = G
        xi = np.zeros((Nx, M - s + 1, N - q + 1), dtype=complex)
        xi[:, 0, 0] = W[:, s, q]
        w = field_tfe_helmholtz_m_and_n_lf(xi, f, pp, gamma_w_bar_p, alpha_bar, gamma_w_bar,
                                           Dz, b, Nx, Nz, N - q, M - s, identy, alpha_bar_p)
        J = dno_tfe_helmholtz_m_and_n_lf(w, f, pp, Dz, b, Nx, Nz, N - q, M - s)
        J_W[:, :M - s + 1, s, :N - q + 1, q] = J

    # n = 0, m = 0
    U[:, 0, 0], W[:, 0, 0] = AInverse(zeta_n_m[:, 0, 0], -psi_n_m[:, 0, 0],
                                      gamma_u_bar_p, gamma_w_bar_p, Nx, tau2)
    store(0, 0)

    for n in range(N + 1):
        if verbose:
            print(f'  two_layer_solve_fast: n = {n}/{N}', flush=True)
        for m in range(M + 1):
            Q = zeta_n_m[:, m, n]
            R = -psi_n_m[:, m, n] - phase_correction(U, W, m, n, alpha_bar, f_x, tau2)
            for r in range(n + 1):
                for s in range(m + 1):
                    if r < n or s < m:
                        R = R - G_U[:, m - s, s, n - r, r] - tau2 * J_W[:, m - s, s, n - r, r]
            if n > 0 or m > 0:
                U[:, m, n], W[:, m, n] = AInverse(Q, R, gamma_u_bar_p, gamma_w_bar_p, Nx, tau2)
                store(n, m)   # (MATLAB recomputes (0,0) here too; result identical)
    ubar, wbar = _traces(U, W, f, pp, gamma_u_bar_p, gamma_w_bar_p, alpha_bar, gamma_u_bar,
                         gamma_w_bar, Dz, a, b, Nx, Nz, N, M, identy, alpha_bar_p)
    return U, W, ubar, wbar
