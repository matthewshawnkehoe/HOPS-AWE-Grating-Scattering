"""Interface data / profiles: setup_xi_u_nu_u_n_m.m, setup_xi_w_nu_w_n_m.m,
setup_zeta_psi_n_m.m, fourier_repr_lipschitz.m, fourier_repr_rough.m.

All (n,m) arrays have the MATLAB layout (Nx, M+1, N+1): [:, m, n].
"""
import numpy as np
from .expansions import gamma_exp, E_exp, E_exp_lf
from .spectral import csqrt
from . import config


def setup_xi_u_nu_u_n_m(A, r, xx, pp, alpha_bar_p, gamma_bar_p, f, f_x, Nx, N, M):
    """Upper-layer manufactured data  u_r = A exp(i p_r x + i gamma_r z)."""
    alpha_bar = alpha_bar_p[0]
    gamma_bar = gamma_bar_p[0]
    k_bar = csqrt(alpha_bar ** 2 + gamma_bar ** 2)
    pp_r = pp[r]
    gamma_r_m = gamma_exp(alpha_bar, alpha_bar_p[r], gamma_bar, gamma_bar_p[r], k_bar, M)
    E = E_exp(gamma_r_m, f, N, M)
    upper = A * np.exp(1j * pp_r * xx)
    xi = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    nu = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    for n in range(N + 1):
        for m in range(M + 1):
            xi[:, m, n] = upper * E[:, m, n]
            for ell in range(m + 1):
                nu[:, m, n] += (-1j * gamma_r_m[m - ell]) * xi[:, ell, n]
            if n > 0:
                nu[:, m, n] += f_x * (1j * pp_r) * xi[:, m, n - 1]
    return xi, nu


def setup_xi_w_nu_w_n_m(B, r, xx, pp, alpha_bar_p, gamma_bar_p, f, f_x, Nx, N, M):
    """Lower-layer manufactured data  w_r = B exp(i p_r x - i gamma_r z)  (uses +d_N w)."""
    alpha_bar = alpha_bar_p[0]
    gamma_bar = gamma_bar_p[0]
    k_bar = csqrt(alpha_bar ** 2 + gamma_bar ** 2)
    pp_r = pp[r]
    gamma_r_m = gamma_exp(alpha_bar, alpha_bar_p[r], gamma_bar, gamma_bar_p[r], k_bar, M)
    E = E_exp_lf(gamma_r_m, f, N, M)
    lower = B * np.exp(1j * pp_r * xx)
    xi = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    nu = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    for n in range(N + 1):
        for m in range(M + 1):
            xi[:, m, n] = lower * E[:, m, n]
            for ell in range(m + 1):
                nu[:, m, n] -= lower * (1j * gamma_r_m[m - ell]) * E[:, ell, n]
            if n > 0:
                nu[:, m, n] -= lower * f_x * (1j * pp_r) * E[:, m, n - 1]
    return xi, nu


def setup_zeta_psi_n_m(xx, pp, alpha_bar, gamma_u_bar, f, f_x, Nx, N, M):
    """Plane-wave incidence data zeta_{n,m}, psi_{n,m} used by refl_map.

    Note (faithful to MATLAB): after ``for ell=0:m`` the loop variable equals m,
    so the n>0 term uses alpha_u_m(m-ell+1) = alpha_u_m(1) = alpha_bar.
    """
    L = max(M, 1) + 1
    gamma_u_m = np.zeros(L, dtype=complex)
    gamma_u_m[0] = gamma_u_bar
    gamma_u_m[1] = gamma_u_bar
    alpha_u_m = np.zeros(L, dtype=complex)
    alpha_u_m[0] = alpha_bar
    alpha_u_m[1] = alpha_bar
    E = E_exp_lf(gamma_u_m, f, N, M)
    zeta = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    psi = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    for n in range(N + 1):
        for m in range(M + 1):
            zeta[:, m, n] = -E[:, m, n]
            for ell in range(m + 1):
                psi[:, m, n] += (1j * gamma_u_m[m - ell]) * E[:, ell, n]
            if n > 0:
                if config.ALPHA_FIX:     # alpha(delta) = alpha_bar (1 + delta): full convolution
                    for ell in range(m + 1):
                        psi[:, m, n] += f_x * (1j * alpha_u_m[m - ell]) * E[:, ell, n - 1]
                else:                    # MATLAB: loop variable ell = m after the loop
                    ell = m
                    psi[:, m, n] += f_x * (1j * alpha_u_m[m - ell]) * E[:, ell, n - 1]
    return zeta, psi


def fourier_repr_lipschitz(P, x, plot=False):
    """Truncated Fourier series of the Lipschitz (sawtooth) profile f_L, eq. (38b)."""
    x = np.asarray(x, dtype=float)
    f = np.zeros_like(x)
    f_x = np.zeros_like(x)
    for k in range(1, int(np.ceil(P / 2)) + 1):
        f += 8 / (np.pi ** 2 * (2 * k - 1) ** 2) * np.cos((2 * k - 1) * x)
        f_x += (-8 * (2 * k - 1)) / (np.pi ** 2 * (2 * k - 1) ** 2) * np.sin((2 * k - 1) * x)
    if plot:
        _plot_profile(x, f, f_x, 'Lipschitz', 3)
    return f, f_x


def fourier_repr_rough(P, x, plot=False):
    """Truncated Fourier series of the rough (C^4) profile f_r, eq. (38a)."""
    x = np.asarray(x, dtype=float)
    f = np.zeros_like(x)
    f_x = np.zeros_like(x)
    for k in range(1, P + 1):
        c = 96 * (2 * k ** 2 * np.pi ** 2 - 21) / (125 * k ** 8)
        f += c * np.cos(k * x)
        f_x += -k * c * np.sin(k * x)
    if plot:
        _plot_profile(x, f, f_x, 'Rough', 4)
    return f, f_x


def _plot_profile(x, f, f_x, name, num):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, num=num, figsize=(9, 3.5))
    ax[0].plot(x, f); ax[0].set_title(f'{name} profile for $f$')
    ax[1].plot(x, f_x); ax[1].set_title(f'{name} profile for $f_x$')
    fig.tight_layout()
