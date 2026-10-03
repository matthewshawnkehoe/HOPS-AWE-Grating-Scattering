"""energy_defect.m -- reflectivity (ru), transmissivity (rl) and energy defect ee = 1 - ru - rl."""
import numpy as np
from .spectral import setup_2d, csqrt
from .summation import fcn_sum_fast


def energy_defect(tau2, ubar_n_m, wbar_n_m, d, alpha_bar, gamma_u_bar, gamma_w_bar, Eps, delta,
                  Nx, N, M, N_Eps, N_delta, SumType, taylor_full_order=False):
    """ubar_n_m, wbar_n_m in the permuted layout (N+1, M+1, Nx)  (MATLAB permute(.,[3 2 1])).

    Returns ee, ru, rl of shape (N_Eps, N_delta).  Vectorised over Eps (all N_Eps
    values are summed at once for each delta) -- same arithmetic as the double loop.
    """
    Eps = np.atleast_1d(np.asarray(Eps, dtype=float))
    delta = np.atleast_1d(np.asarray(delta, dtype=float))
    ru = np.zeros((N_Eps, N_delta), dtype=complex)
    rl = np.zeros((N_Eps, N_delta), dtype=complex)
    k_u_bar = csqrt(alpha_bar ** 2 + gamma_u_bar ** 2)
    k_w_bar = csqrt(alpha_bar ** 2 + gamma_w_bar ** 2)
    # MATLAB keeps PropMode_* across iterations (only ever set to 1); reproduced here.
    PropMode_u = np.zeros(Nx)
    PropMode_w = np.zeros(Nx)
    # Sum the series for all (eps, delta) at once (in chunks of delta to bound memory).
    E2 = Eps[:N_Eps, None]
    ubar_all = np.empty((N_Eps, N_delta, Nx), dtype=complex)
    wbar_all = np.empty((N_Eps, N_delta, Nx), dtype=complex)
    chunk = max(1, int(2e5 // max(1, N_Eps * Nx)))
    for i0 in range(0, N_delta, chunk):
        dl = delta[None, i0:i0 + chunk]
        ubar_all[:, i0:i0 + chunk] = fcn_sum_fast(SumType, ubar_n_m, E2, dl, Nx, N, M, taylor_full_order)
        wbar_all[:, i0:i0 + chunk] = fcn_sum_fast(SumType, wbar_n_m, E2, dl, Nx, N, M, taylor_full_order)
    for ell in range(N_delta):
        alpha = (1 + delta[ell]) * alpha_bar
        gamma_u = (1 + delta[ell]) * gamma_u_bar
        gamma_w = (1 + delta[ell]) * gamma_w_bar
        k_u = (1 + delta[ell]) * k_u_bar
        k_w = (1 + delta[ell]) * k_w_bar
        _, _, alpha_p, gamma_u_p, _, _ = setup_2d(Nx, d, alpha, gamma_u)
        _, _, alpha_p, gamma_w_p, _, _ = setup_2d(Nx, d, alpha, gamma_w)
        gamma_i = gamma_u_p[0]
        ubar = ubar_all[:, ell]
        wbar = wbar_all[:, ell]
        ubarhat = np.fft.fft(ubar, axis=-1) / Nx          # (N_Eps, Nx)
        wbarhat = np.fft.fft(wbar, axis=-1) / Nx
        alpha_p_2 = alpha_p ** 2
        # MATLAB '<' on complex numbers compares real parts only
        PropMode_u[np.real(alpha_p_2) < np.real(k_u ** 2)] = 1
        PropMode_w[np.real(alpha_p_2) < np.real(k_w ** 2)] = 1
        Bq = gamma_u_p / gamma_i * np.abs(ubarhat) ** 2
        Cq = gamma_w_p / gamma_i * np.abs(wbarhat) ** 2
        ru[:, ell] = np.sum(PropMode_u * Bq, axis=-1)
        rl[:, ell] = tau2 * np.sum(PropMode_w * Cq, axis=-1)
    ee = 1.0 - ru - rl
    return ee, ru, rl
