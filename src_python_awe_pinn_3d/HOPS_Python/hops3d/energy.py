"""Reflectivity, transmissivity and energy defect in 3D (analogue of energy_defect.m).

For the doubly periodic grating the reflected / transmitted fields are Rayleigh expansions over
the 2D lattice of Fourier modes (p, q):
    R = sum_{propagating (p,q)} (gamma^u_pq / gamma^u_00) |ubar^_pq|^2,
    T = tau^2 sum_{propagating (p,q)} (gamma^w_pq / gamma^u_00) |wbar^_pq|^2,
    D = 1 - R - T      (= 0 for lossless layers: conservation of energy),
with ubar^_pq the Fourier coefficients of the (eps, delta)-summed trace ubar(x, y) at z = a
(the gamma e^{i gamma a} phase factors have modulus 1 for propagating modes).
A mode is propagating when Re(alpha_p^2 + beta_q^2) < Re(k^2)  (the 2D code's rule).
"""
import numpy as np
from .grid import fft2, outgoing_sqrt
from hops.summation import fcn_sum_fast


def sum_interface_3d(SumType, c_n_m, Eps, delta, N, M, taylor_full_order=False):
    """Sum an interface series (Nx, Ny, M+1, N+1) on an (eps, delta) grid.

    Eps (N_Eps,), delta (N_delta,)  ->  (N_Eps, N_delta, Nx, Ny).  Uses the 2D summation engine
    (polar re-summation + Taylor / Pade) pointwise on the flattened (x, y) grid.
    """
    Nx, Ny = c_n_m.shape[:2]
    K = Nx * Ny
    c = np.transpose(c_n_m.reshape(K, *c_n_m.shape[2:])[:, :M + 1, :N + 1], (2, 1, 0))   # (N+1, M+1, K)
    Eps = np.atleast_1d(np.asarray(Eps, dtype=float))
    delta = np.atleast_1d(np.asarray(delta, dtype=float))
    out = np.empty((Eps.size, delta.size, K), dtype=complex)
    chunk = max(1, int(4e5 // max(1, Eps.size * K)))
    for i0 in range(0, delta.size, chunk):
        out[:, i0:i0 + chunk] = fcn_sum_fast(SumType, c, Eps[:, None], delta[None, i0:i0 + chunk], K, N, M,
                                             taylor_full_order)
    return out.reshape(Eps.size, delta.size, Nx, Ny)


def _mode_arrays(kx, ky, alpha_bar, beta_bar, gamma_u_bar, gamma_w_bar, delta):
    k_u2 = alpha_bar ** 2 + beta_bar ** 2 + gamma_u_bar ** 2
    k_w2 = alpha_bar ** 2 + beta_bar ** 2 + gamma_w_bar ** 2
    s = (1 + delta)                                               # (N_delta,)
    ap = alpha_bar * s[None, None, :] + kx[:, None, None]         # (Nx, 1, N_delta)
    bq = beta_bar * s[None, None, :] + ky[None, :, None]          # (1, Ny, N_delta)
    lat2 = ap ** 2 + bq ** 2
    gu = outgoing_sqrt(k_u2 * s ** 2 - lat2)                      # (Nx, Ny, N_delta)
    gw = outgoing_sqrt(k_w2 * s ** 2 - lat2)
    prop_u = np.real(lat2) < np.real(k_u2 * s ** 2)
    prop_w = np.real(lat2) < np.real(k_w2 * s ** 2)
    return gu, gw, prop_u, prop_w


def energy_defect_3d(tau2, ubar_n_m, wbar_n_m, kx, ky, alpha_bar, beta_bar, gamma_u_bar, gamma_w_bar,
                     Eps, delta, N, M, SumType, taylor_full_order=False, return_modes=False,
                     sum_domain='fourier'):
    """ee, ru, rl of shape (N_Eps, N_delta)  (ee = 1 - ru - rl).

    sum_domain = 'physical': as energy_defect.m -- sum ubar(x, y; eps, delta) pointwise on the
                 Nx x Ny grid (Taylor / Pade), then FFT.  Cost ~ Nx Ny series per (eps, delta).
                 'fourier' (default): sum the series of the Fourier amplitudes ubar^_pq of the modes
                 that propagate somewhere in the window (a handful), directly.  For Taylor both are
                 identical (summation is linear); for Pade they are two valid re-summations that agree
                 to the accuracy of the method, and 'fourier' is ~Nx Ny / (#propagating) times cheaper.
    return_modes=True also returns per-mode efficiencies (N_Eps, N_delta, Nx, Ny) for R and T.
    """
    Nx, Ny = ubar_n_m.shape[:2]
    delta = np.atleast_1d(np.asarray(delta, dtype=float))
    Eps = np.atleast_1d(np.asarray(Eps, dtype=float))
    gu, gw, prop_u, prop_w = _mode_arrays(kx, ky, alpha_bar, beta_bar, gamma_u_bar, gamma_w_bar, delta)
    if sum_domain == 'physical':
        ub = sum_interface_3d(SumType, ubar_n_m, Eps, delta, N, M, taylor_full_order)
        wb = sum_interface_3d(SumType, wbar_n_m, Eps, delta, N, M, taylor_full_order)
        uh = fft2(np.moveaxis(ub, (2, 3), (0, 1))) / (Nx * Ny)       # (Nx, Ny, N_Eps, N_delta)
        wh = fft2(np.moveaxis(wb, (2, 3), (0, 1))) / (Nx * Ny)
    else:
        uh = np.zeros((Nx, Ny, Eps.size, delta.size), dtype=complex)
        wh = np.zeros_like(uh)
        for c_n_m, prop, out in ((ubar_n_m, prop_u, uh), (wbar_n_m, prop_w, wh)):
            sel = np.nonzero(prop.any(axis=-1))
            if sel[0].size == 0:
                continue
            ch = fft2(c_n_m)[sel] / (Nx * Ny)                         # (K, M+1, N+1)
            vals = sum_interface_3d(SumType, ch[:, None], Eps, delta, N, M, taylor_full_order)
            out[sel] = np.moveaxis(vals[:, :, :, 0], 2, 0)
    g_i = gu[0, 0]                                                # (N_delta,)
    Ru = np.where(prop_u[:, :, None, :], (gu / g_i)[:, :, None, :] * np.abs(uh) ** 2, 0.0)
    Tw = np.where(prop_w[:, :, None, :], tau2 * (gw / g_i)[:, :, None, :] * np.abs(wh) ** 2, 0.0)
    ru = Ru.sum(axis=(0, 1))
    rl = Tw.sum(axis=(0, 1))
    ee = 1.0 - ru - rl
    if return_modes:
        return ee, ru, rl, np.moveaxis(Ru, (0, 1), (2, 3)), np.moveaxis(Tw, (0, 1), (2, 3))
    return ee, ru, rl
