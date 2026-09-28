"""Method of Manufactured Solutions in 3D (shared by mms_error_3D.py, test_single_eps_delta_3D.py and
the tests).  Analogue of the MMS set-up in mms_error.m / test_single_eps_delta.m:

    u_rs = A exp(i p_r x + i q_s y + i gamma^u_rs z)    (upper, outgoing up),
    w_rs = B exp(i p_r x + i q_s y - i gamma^w_rs z)    (lower, outgoing down),

with gamma_rs(delta)^2 = k(delta)^2 - (alpha(delta) + p_r)^2 - (beta(delta) + q_s)^2.  Their interface
traces and normal derivatives are fed to the two-layer solver, which must return them.
"""
import numpy as np
from .setup_data import setup_xi_u_nu_u_n_m_3d, setup_xi_w_nu_w_n_m_3d
from .two_layer import phase_correction_3d


def manufactured_data(P, A, B, rs):
    """(xi_u, nu_u, xi_w, nu_w, zeta, psi) of the manufactured solution, layout (Nx, Ny, M+1, N+1).
    psi contains the Bloch-phase correction (the solver's unknowns are phase-removed)."""
    N, M = P.N, P.M
    xi_u, nu_u = setup_xi_u_nu_u_n_m_3d(A, rs, P.xx, P.yy, P.kx, P.ky, P.alphap, P.betap, P.gammap_u,
                                        P.f, P.f_x, P.f_y, N, M)
    xi_w, nu_w = setup_xi_w_nu_w_n_m_3d(B, rs, P.xx, P.yy, P.kx, P.ky, P.alphap, P.betap, P.gammap_w,
                                        P.f, P.f_x, P.f_y, N, M)
    zeta = xi_u - xi_w
    psi = -nu_u - P.tau2 * nu_w
    for n in range(N + 1):
        for m in range(M + 1):
            psi[:, :, m, n] -= phase_correction_3d(xi_u, xi_w, m, n, P.alpha_bar, P.beta_bar,
                                                   P.f_x, P.f_y, P.tau2)
    return xi_u, nu_u, xi_w, nu_w, zeta, psi


def exact_values(P, A, B, rs, Eps, delta, volume=False):
    """Exact interface quantities at (Eps, delta); Eps may be an array (-> leading axis).

    Returns dict with xi_u, nu_u (= G), xi_w, nu_w (= J), zeta, psi, U, W, G, J, ubar, wbar and, if
    volume, the fields u, w in the transformed coordinates (Nx, Ny, Nz+1)."""
    r, s = rs
    Eps = np.asarray(Eps, dtype=float)
    E = Eps[..., None, None]
    ar = P.alphap[r] + delta * P.alpha_bar
    bs = P.betap[s] + delta * P.beta_bar
    k_u = (1 + delta) * P.n_u * P.omega_bar
    k_w = (1 + delta) * P.n_w * P.omega_bar
    g_u = np.lib.scimath.sqrt(k_u ** 2 - ar ** 2 - bs ** 2)
    g_w = np.lib.scimath.sqrt(k_w ** 2 - ar ** 2 - bs ** 2)
    ph = np.exp(1j * P.kx[r] * P.xx)[:, None] * np.exp(1j * P.ky[s] * P.yy)[None, :]
    lat = 1j * (P.kx[r] * P.f_x + P.ky[s] * P.f_y)
    ex = {}
    ex['xi_u'] = A * ph * np.exp(1j * g_u * E * P.f)
    ex['nu_u'] = (-1j * g_u + lat * E) * ex['xi_u']
    ex['xi_w'] = B * ph * np.exp(-1j * g_w * E * P.f)
    ex['nu_w'] = (-1j * g_w - lat * E) * ex['xi_w']
    ex['zeta'] = ex['xi_u'] - ex['xi_w']
    # physical psi = -nu_u - tau2 nu_w - i(alpha g_x + beta g_y)(xi_u - tau2 xi_w)
    phase = 1j * (P.alpha_bar * (1 + delta) * P.f_x + P.beta_bar * (1 + delta) * P.f_y) * E
    ex['psi'] = -ex['nu_u'] - P.tau2 * ex['nu_w'] - phase * (ex['xi_u'] - P.tau2 * ex['xi_w'])
    ex['U'], ex['W'], ex['G'], ex['J'] = ex['xi_u'], ex['xi_w'], ex['nu_u'], ex['nu_w']
    ex['ubar'] = A * ph * np.exp(1j * g_u * P.a) * np.ones_like(E)
    ex['wbar'] = B * ph * np.exp(1j * g_w * P.b) * np.ones_like(E)
    if volume:
        Nz = P.Nz
        tz = np.cos(np.pi * np.arange(Nz + 1) / Nz)
        g = float(Eps) * P.f[:, :, None]
        zp = (P.a / 2.0) * (tz - 1.0) + P.a
        Z = (P.a - g) * zp / P.a + g
        ex['u'] = A * ph[:, :, None] * np.exp(1j * g_u * Z)
        zp = (P.b / 2.0) * (tz - 1.0)
        Z = (P.b + g) * zp / P.b + g
        ex['w'] = B * ph[:, :, None] * np.exp(-1j * g_w * Z)
    return ex
