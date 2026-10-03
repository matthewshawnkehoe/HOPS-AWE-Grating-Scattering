"""Interface data in 3D: analogues of setup_xi_u_nu_u_n_m.m, setup_xi_w_nu_w_n_m.m,
setup_zeta_psi_n_m.m, plus doubly periodic grating profiles.

All (n, m) arrays have the layout (Nx, Ny, M+1, N+1): [:, :, m, n].

Sign conventions are those of the 2D code (see hops/two_layer.py):
    G = -d_N u (upper DNO),  J = +d_N w (lower DNO),  N = (-g_x, -g_y, 1),
and the manufactured / computed DNOs use the PHASE-REMOVED derivatives (d_x, d_y act on
u~ = e^{-i(alpha x + beta y)} u).  The Bloch-phase part of the physical normal derivative,
i(alpha g_x + beta g_y)(U - tau^2 W), is added in the interface equation
(two_layer.phase_correction_3d), exactly as config.ALPHA_FIX item 3 does in 2D.
"""
import numpy as np
from .expansions import gamma_exp_3d, E_exp_3d
from .grid import csqrt
from hops.setup_data import fourier_repr_lipschitz, fourier_repr_rough


def _mode_data(r, s, kx, ky, alphap, betap, gammap):
    alpha_bar, beta_bar = alphap[0], betap[0]
    gamma_bar = gammap[0, 0]
    k_bar = csqrt(alpha_bar ** 2 + beta_bar ** 2 + gamma_bar ** 2)
    return alpha_bar, beta_bar, gamma_bar, k_bar


def setup_xi_u_nu_u_n_m_3d(A, rs, xx, yy, kx, ky, alphap, betap, gammap, f, f_x, f_y, N, M):
    """Upper-layer manufactured data  u_rs = A exp(i p_r x + i q_s y + i gamma_rs z)  (outgoing up).

    xi = u_rs(x, y, g),   nu = -d_N u_rs = (-i gamma_rs + i (p_r g_x + q_s g_y)) xi.
    """
    r, s = rs
    ab, bb, gb, kb = _mode_data(r, s, kx, ky, alphap, betap, gammap)
    g_m = gamma_exp_3d(ab, bb, alphap[r], betap[s], gb, gammap[r, s], kb, M)
    E = E_exp_3d(g_m, f, N, M, +1.0)
    ph = A * np.exp(1j * kx[r] * xx)[:, None] * np.exp(1j * ky[s] * yy)[None, :]
    lat = 1j * (kx[r] * f_x + ky[s] * f_y)
    xi = ph[..., None, None] * E
    nu = np.zeros_like(xi)
    for n in range(N + 1):
        for m in range(M + 1):
            for ell in range(m + 1):
                nu[..., m, n] += (-1j * g_m[m - ell]) * xi[..., ell, n]
            if n > 0:
                nu[..., m, n] += lat * xi[..., m, n - 1]
    return xi, nu


def setup_xi_w_nu_w_n_m_3d(B, rs, xx, yy, kx, ky, alphap, betap, gammap, f, f_x, f_y, N, M):
    """Lower-layer manufactured data  w_rs = B exp(i p_r x + i q_s y - i gamma_rs z)  (outgoing down).

    xi = w_rs(x, y, g),   nu = (-i gamma_rs - i (p_r g_x + q_s g_y)) xi   (= +d_N w; the 2D code's
    'nu_w' convention).
    """
    r, s = rs
    ab, bb, gb, kb = _mode_data(r, s, kx, ky, alphap, betap, gammap)
    g_m = gamma_exp_3d(ab, bb, alphap[r], betap[s], gb, gammap[r, s], kb, M)
    E = E_exp_3d(g_m, f, N, M, -1.0)
    ph = B * np.exp(1j * kx[r] * xx)[:, None] * np.exp(1j * ky[s] * yy)[None, :]
    lat = 1j * (kx[r] * f_x + ky[s] * f_y)
    xi = ph[..., None, None] * E
    nu = np.zeros_like(xi)
    for n in range(N + 1):
        for m in range(M + 1):
            for ell in range(m + 1):
                nu[..., m, n] -= (1j * g_m[m - ell]) * xi[..., ell, n]
            if n > 0:
                nu[..., m, n] -= lat * xi[..., m, n - 1]
    return xi, nu


def setup_zeta_psi_n_m_3d(alpha_bar, beta_bar, gamma_u_bar, f, f_x, f_y, N, M):
    """Plane-wave incidence  u_inc = exp(i alpha x + i beta y - i gamma^u z)  (analogue of
    setup_zeta_psi_n_m.m with the full alpha(delta) convolution, i.e. config.ALPHA_FIX item 2).

    zeta = -u~_inc(g) = -exp(-i gamma eps f),
    psi  = d_N u_inc (physical)  = (i gamma + i alpha g_x + i beta g_y) exp(-i gamma eps f),
    with gamma(delta) = gamma_bar (1+delta), alpha(delta) = alpha_bar (1+delta), beta likewise.
    """
    L = max(M, 1) + 1
    g_m = np.zeros(L, dtype=complex); g_m[:2] = gamma_u_bar
    a_m = np.zeros(L, dtype=complex); a_m[:2] = alpha_bar
    b_m = np.zeros(L, dtype=complex); b_m[:2] = beta_bar
    E = E_exp_3d(g_m, f, N, M, -1.0)
    zeta = -E.copy()
    psi = np.zeros_like(E)
    for n in range(N + 1):
        for m in range(M + 1):
            for ell in range(m + 1):
                psi[..., m, n] += (1j * g_m[m - ell]) * E[..., ell, n]
                if n > 0:
                    psi[..., m, n] += (f_x * (1j * a_m[m - ell]) + f_y * (1j * b_m[m - ell])) * E[..., ell, n - 1]
    return zeta, psi


# ----------------------------------------------------------------------------
# doubly periodic profiles f(x, y)  ->  f, f_x, f_y  on the (Nx, Ny) grid
# ----------------------------------------------------------------------------
PROFILES_3D = ('cosx', 'fs1d', 'cosxcosy', 'cosx+cosy', 'cos4x', 'cos4xcos4y', 'cos4x+cos4y', 'fs', 'fs_sum',
               'sin4xsin4y', 'cos2xcosy', 'egg', 'rough', 'lipschitz', 'rough120', 'lipschitz120', 'rough:<P>', 'lipschitz:<P>',
               'expr:<numpy expression in x, y>')


def profile_fn_3d(name, xx, yy):
    """Grating profiles f(x, y) with spectral derivatives.  1D names ('cosx', 'cos4x', ...) give
    y-invariant (classical, 'singly periodic') gratings -> reduce exactly to the 2D code."""
    X, Y = np.meshgrid(xx, yy, indexing='ij')
    c, s = np.cos, np.sin
    tab = {
        'cosx': (c(X), -s(X), 0 * X),
        'cos4x': (c(4 * X), -4 * s(4 * X), 0 * X),
        'fs1d': (c(4 * X) / 4, -s(4 * X), 0 * X),                       # paper f_s, y-invariant
        'cosxcosy': (c(X) * c(Y), -s(X) * c(Y), -c(X) * s(Y)),
        'cosx+cosy': ((c(X) + c(Y)) / 2, -s(X) / 2, -s(Y) / 2),
        'cos4xcos4y': (c(4 * X) * c(4 * Y), -4 * s(4 * X) * c(4 * Y), -4 * c(4 * X) * s(4 * Y)),
        'cos4x+cos4y': ((c(4 * X) + c(4 * Y)) / 2, -2 * s(4 * X), -2 * s(4 * Y)),
        'fs': (c(4 * X) * c(4 * Y) / 4, -s(4 * X) * c(4 * Y), -c(4 * X) * s(4 * Y)),   # 3D f_s
        'fs_sum': ((c(4 * X) + c(4 * Y)) / 8, -s(4 * X) / 2, -s(4 * Y) / 2),
        'sin4xsin4y': (s(4 * X) * s(4 * Y), 4 * c(4 * X) * s(4 * Y), 4 * s(4 * X) * c(4 * Y)),
        'cos2xcosy': (c(2 * X) * c(Y), -2 * s(2 * X) * c(Y), -c(2 * X) * s(Y)),
        'egg': ((c(X) + c(Y) + c(X) * c(Y)) / 3, -(s(X) + s(X) * c(Y)) / 3, -(s(Y) + c(X) * s(Y)) / 3),
    }
    if name in tab:
        return tab[name]
    cases = [('rough', 40), ('lipschitz', 40), ('rough120', 120), ('lipschitz120', 120)]
    if name.startswith(('rough:', 'lipschitz:')):          # e.g. 'rough:20' = f_{r,P} with P = 20 terms
        cases = [(name, int(name.split(':')[1]))]
    for base, P in cases:
        if name == base:
            fn = fourier_repr_rough if base.startswith('rough') else fourier_repr_lipschitz
            fx1, dfx1 = fn(P, xx)
            fy1, dfy1 = fn(P, yy)
            # 'sum' construction keeps the regularity class of the 1D profile (paper eq. (38))
            return (0.5 * (fx1[:, None] + fy1[None, :]), 0.5 * dfx1[:, None] + 0 * Y,
                    0.5 * dfy1[None, :] + 0 * X)
    if name.startswith('expr:'):
        ns = {k: getattr(np, k) for k in ('sin', 'cos', 'exp', 'tanh', 'sinh', 'cosh', 'abs', 'pi', 'sqrt')}
        ns.update(x=X, y=Y)
        f = np.asarray(eval(name[5:], {'__builtins__': {}}, ns), dtype=float) * np.ones_like(X)
        f = f - f.mean()
        from .grid import grad_profile
        kx = np.fft.fftfreq(xx.size, d=(xx[1] - xx[0]) if xx.size > 1 else 1.0) * 2 * np.pi
        ky = np.fft.fftfreq(yy.size, d=(yy[1] - yy[0]) if yy.size > 1 else 1.0) * 2 * np.pi
        fx, fy = grad_profile(f, kx, ky)
        return f, fx, fy
    raise ValueError(f'unknown 3D profile {name}; choose one of {PROFILES_3D}')
