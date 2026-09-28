"""Transformed Field Expansion (TFE) field solvers and DNOs:
field_tfe_helmholtz_m_and_n.m, field_tfe_helmholtz_m_and_n_lf.m,
dno_tfe_helmholtz_m_and_n.m, dno_tfe_helmholtz_m_and_n_lf.m.

Layouts:  fields u_n_m, w_n_m  -> (Nx, Nz+1, M+1, N+1)   i.e. [:, ell, m, n]
          interface data / DNOs -> (Nx, M+1, N+1)          i.e. [:, m, n]
Chebyshev index ell=0 is t=+1 (top of the layer), ell=Nz is t=-1 (bottom).
"""
import numpy as np
from .spectral import dx, dz
from .expansions import T_dno
from . import config
from .bvp import solvebvp_colloc_fast, solvebvp_colloc_fast_lf

fft = lambda u: np.fft.fft(u, axis=0)
ifft = lambda u: np.fft.ifft(u, axis=0)


def _geometry(Dz, z_min, z_max, Nz):
    ll = np.arange(Nz + 1)
    D = (2.0 / (z_max - z_min)) * Dz
    D2 = D @ D
    tilde_z = np.cos(np.pi * ll / Nz)
    z = ((z_max - z_min) / 2.0) * (tilde_z - 1.0) + z_max
    return D, D2, D[0, :], D[-1, :], z


def _volume_rhs(unm, m, n, p, Dz, L, alpha, gamma, C):
    """Right-hand side F_{n,m} of the TFE recursion (shared by both layers)."""
    A1_xx, A1_xz, A1_zx, B1_x, A2_xx, A2_xz, A2_zx, A2_zz, B2_x, B2_z, S1, S2 = C
    Fnm = np.zeros(unm.shape[:2], dtype=complex)
    g2 = gamma ** 2
    # config.ALPHA_FIX item 4: the Bloch term 2 i alpha d_x also acts through the TFE change of
    # variables, d_x -> d_x' + (dz'/dx) d_z', giving 2 i alpha A_xz d_z' (missing in MATLAB).
    fix = config.ALPHA_FIX and alpha != 0
    if n >= 1:
        u = unm[:, :, m, n - 1]
        u_x = dx(u, p)
        Fnm -= dx(A1_xx * u_x, p)
        Fnm -= dz(A1_zx * u_x, Dz, L)
        Fnm -= B1_x * u_x
        u_z = dz(u, Dz, L)
        Fnm -= dx(A1_xz * u_z, p)
        Fnm -= 2 * 1j * alpha * S1 * u_x
        Fnm -= g2 * S1 * u
        if fix:
            Fnm -= 2 * 1j * alpha * A1_xz * u_z
    if m >= 1:
        u = unm[:, :, m - 1, n]
        u_x = dx(u, p)
        Fnm -= 2 * 1j * alpha * u_x
        Fnm -= 2 * g2 * u
    if n >= 1 and m >= 1:
        u = unm[:, :, m - 1, n - 1]
        u_x = dx(u, p)
        Fnm -= 2 * 1j * alpha * S1 * u_x
        Fnm -= 2 * g2 * S1 * u
        if fix:
            Fnm -= 2 * 1j * alpha * A1_xz * dz(u, Dz, L)
    if n >= 2:
        u = unm[:, :, m, n - 2]
        u_x = dx(u, p)
        Fnm -= dx(A2_xx * u_x, p)
        Fnm -= dz(A2_zx * u_x, Dz, L)
        Fnm -= B2_x * u_x
        u_z = dz(u, Dz, L)
        Fnm -= dx(A2_xz * u_z, p)
        Fnm -= dz(A2_zz * u_z, Dz, L)
        Fnm -= B2_z * u_z
        Fnm -= 2 * 1j * alpha * S2 * u_x
        Fnm -= g2 * S2 * u
        if fix:
            Fnm -= 2 * 1j * alpha * A2_xz * u_z
    if m >= 2:
        Fnm -= g2 * unm[:, :, m - 2, n]
    if n >= 1 and m >= 2:
        Fnm -= g2 * S1 * unm[:, :, m - 2, n - 1]
    if n >= 2 and m >= 1:
        u = unm[:, :, m - 1, n - 2]
        u_x = dx(u, p)
        Fnm -= 2 * 1j * alpha * S2 * u_x
        Fnm -= 2 * g2 * S2 * u
        if fix:
            Fnm -= 2 * 1j * alpha * A2_xz * dz(u, Dz, L)
    if n >= 2 and m >= 2:
        Fnm -= g2 * S2 * unm[:, :, m - 2, n - 2]
    return Fnm


def field_tfe_helmholtz_m_and_n(xi_n_m, f, p, gammap, alpha, gamma, Dz, a, Nx, Nz, N, M,
                                identy, alphap):
    """Upper-layer field u_{n,m}(x,z) in transformed coordinates, z in [0, a]."""
    unm = np.zeros((Nx, Nz + 1, M + 1, N + 1), dtype=complex)
    k2 = alphap[0] ** 2 + gammap[0] ** 2
    ell_top = 0
    xi_hat = fft(xi_n_m)
    f = np.asarray(f)
    f_x = np.real(ifft(1j * p * fft(f)))
    D, D2, D_start, D_end, z = _geometry(Dz, 0.0, a, Nz)

    f_full = np.tile(f[:, None], (1, Nz + 1))
    f_x_full = np.tile(f_x[:, None], (1, Nz + 1))
    amz = np.tile((a - z)[None, :], (Nx, 1))

    Tu = T_dno(alpha, alphap if config.ALPHA_FIX else p, gamma, gammap, k2, Nx, M)

    unm[:, :, 0, 0] = ifft(np.exp(1j * np.outer(gammap, z)) * xi_hat[:, 0, 0][:, None])

    A1_xx = -(2.0 / a) * f_full
    A1_xz = -(1.0 / a) * amz * f_x_full
    A1_zx = A1_xz
    A2_xx = (1.0 / a ** 2) * f_full ** 2
    A2_xz = (1.0 / a ** 2) * amz * (f_full * f_x_full)
    A2_zx = A2_xz
    A2_zz = (1.0 / a ** 2) * amz ** 2 * f_x_full ** 2
    B1_x = (1.0 / a) * f_x_full
    B2_x = -(1.0 / a ** 2) * f_full * f_x_full
    B2_z = -(1.0 / a ** 2) * amz * f_x_full ** 2
    S1 = -(2.0 / a) * f_full
    S2 = (1.0 / a ** 2) * f_full ** 2
    C = (A1_xx, A1_xz, A1_zx, B1_x, A2_xx, A2_xz, A2_zx, A2_zz, B2_x, B2_z, S1, S2)

    gammagamma = k2 - alphap ** 2
    for n in range(N + 1):
        for m in range(M + 1):
            if n == 0 and m == 0:
                continue          # MATLAB solves but discards the (0,0) result
            Fnm = _volume_rhs(unm, m, n, p, Dz, a, alpha, gamma, C)
            Jnm = np.zeros(Nx, dtype=complex)
            for r in range(m):
                Jnm += ifft(Tu[:, m - r] * fft(unm[:, ell_top, r, n]))
            if n >= 1:
                for r in range(m + 1):
                    Snm = ifft(Tu[:, m - r] * fft(unm[:, ell_top, r, n - 1]))
                    Jnm -= (1.0 / a) * f * Snm
            Uhat = solvebvp_colloc_fast(fft(Fnm).T, 1.0, 0.0, gammagamma,
                                        1.0, 0.0, xi_hat[:, m, n],
                                        -1j * gammap, 1.0, fft(Jnm),
                                        Nx, identy, D, D2, D_start, D_end)
            unm[:, :, m, n] = ifft(Uhat)
    return unm


def field_tfe_helmholtz_m_and_n_lf(xi_lf_n_m, f, p, gammapw, alpha, gammaw, Dz, b, Nx, Nz, N, M,
                                   identy, alphap):
    """Lower-layer field w_{n,m}(x,z) in transformed coordinates, z in [-b, 0]."""
    wmn = np.zeros((Nx, Nz + 1, M + 1, N + 1), dtype=complex)
    k2 = alphap[0] ** 2 + gammapw[0] ** 2
    ell_bottom = Nz
    xi_hat = fft(xi_lf_n_m)
    f = np.asarray(f)
    f_x = np.real(ifft(1j * p * fft(f)))
    D, D2, D_start, D_end, z = _geometry(Dz, -b, 0.0, Nz)

    f_full = np.tile(f[:, None], (1, Nz + 1))
    f_x_full = np.tile(f_x[:, None], (1, Nz + 1))
    bpz = np.tile((b + z)[None, :], (Nx, 1))

    Tw = T_dno(alpha, alphap if config.ALPHA_FIX else p, gammaw, gammapw, k2, Nx, M)

    wmn[:, :, 0, 0] = ifft(np.exp(-1j * np.outer(gammapw, z)) * xi_hat[:, 0, 0][:, None])

    A1_xx = (2.0 / b) * f_full
    A1_xz = -(1.0 / b) * bpz * f_x_full
    A1_zx = A1_xz
    A2_xx = (1.0 / b ** 2) * f_full ** 2
    A2_xz = -(1.0 / b ** 2) * bpz * (f_full * f_x_full)
    A2_zx = A2_xz
    A2_zz = (1.0 / b ** 2) * bpz ** 2 * f_x_full ** 2
    B1_x = -(1.0 / b) * f_x_full
    B2_x = -(1.0 / b ** 2) * f_full * f_x_full
    B2_z = (1.0 / b ** 2) * bpz * f_x_full ** 2
    S1 = (2.0 / b) * f_full
    S2 = (1.0 / b ** 2) * f_full ** 2
    C = (A1_xx, A1_xz, A1_zx, B1_x, A2_xx, A2_xz, A2_zx, A2_zz, B2_x, B2_z, S1, S2)

    gammagamma = k2 - alphap ** 2
    for n in range(N + 1):
        for m in range(M + 1):
            if n == 0 and m == 0:
                continue
            Fnm = _volume_rhs(wmn, m, n, p, Dz, b, alpha, gammaw, C)
            Qnm = np.zeros(Nx, dtype=complex)
            for r in range(m):
                Qnm -= ifft(Tw[:, m - r] * fft(wmn[:, ell_bottom, r, n]))
            if n >= 1:
                for r in range(m + 1):
                    Snm = ifft(Tw[:, m - r] * fft(wmn[:, ell_bottom, r, n - 1]))
                    Qnm -= (1.0 / b) * f * Snm
            What = solvebvp_colloc_fast_lf(fft(Fnm).T, 1.0, 0.0, gammagamma,
                                           1j * gammapw, 1.0, fft(Qnm),
                                           1.0, 0.0, xi_hat[:, m, n],
                                           Nx, identy, D, D2, D_start, D_end)
            wmn[:, :, m, n] = ifft(What)
    return wmn


def dno_tfe_helmholtz_m_and_n(unm, f, p, Dz, a, Nx, Nz, N, M):
    """Upper-layer DNO  G_{n,m} = -d_N u  at the interface; shape (Nx, M+1, N+1)."""
    Gnm = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    lb = Nz
    f = np.asarray(f)
    f_x = ifft(1j * p * fft(f))
    for n in range(N + 1):
        for m in range(M + 1):
            u_z = dz(unm[:, :, m, n], Dz, a)
            Gnm[:, m, n] = -u_z[:, lb]
            if n >= 1:
                u_x = dx(unm[:, :, m, n - 1], p)
                Gnm[:, m, n] += f_x * u_x[:, lb]
                Gnm[:, m, n] += (1.0 / a) * (f * Gnm[:, m, n - 1])
            if n >= 2:
                u_x = dx(unm[:, :, m, n - 2], p)
                Gnm[:, m, n] -= (1.0 / a) * (f * (f_x * u_x[:, lb]))
                u_z = dz(unm[:, :, m, n - 2], Dz, a)
                Gnm[:, m, n] -= f_x * (f_x * u_z[:, lb])
    return Gnm


def dno_tfe_helmholtz_m_and_n_lf(wnm, f, p, Dz, b, Nx, Nz, N, M):
    """Lower-layer DNO  J_{n,m} = +d_N w  at the interface; shape (Nx, M+1, N+1)."""
    Gnm = np.zeros((Nx, M + 1, N + 1), dtype=complex)
    lt = 0
    f = np.asarray(f)
    f_x = ifft(1j * p * fft(f))
    for n in range(N + 1):
        for m in range(M + 1):
            w_z = dz(wnm[:, :, m, n], Dz, b)
            Gnm[:, m, n] = w_z[:, lt]
            if n >= 1:
                w_x = dx(wnm[:, :, m, n - 1], p)
                Gnm[:, m, n] -= f_x * w_x[:, lt]
                Gnm[:, m, n] -= (1.0 / b) * (f * Gnm[:, m, n - 1])
            if n >= 2:
                w_x = dx(wnm[:, :, m, n - 2], p)
                Gnm[:, m, n] -= (1.0 / b) * (f * (f_x * w_x[:, lt]))
                w_z = dz(wnm[:, :, m, n - 2], Dz, b)
                Gnm[:, m, n] += f_x * (f_x * w_z[:, lt])
    return Gnm
