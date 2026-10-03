"""Truth references for the 3D error analysis: HOPS solved POINTWISE in frequency (delta = 0, M = 0: no
frequency expansion), Pade in eps with N = 24, on a refined grid (2 Nx, 2 Ny, Nz + 16; an invariant
direction Ny = 1 stays 1), and the same at the scenario's own grid (summation-only error)."""
import warnings

import numpy as np

from .core import rm3                  # noqa: F401  (first: puts HOPS_Python on sys.path)
import hops3d as h3
from hops3d.grid import fft2


def hops_pointwise(sc_info, n_u, n_w, omega, alpha, beta, Eps, N=24, Nx=None, Ny=None, Nz=None):
    """R, T, D (len(Eps)) of HOPS at the single frequency omega (alpha, beta: their values AT omega)."""
    Nx = Nx or sc_info['Nx']
    Ny = Ny or sc_info['Ny']
    Nz = Nz or sc_info['Nz']
    xx = (2 * np.pi / Nx) * np.arange(Nx)
    yy = (2 * np.pi / Ny) * np.arange(Ny)
    f, fx, fy = h3.profile_fn_3d(sc_info['profile'], xx, yy)
    Mode = 1 if sc_info['mode'] == 'TE' else 2
    P = h3.make_problem(Nx, Ny, Nz, N, 0, n_u, n_w, omega, alpha, beta, f=f, f_x=fx, f_y=fy,
                        a=sc_info['a'], b=sc_info['b'], Mode=Mode)
    zeta, psi = h3.setup_zeta_psi_n_m_3d(alpha, beta, P.gamma_u_bar, P.f, P.f_x, P.f_y, N, 0)
    U, W, ubar, wbar = h3.two_layer_solve_3d_coupled(P, zeta, psi)
    return energy_eps(P, ubar[:, :, 0, :], wbar[:, :, 0, :], alpha, beta, Eps)


def energy_eps(P, ub, wb, alpha, beta, Eps, kind='pade'):
    """R, T, D from eps-series of the traces (Nx, Ny, N+1): [N/2, N/2] Pade (padesum.m) of every
    propagating Fourier amplitude (the M = 0 case, where the polar (eps, delta) summation would keep
    only order 0)."""
    from hops.summation import padesum
    Nx, Ny, N1 = ub.shape
    Eps = np.atleast_1d(np.asarray(Eps, float))
    k_u2 = alpha ** 2 + beta ** 2 + P.gamma_u_bar ** 2
    k_w2 = alpha ** 2 + beta ** 2 + P.gamma_w_bar ** 2
    ap, bq = alpha + P.kx[:, None], beta + P.ky[None, :]
    lat2 = ap ** 2 + bq ** 2
    gu, gw = h3.outgoing_sqrt(k_u2 - lat2), h3.outgoing_sqrt(k_w2 - lat2)
    pu, pw = np.real(lat2) < np.real(k_u2), np.real(lat2) < np.real(k_w2)
    R = np.zeros(Eps.size); T = np.zeros(Eps.size)
    for c3, prop, g, out, fac in ((ub, pu, gu, R, 1.0), (wb, pw, gw, T, P.tau2)):
        ch = fft2(c3) / (Nx * Ny)
        for i, j in zip(*np.nonzero(prop)):
            if kind == 'pade':
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    v = np.array([padesum(ch[i, j], float(e), (N1 - 1) // 2)[0] for e in Eps]).ravel()
            else:
                v = np.polynomial.polynomial.polyval(Eps, ch[i, j])
            out += np.real(fac * g[i, j] / gu[0, 0]) * np.abs(v) ** 2
    return R, T, 1.0 - R - T


def refined(info, fine=True):
    Nx, Ny, Nz = info['Nx'], info['Ny'], info['Nz']
    if not fine:
        return Nx, Ny, Nz
    return (2 * Nx if Nx > 1 else 1), (2 * Ny if Ny > 1 else 1), Nz + 16
