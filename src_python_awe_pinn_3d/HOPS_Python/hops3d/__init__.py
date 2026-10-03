"""hops3d -- the 3D (doubly periodic, two-layer, scalar Helmholtz) analogue of the HOPS/AWE
joint boundary/frequency perturbation method of

    M. Kehoe & D. P. Nicholls, "A stable HOPS/AWE method for the numerical solution of
    grating scattering problems", J. Sci. Comput. 100:9 (2024)   (2D: package ``hops``).

The field is expanded jointly in the interface height eps (g = eps f(x, y)) and the frequency
perturbation delta (omega = omega_bar (1 + delta)):

    u(x, y, z; eps, delta) = sum_{n,m} u_{n,m}(x, y, z) eps^n delta^m,
    R(eps, delta)          = sum_{n,m} R_{n,m} eps^n delta^m,

and the series are summed with the same polar Taylor / Pade re-summation as in 2D
(hops.summation).  Physics: two homogeneous layers (indices n^u above, n^w below) separated by
a crossed grating, plane wave incidence exp(i(alpha x + beta y - gamma^u z)); tau^2 = 1 or
(n^u/n^w)^2 in the transmission condition ("TE"/"TM" in 2D; for crossed gratings this is the
scalar model -- full vector Maxwell requires coupled polarisations).
"""
import numpy as np

from .grid import (setup_3d, dx3, dy3, dz3, grad_profile, cheb, csqrt, outgoing_sqrt, fft2, ifft2,
                   wavenumbers)
from .expansions import gamma_exp_3d, T_dno_3d, E_exp_3d
from .setup_data import (setup_xi_u_nu_u_n_m_3d, setup_xi_w_nu_w_n_m_3d, setup_zeta_psi_n_m_3d,
                         profile_fn_3d, PROFILES_3D)
from .layer import layer_setup, layer_run, LayerTFE, field_tfe_helmholtz_3d, dno_tfe_helmholtz_3d
from .two_layer import (AInverse_3d, phase_correction_3d, Problem3D, two_layer_solve_3d,
                        two_layer_solve_3d_operator, two_layer_solve_3d_coupled, two_layer_solve_3d_auto, layer_operators_3d, SOLVERS_3D)
from .energy import energy_defect_3d, sum_interface_3d

__version__ = '1.0.0'


def make_problem(Nx, Ny, Nz, N, M, n_u, n_w, omega_bar, alpha_bar=0.0, beta_bar=0.0, f=None,
                 f_x=None, f_y=None, profile=None, a=1.0, b=1.0, Mode=2, dx_period=2 * np.pi,
                 dy_period=2 * np.pi, c_0=1.0):
    """Assemble a Problem3D for one frequency window (omega_bar, lateral wavenumbers alpha_bar,
    beta_bar of the incident wave, i.e. k^u_bar sin(theta) cos(phi), k^u_bar sin(theta) sin(phi))."""
    k_u = n_u * omega_bar / c_0
    k_w = n_w * omega_bar / c_0
    gamma_u_bar = csqrt(k_u ** 2 - alpha_bar ** 2 - beta_bar ** 2)
    gamma_w_bar = csqrt(k_w ** 2 - alpha_bar ** 2 - beta_bar ** 2)
    xx, yy, kx, ky, alphap, betap, gammap_u = setup_3d(Nx, Ny, dx_period, dy_period, alpha_bar, beta_bar,
                                                       gamma_u_bar)
    _, _, _, _, _, _, gammap_w = setup_3d(Nx, Ny, dx_period, dy_period, alpha_bar, beta_bar, gamma_w_bar)
    if f is None:
        f, f_x, f_y = profile_fn_3d(profile, xx, yy)
    Dz, _ = cheb(Nz)
    tau2 = 1.0 if Mode == 1 else (n_u / n_w) ** 2
    P = Problem3D(f, f_x, f_y, kx, ky, alphap, betap, gammap_u, gammap_w, alpha_bar, beta_bar,
                  gamma_u_bar, gamma_w_bar, Dz, a, b, Nz, N, M, tau2)
    P.xx, P.yy, P.n_u, P.n_w, P.omega_bar = xx, yy, n_u, n_w, omega_bar
    return P
