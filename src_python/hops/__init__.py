"""HOPS/AWE grating-scattering solver -- Python port of the MATLAB code of
M. Kehoe & D. P. Nicholls, "A stable HOPS/AWE method for the numerical solution
of grating scattering problems", J. Sci. Comput. 100:9 (2024).

Every MATLAB function has a Python function with the same name.
"""
from .spectral import cheb, dx, dz, setup_2d, csqrt
from .expansions import gamma_exp, E_exp, E_exp_lf, A_exp, T_dno
from .setup_data import (setup_xi_u_nu_u_n_m, setup_xi_w_nu_w_n_m, setup_zeta_psi_n_m,
                         fourier_repr_lipschitz, fourier_repr_rough)
from .bvp import solvebvp_colloc_fast, solvebvp_colloc_fast_lf
from .fields import (field_tfe_helmholtz_m_and_n, field_tfe_helmholtz_m_and_n_lf,
                     dno_tfe_helmholtz_m_and_n, dno_tfe_helmholtz_m_and_n_lf)
from .two_layer import AInverse, two_layer_solve, two_layer_solve_fast
from .coupled import two_layer_solve_coupled
from .summation import (taylorsum, taylorsum2, taylorsum_2_coeff, padesum, padesum_safe,
                        padesum2, padesum2_safe, padeapprox, padesum_robust,
                        fcn_sum, fcn_sum_fast, vol_fcn_sum, sum_series)
from .energy import energy_defect
from .plotting import plot_errors

__version__ = '1.0.0'
