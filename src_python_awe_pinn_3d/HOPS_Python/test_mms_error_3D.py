"""test_mms_error_3D.py -- rigorous validation of the 3D HOPS/AWE code (package hops3d); 3D analogue of
test_mms_error.py.  Runs every check in tests_3d/test_hops3d.py and prints a PASS/FAIL table.

  E1  y-invariant gratings reduce EXACTLY (round-off) to the validated 2D solver, alpha = 0 and 0.1
  E2  y-invariant reflectivity / transmissivity == 2D energy_defect
  E3  3D solvers agree: coupled (default) == lean (= two_layer_solve_fast ordering) == operator
  E4  manufactured doubly periodic solution recovered: crossed profile, oblique alpha, beta != 0
  E5  test_single_eps_delta_3D (all fields, DNOs, U, W, traces; Taylor/Pade/Pade-safe) ~1e-13
  E6  flat interface: Fresnel coefficients and R + T = 1 (oblique, TE/TM)
  E7  energy conservation for a crossed dielectric grating, |D| < 1e-11 (Pade)
  E8  x <-> y symmetry for a symmetric profile and alpha = beta
  E9  frequency expansion of i gamma_pq(delta) (T_dno_3d)
  E10 Rayleigh-frequency windows (3D lattice, fixed angle; Ny = 1 gives the 2D bands)
  E11 Fourier- and physical-domain summation of R agree
  E12 conical incidence on a y-invariant grating: energy conservation

Equivalent:  pytest -v tests_3d          (PyCharm: right-click this file -> Run 'pytest in ...')
"""
import time
import traceback

import matplotlib
matplotlib.use('Agg')

from tests_3d import test_hops3d as T
from tests_3d.test_hops3d import *  # noqa: F401,F403,E402   (pytest collects these when run on this file)

CHECKS = [
    ('E1 y-invariant 3D == 2D solver (alpha = 0, 0.1; Ny = 1, 4)', T.test_y_invariant_reduces_to_2d),
    ('E2 y-invariant 3D reflectivity == 2D energy_defect', T.test_y_invariant_energy_equals_2d),
    ('E3 coupled == lean == operator 3D solvers', T.test_solvers_agree),
    ('E4 manufactured solution, crossed egg-crate profile, alpha = 0.13, beta = 0.21', T.test_mms_crossed_oblique),
    ('E5 test_single_eps_delta_3D: fields, DNOs, U, W, ubar, wbar', T.test_single_eps_delta_3d_fields),
    ('E6 flat interface: Fresnel, R + T = 1', T.test_flat_interface_fresnel),
    ('E7 energy conservation, crossed dielectric grating (Pade)', T.test_energy_conservation_crossed_dielectric),
    ('E8 x <-> y symmetry', T.test_xy_symmetry),
    ('E9 T_dno_3d frequency expansion', T.test_T_dno_3d_expansion),
    ('E10 Rayleigh windows (3D lattice, angles, 2D bands)', T.test_rayleigh_windows),
    ('E11 Fourier- vs physical-domain summation', T.test_fourier_vs_physical_summation),
    ('E12 conical incidence (beta != 0) on a 1D grating: energy conservation', T.test_conical_incidence_energy),
]

if __name__ == '__main__':
    results = []
    for name, fn in CHECKS:
        t0 = time.time()
        try:
            fn()
            status = 'PASS'
        except Exception:
            status = 'FAIL'
            traceback.print_exc()
        results.append((status, name, time.time() - t0))
        print(f'[{status}] {name}  ({time.time() - t0:.1f} s)', flush=True)
    n_fail = sum(r[0] == 'FAIL' for r in results)
    print(f'\n{len(results) - n_fail}/{len(results)} checks passed' + ('' if n_fail == 0 else f', {n_fail} FAILED'))
