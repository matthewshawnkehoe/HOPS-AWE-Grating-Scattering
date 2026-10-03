"""test_mms_error.py -- rigorous validation of the Python HOPS/AWE port.

(The GitHub README refers to a ``test_mms_error.m``; that file is not in the
repository -- its role is played by ``test_single_eps_delta.m``, ported as
``test_single_eps_delta.py``.  This script is the corresponding validation driver.)

It runs every check in tests/test_against_matlab.py and prints a PASS/FAIL table:
  Group A  bit-level comparison with the ORIGINAL MATLAB code (run under Octave)
           for every building block, the full test_single_eps_delta pipeline,
           a reduced mms_error.m run and reduced refl_map.m runs (silver + dielectric).
  Group B  mathematical checks: exact Chebyshev derivative (no sign flip), fast ==
           slow two-layer solver, MMS recovery to ~1e-12, energy conservation,
           Pade robustness for degenerate coefficients.

Equivalent:  pytest -v tests
"""
import time
import traceback

import matplotlib
matplotlib.use('Agg')

from tests import test_against_matlab as T
# When PyCharm/pytest runs this file (files named test_*.py are collected by pytest),
# expose every check of the suite as a pytest test:
from tests.test_against_matlab import *  # noqa: F401,F403,E402

CHECKS = [
    ('A1 building blocks vs MATLAB (cheb, dx, dz, setup_2d, expansions, T_dno, AInverse, sums, Pade)',
     T.test_units_vs_matlab),
    ('A2 fields/DNOs/two-layer solve vs MATLAB (test_single_eps_delta config)', T.test_single_pipeline_vs_matlab),
    ('A3 test_single_eps_delta error tables vs MATLAB (all n,m, 3 summations)', T.test_single_errors_vs_matlab),
    ('A4 mms_error.m (reduced grid) vs MATLAB', T.test_mms_vs_matlab),
    ('A5 refl_map.m silver (reduced grid) vs MATLAB', T.test_refl_silver_vs_matlab),
    ('A6 refl_map dielectric / test_scenarios.m (reduced grid) vs MATLAB', T.test_refl_dielectric_vs_matlab),
    ('B1 Chebyshev dz is exact; negating D is wrong', T.test_dz_exact_derivative),
    ('B2 two_layer_solve == two_layer_solve_fast', T.test_fast_equals_slow_two_layer),
    ('B3 manufactured solution recovered (U,W,G,J,ubar,wbar,u,w)', T.test_mms_recovers_exact_solution),
    ('B4 energy conservation R+T=1 for a dielectric', T.test_energy_conservation_dielectric),
    ('B5 Pade robustness (tiny / zero coefficients)', T.test_pade_handles_tiny_and_zero_coefficients),
    ('B6 operator solver == two_layer_solve_fast (dense + FFT paths, alpha = 0 and 0.1)',
     T.test_operator_solver_equals_classic),
    ('B7 coupled solver == two_layer_solve_fast (alpha = 0 and 0.1)', T.test_coupled_solver_equals_classic),
    ('A9 refl_map gold (paper Fig. 10b, reduced grid) vs MATLAB; D = 1 - R for a metal', T.test_refl_gold_vs_matlab),
    ('A8 FULL paper Fig. 9 map (6 bands, 100x100): R and D vs MATLAB/Octave', T.test_refl_dielectric_full_map_vs_matlab),
    ('A7 refl_map.m silver with the classic (MATLAB) solver vs MATLAB', T.test_refl_silver_vs_matlab_classic_solver),
    ('C1 refractiveindex.info reader (Sellmeier + tabulated n,k)', T.test_refractiveindex_reader),
    ('C2 grating-coupled surface plasmon position (Ag, P = 0.5 um) vs. SPP condition', T.test_surface_plasmon_position),
    ('D1 oblique incidence (alpha = 0.1): manufactured solution, both layers', T.test_mms_oblique_incidence),
    ('D2 oblique incidence (alpha = 0.01, paper Fig. 14): energy defect as for alpha = 0', T.test_energy_defect_oblique_dielectric),
]

if __name__ == '__main__':
    results = []
    for name, fn in CHECKS:
        t0 = time.time()
        try:
            fn()
            status = 'PASS'
        except Exception as exc:                       # includes pytest.skip
            status = 'SKIP' if type(exc).__name__ == 'Skipped' else 'FAIL'
            if status == 'FAIL':
                traceback.print_exc()
        results.append((status, name, time.time() - t0))
        print(f'[{status}] {name}  ({time.time() - t0:.1f} s)', flush=True)
    n_fail = sum(r[0] == 'FAIL' for r in results)
    print(f'\n{len(results) - n_fail}/{len(results)} checks passed' + ('' if n_fail == 0 else f', {n_fail} FAILED'))
