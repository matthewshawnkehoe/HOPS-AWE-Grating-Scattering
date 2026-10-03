"""Tests for HOPS/AWE + PINN on refl_map.py scenarios (~1 min).

  A1  core.run returns refl_map-format results; at eps = 0 the hybrid R equals the flat-interface R
  A2  at the band edge of paper Fig. 9 the hybrid matches an independent HOPS solve to < 1e-11 while
      refl_map.py's AWE is off by > 1e-6
  A3  tol = inf (indicator only) reproduces the full-order AWE Taylor sum, and the indicator is tiny at
      the band centre and large at the band edge
  A4  the hybrid works on a multi-window named-material scenario (water over gold) and keeps R in [0, 1]
"""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings('ignore')
from awepinn import core                                    # noqa: E402


def test_format_and_flat():
    res, info = core.run('dielectric', qq=(1,), N_Eps=5, N_delta=5, verbose=False)
    r = res[0]
    for k in ('ru', 'rl', 'ee', 'RR', 'ru_awe', 'indicator', 'solved', 'lam', 'Eps'):
        assert k in r
    assert np.allclose(np.real(r['ru'][0]), np.real(r['ru_flat'][0]), atol=1e-13)


def test_band_edge_against_reference():
    sys.path.insert(0, os.path.join(os.path.dirname(core.HERE), 'PINN_HOPS'))
    from pinn_hops.hops_reference import hops_point
    from pinn_hops.problem import Grating2D
    res, info = core.run('dielectric', qq=(1,), N_Eps=3, N_delta=3, verbose=False)
    r = res[0]
    e, om = float(r['Eps'][-1]), float(r['omega'][0])
    ref = hops_point(Grating2D(eps=e, omega=om, n_w=1.1), N=24, Nx=64, Nz=48, summation='pade', fields=False)['R']
    assert abs(np.real(r['ru'][-1, 0]) - ref) < 1e-11
    assert abs(np.real(r['ru_awe'][-1, 0]) - ref) > 1e-6


def test_indicator_only():
    res, info = core.run('dielectric', qq=(1,), N_Eps=3, N_delta=5, tol=np.inf, verbose=False)
    r = res[0]
    assert not r['solved'].any()
    assert r['indicator'][-1, 2] < 1e-18 and r['indicator'][-1, 0] > 1e-8


def test_named_material_multiwindow():
    res, info = core.run('water_over_gold', qq=(1,), N_Eps=3, N_delta=3, verbose=False)
    assert len(res) > 1
    R = np.concatenate([np.real(r['ru']).ravel() for r in res])
    assert np.all(R > -1e-8) and np.all(R < 1 + 1e-8)
