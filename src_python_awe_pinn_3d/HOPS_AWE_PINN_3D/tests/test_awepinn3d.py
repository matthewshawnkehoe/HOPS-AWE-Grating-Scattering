"""Tests of the 3D HOPS/AWE + PINN (physics-informed summation).   python -m pytest tests -q"""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings('ignore')
from awepinn3d import core, reference as ref        # noqa: E402
import hops3d as h3                                  # noqa: E402


def _one(scenario, n=3, **kw):
    res, info = core.run(scenario, qq=(1,), N_Eps=n, N_delta=n, keep_window=True, tol=-1, **kw)
    return res[0], info


def test_flat_interface_and_taylor_consistency():
    """eps = 0: hybrid R = AWE R (Fresnel); the indicator's Taylor-sum R (from the volume fields) equals the
    full double Taylor sum of the solver's traces ubar_{n,m} (refl_map_3D energy formula)"""
    r, info = _one('dielectric')
    assert np.max(np.abs(np.real(r['ru'][0]) - np.real(r['ru_flat'][0]))) < 1e-12
    W = r['W']
    e, d = 0.1, 0.5 * r['dmax']
    s = W.indicator(e, r['omega_bar'] * (1 + d))
    tw = W.taylor_weights(e, d).reshape(info['M'] + 1, info['N'] + 1)
    uh = h3.fft2(np.einsum('xymn,mn->xy', r['ubar_n_m'], tw)) / (W.Nx * W.Ny)
    gu, pu = W._gamma_pq('u', 1 + d)
    R_direct = float(np.sum(np.where(pu, np.real(gu / gu[0, 0]) * np.abs(uh) ** 2, 0)))
    assert abs(s['R_taylor'] - R_direct) < 1e-12


def test_dielectric_energy_conservation():
    """lossless crossed grating: |D| of the hybrid at round-off level, far below the AWE map's"""
    r, _ = _one('dielectric', n=4)
    assert np.max(np.abs(r['ee'])) < 1e-11
    assert np.max(np.abs(r['ee_awe'])) > 1e3 * np.max(np.abs(r['ee']))


def test_discrete_form_no_aliasing():
    """cos4x cos4y on 32 x 32 (under-resolved products): the Taylor sum has a tiny residual in the discrete
    (flux) form but not in the chain-rule form -- the reason the discrete form is the default"""
    r, _ = _one('gold', n=3)
    W = r['W']
    ind = W.indicator(0.03, r['omega_bar'])['loss_awe_taylor']
    assert ind < 1e-15
    res_c, _ = core.run('gold', qq=(1,), N_Eps=3, N_delta=3, keep_window=True, tol=-1, form='chain')
    ind_c = res_c[0]['W'].indicator(0.03, r['omega_bar'])['loss_awe_taylor']
    assert ind_c > 1e3 * ind


def test_summation_error_vs_same_grid_reference():
    """y-invariant silver (Ny = 1): at eps = 0.2 the hybrid is closer to the same-grid pointwise HOPS than AWE"""
    r, info = _one('silver_1d', n=3)
    j = 2
    om = float(r['omega'][j])
    R, T, D = ref.hops_pointwise(info, r['n_u'], r['n_w'], om, 0.0, 0.0, r['Eps'])
    e_awe = abs(np.real(r['ru_awe'][-1, j]) - R[-1])
    e_hyb = abs(np.real(r['ru'][-1, j]) - R[-1])
    assert e_hyb < 1e-9 and e_hyb < 0.01 * e_awe
