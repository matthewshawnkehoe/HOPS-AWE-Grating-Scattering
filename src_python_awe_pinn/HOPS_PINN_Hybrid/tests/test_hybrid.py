"""Tests for the HOPS/AWE + least-squares-PINN hybrid.

  H1  the basis derivatives (FFT in x, Chebyshev in zeta) and the chain rule through the TFE map are
      exact for a smooth test field (compared with the analytic derivatives)
  H2  at the band centre (delta = 0, small eps) the AWE Taylor sum already satisfies the equations to
      round-off (PINN loss ~1e-23): the hybrid's residual rows and HOPS discretise the same problem
  H3  the physics-informed sum is never worse than the Taylor sum in PINN loss, and at the band edge of
      the paper's Fig. 9 band it reduces the error in R from ~1e-5 to ~1e-14 (reference: HOPS at delta = 0)
  H4  the adaptive map equals the full hybrid map where it solves, and the Taylor sum elsewhere
  H5  POD compression of the basis does not change R
"""
import os
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings('ignore')

from hybrid import AWEBand, PISum                     # noqa: E402
from pinn_hops import Grating2D                       # noqa: E402
from pinn_hops.hops_reference import hops_point       # noqa: E402


@pytest.fixture(scope='module')
def band():
    return AWEBand('dielectric', q=1, N_Eps=5, N_delta=5)


def test_chain_rule(band):
    B = band
    S = PISum(B)
    X, Z = np.meshgrid(B.x, B.zeta, indexing='ij')
    # test field Phi(x, zeta) = exp(i 2x) (zeta^3 + zeta), mapped to physical coordinates of the upper layer
    Phi = np.exp(2j * X) * (Z ** 3 + Z)
    S.feat = {layer: B._derivs(Phi[:, :, None]) for layer in ('u', 'w')}
    eps = 0.15
    d = S._phys('u', eps)
    g = eps * B.f[:, None]
    gx = eps * B.fx[:, None]
    L = B.a - g
    z = g + Z * L
    # analytic: phi(x, z) = Phi(x, (z - g)/(a - g)); check d/dz and d2/dz2 exactly, d/dx by finite differences
    zeta = lambda x_, z_, g_: (z_ - g_) / (B.a - g_)
    assert np.allclose(d['z'][:, :, 0], np.exp(2j * X) * (3 * Z ** 2 + 1) / L, atol=1e-10)
    assert np.allclose(d['zz'][:, :, 0], np.exp(2j * X) * 6 * Z / L ** 2, atol=1e-9)
    h = 1e-6
    gp = eps * np.cos(B.x + h)[:, None]
    gm = eps * np.cos(B.x - h)[:, None]
    fp = np.exp(2j * (X + h)) * (zeta(X + h, z, gp) ** 3 + zeta(X + h, z, gp))
    fm = np.exp(2j * (X - h)) * (zeta(X - h, z, gm) ** 3 + zeta(X - h, z, gm))
    assert np.allclose(d['x'][:, :, 0], (fp - fm) / (2 * h), atol=1e-7)
    del gx


def test_taylor_sum_is_consistent(band):
    s = PISum(band).indicator(0.05, band.omega_bar)
    assert s['loss_awe_taylor'] < 1e-20, s


def test_band_edge_improvement(band):
    S = PISum(band, compress=1e-13)
    e, o = 0.2, band.omega_bar * (1 - 0.99 / 3)
    ref = hops_point(Grating2D(eps=e, omega=o, n_w=1.1), N=24, Nx=64, Nz=48, summation='pade', fields=False)
    s = S.solve(e, o)
    assert s['loss'] <= s['loss_awe_taylor']
    assert abs(s['R_taylor'] - ref['R']) > 1e-6          # the AWE Taylor sum is off at the band edge ...
    assert abs(s['R'] - ref['R']) < 1e-12                # ... the physics-informed sum is not
    assert abs(s['D']) < 1e-12


def test_adaptive_map(band):
    S = PISum(band, compress=1e-13)
    full = S.map()
    ad = S.map(adaptive_tol=1e-18)
    assert 0 < ad['n_solved'] < band.Eps.size * band.delta.size
    assert np.nanmax(np.abs(ad['R'] - full['R'])) < 1e-12


def test_compression(band):
    e, o = 0.2, band.omega_bar * 1.3
    R0 = PISum(band).solve(e, o)['R']
    R1 = PISum(band, compress=1e-13).solve(e, o)['R']
    assert abs(R0 - R1) < 1e-13
