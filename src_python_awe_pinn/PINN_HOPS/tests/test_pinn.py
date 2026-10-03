"""Tests for the PINN grating solver and its comparison with HOPS/AWE.

  P1  the loss is exactly the governing equations (6a)-(6h): the exact flat-interface (Fresnel) solution
      has zero loss
  P2  ... and so does the HOPS/AWE solution of the corrugated problem (dielectric, metal, oblique incidence):
      PINN and HOPS solve the SAME boundary value problem (loss of the HOPS field ~1e-24)
  P3  the Taylor-mode (forward) derivatives of the networks equal reverse-mode autograd
  P4  the R, T, D post-processing equals the Fresnel formulas for a flat interface (D = 0)
  P5  the flat-interface ansatz reproduces eps = 0 exactly (zero loss for any network weights)
  P6  a short PINN training on the paper's Fig. 9 case (n_w = 1.1, cos x, eps = 0.1, omega = 1.5)
      reaches R within 2 % and |D| < 5e-3 of HOPS (a full run gets ~1e-3 and ~2e-5, see README)

    pytest -q tests                  (P6 takes ~2 min on 2 CPU cores; skip it with -k "not training")
"""
import os
import sys
import warnings

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
warnings.filterwarnings('ignore')

from pinn_hops import Grating2D, GratingPINN                              # noqa: E402
from pinn_hops.pinn import FieldNet, _grads                               # noqa: E402
from pinn_hops.hops_reference import hops_point, HOPSFieldTorch            # noqa: E402


def _flat_exact(G):
    r, t = G.flat_solution()
    gu, gw = complex(G.gamma_u), complex(G.gamma_w)

    def field(layer, x, z, e, o):
        v = (r * torch.exp(1j * gu * z) if layer == 'u' else t * torch.exp(-1j * gw * z)) + 0 * x
        return v.real, v.imag
    return field


def test_loss_zero_for_exact_flat_solution():
    for G in (Grating2D(eps=0.0), Grating2D(eps=0.0, n_w=0.05 + 2.275j), Grating2D(eps=0.0, alpha=0.3)):
        P = GratingPINN(G)
        P.override = _flat_exact(G)
        P.sample()
        L, _ = P.loss()
        assert float(L) < 1e-28, float(L)


def test_loss_zero_for_hops_solution():
    for G, K in ((Grating2D(eps=0.1), 6), (Grating2D(eps=0.1, n_w=0.05 + 2.275j), 6),
                 (Grating2D(eps=0.2, alpha=0.1), 6), (Grating2D(eps=0.1, n_w=1.48 + 1.883j, profile='cos4x'), 12)):
        Nx = 128 if G.profile == 'cos4x' else 32
        ref = hops_point(G, summation='taylor', Nx=Nx, Nz=48)
        P = GratingPINN(G, K=K)
        P.override = HOPSFieldTorch(ref, G)
        P.sample()
        L, parts = P.loss()
        assert float(L) < 1e-12, (G, parts)


def test_forward_mode_derivatives():
    net = FieldNet(K=5, width=16, depth=3, param_ranges=[(0, 0.2), (1, 2)])
    x = (torch.rand(40) * 6).requires_grad_()
    z = (torch.rand(40) * 2 - 1).requires_grad_()
    e, o = torch.rand(40) * 0.2, 1 + torch.rand(40)
    D = net.derivs(x, z, (e, o))
    fr, fi = net(x, z, (e, o))
    A = _grads(fr, fi, x, z)
    for k in ('x', 'z', 'xx', 'zz'):
        assert torch.max(torch.abs(torch.complex(D[k][:, 0], D[k][:, 1]) - A[k])) < 1e-12


def test_energy_formula_flat():
    for G in (Grating2D(eps=0.0), Grating2D(eps=0.0, n_w=1.5, alpha=0.4, mode='TE')):
        r, t = G.flat_solution()
        x = 2 * np.pi * np.arange(64) / 64
        R, T, D = G.energy(r * np.exp(1j * G.gamma_u * G.a) + 0 * x, t * np.exp(1j * G.gamma_w * G.b) + 0 * x)
        assert abs(R - abs(r) ** 2) < 1e-14 and abs(D) < 1e-14


def test_flat_ansatz_exact_at_eps0():
    P = GratingPINN(Grating2D(eps=0.0), K=4, width=16, depth=2)
    P.sample()
    L, _ = P.loss()
    assert float(L) < 1e-28


def test_short_training_dielectric():
    G = Grating2D(eps=0.1)
    ref = hops_point(G, summation='pade')
    P = GratingPINN(G, K=6, width=48, depth=4, seed=0)
    P.train(adam_iters=1500, lbfgs_iters=500, verbose=False)
    R, T, D = P.energy()
    assert abs(R - ref['R']) / ref['R'] < 2e-2, (R, ref['R'])
    assert abs(D) < 5e-3, D


# ---------------------------------------------------------------------------------------------
# improved variants (see README "Better PINN variants")
# ---------------------------------------------------------------------------------------------
def test_activation_derivatives():
    """P7  Taylor-mode derivatives are exact for every activation (tanh, sin, adaptive LAAF)"""
    for act in ('tanh', 'sin', 'laaf'):
        net = FieldNet(K=4, width=12, depth=3, activation=act)
        x = (torch.rand(30) * 6).requires_grad_()
        z = (torch.rand(30) * 2 - 1).requires_grad_()
        D = net.derivs(x, z)
        A = _grads(*net(x, z), x, z)
        for k in ('x', 'z', 'xx', 'zz'):
            assert torch.max(torch.abs(torch.complex(D[k][:, 0], D[k][:, 1]) - A[k])) < 1e-12, (act, k)


def test_lsq_features_derivatives():
    """P8  analytic derivatives of the least-squares PINN features match finite differences"""
    from pinn_hops.lsq_pinn import Features
    for basis in ('rfm', 'separable'):
        for act in ('tanh', 'sin', 'gauss'):
            F = Features(40, 3, act, basis)
            x, z, h = np.array([0.3, 2.0]), np.array([-0.4, 0.7]), 1e-4
            d = F(x, z)
            fd = lambda dx, dz: F(x + dx, z + dz, 0)['f']
            assert np.abs((fd(h, 0) - fd(-h, 0)) / (2 * h) - d['x']).max() < 1e-6
            assert np.abs((fd(0, h) - fd(0, -h)) / (2 * h) - d['z']).max() < 1e-6
            assert np.abs((fd(h, 0) - 2 * d['f'] + fd(-h, 0)) / h ** 2 - d['xx']).max() < 1e-4
            assert np.abs((fd(0, h) - 2 * d['f'] + fd(0, -h)) / h ** 2 - d['zz']).max() < 1e-4


def test_lsq_pinn_flat_and_dielectric():
    """P9  least-squares interface PINN: exact at eps = 0, and matches HOPS/AWE to ~1e-12 for the Fig. 9 case"""
    from pinn_hops.lsq_pinn import LSQPinn
    G0 = Grating2D(eps=0.0)
    R, T, D = LSQPinn(G0, K=4, n_feat=9 * 8, basis='separable').solve().energy()
    assert abs(R - abs(G0.flat_solution()[0]) ** 2) < 1e-14 and abs(D) < 1e-14
    G = Grating2D(eps=0.1)
    ref = hops_point(G, summation='pade', fields=False)
    P = LSQPinn(G, K=10, n_feat=21 * 24, basis='separable', Ms=32, rz=3.0).solve()
    R, T, D = P.energy()
    assert abs(R - ref['R']) / ref['R'] < 1e-10, (R, ref['R'])
    assert abs(D) < 1e-12, D
    # same with the features in the flattened (TFE) coordinates, sin activation
    P = LSQPinn(G, K=10, n_feat=21 * 12, basis='separable', Ms=32, activation='sin', coords='tfe').solve()
    R, T, D = P.energy()
    assert abs(R - ref['R']) / ref['R'] < 1e-9 and abs(D) < 1e-11, (R, D)


def test_tfe_feature_derivatives():
    """P10 chain rule of the flattened-coordinate features (both layers, steep cos 4x interface)"""
    from pinn_hops.lsq_pinn import Features, MappedFeatures
    G = Grating2D(eps=0.2, profile='cos4x')
    for layer, zz in (('u', np.array([0.5, 0.7])), ('w', np.array([-0.5, -0.7]))):
        M = MappedFeatures(Features(40, 3, 'tanh', 'separable', (0.0, 1.0)), layer, G)
        x, h = np.array([0.3, 2.0]), 1e-4
        d = M(x, zz)
        f = lambda a, b: M(x + a, zz + b, 0)['f']
        sc = np.abs(d['f']).max()
        assert np.abs((f(h, 0) - f(-h, 0)) / (2 * h) - d['x']).max() < 1e-5 * sc * 100
        assert np.abs((f(0, h) - f(0, -h)) / (2 * h) - d['z']).max() < 1e-5 * sc * 100
        assert np.abs((f(h, 0) - 2 * d['f'] + f(-h, 0)) / h ** 2 - d['xx']).max() < 1e-4 * sc * 100
        assert np.abs((f(0, h) - 2 * d['f'] + f(0, -h)) / h ** 2 - d['zz']).max() < 1e-4 * sc * 100
