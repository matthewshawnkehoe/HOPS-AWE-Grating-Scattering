"""Validation suite for the 3D HOPS/AWE code (package hops3d).

Group E  checks the 3D joint boundary/frequency perturbation method:
  E1  y-invariant gratings: the 3D solver reduces to the validated 2D solver (to round-off)
  E2  y-invariant energy defect / reflectivity == 2D energy_defect
  E3  the three 3D two-layer solvers (coupled, lean, operator) agree
  E4  manufactured doubly periodic solution recovered (crossed profile, oblique alpha, beta != 0)
  E5  volume fields u, w recovered (test_single_eps_delta_3D, RunNumber 1)
  E6  flat interface: Fresnel reflection / transmission (oblique, TE and TM)
  E7  energy conservation R + T = 1 for a crossed dielectric grating (Pade, |D| ~ 1e-12)
  E8  x <-> y symmetry for a symmetric profile and alpha = beta
  E9  frequency expansion of i gamma_pq(delta) (T_dno_3d) vs. direct evaluation
  E10 Rayleigh-frequency windows (3D lattice; Ny = 1 gives the 2D bands [q, q+1])
  E11 Fourier-domain and physical-domain summation of the reflectivity agree (identical for Taylor)
  E12 conical incidence on a y-invariant grating (beta != 0): energy conservation

Run with   pytest -q tests_3d     or   python test_mms_error_3D.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import hops                                              # noqa: E402
import hops3d as h3                                      # noqa: E402
from hops3d.mms import manufactured_data, exact_values   # noqa: E402
from hops3d.windows import make_windows, rayleigh_3d     # noqa: E402
from hops.summation import fcn_sum                       # noqa: E402


def _rel(a, b):
    return np.max(np.abs(np.asarray(a) - np.asarray(b))) / max(np.max(np.abs(b)), 1e-300)


def _solve_2d(Nx, Nz, N, M, n_u, n_w, omega, alpha, a=1.0, b=1.0):
    from hops.operators import two_layer_solve_operator
    d = 2 * np.pi
    gu = hops.csqrt((n_u * omega) ** 2 - alpha ** 2)
    gw = hops.csqrt((n_w * omega) ** 2 - alpha ** 2)
    xx, pp, abp, gup, _, _ = hops.setup_2d(Nx, d, alpha, gu)
    _, _, _, gwp, _, _ = hops.setup_2d(Nx, d, alpha, gw)
    f, fx = np.cos(xx), -np.sin(xx)
    z2, p2 = hops.setup_zeta_psi_n_m(xx, pp, alpha, gu, f, fx, Nx, N, M)
    Dz, _ = hops.cheb(Nz)
    tau2 = (n_u / n_w) ** 2
    out = two_layer_solve_operator(tau2, z2, p2, gup, gwp, N, Nx, f, fx, pp, alpha, gu, gw, Dz, a, b, Nz, M,
                                   np.eye(Nz + 1), abp)
    return out, (tau2, d, alpha, gu, gw)


def _solve_3d(Nx, Ny, Nz, N, M, n_u, n_w, omega, alpha, beta, profile, solver='coupled', **kw):
    P = h3.make_problem(Nx, Ny, Nz, N, M, n_u, n_w, omega, alpha, beta, profile=profile, **kw)
    z, p = h3.setup_zeta_psi_n_m_3d(alpha, beta, P.gamma_u_bar, P.f, P.f_x, P.f_y, N, M)
    return P, h3.SOLVERS_3D[solver](P, z, p)


# ---------------------------------------------------------------------------------- E
def test_y_invariant_reduces_to_2d():
    Nx, Nz, N, M = 16, 16, 6, 6
    for alpha in (0.0, 0.1):
        (U2, W2, ub2, wb2), _ = _solve_2d(Nx, Nz, N, M, 1.0, 1.1, 1.5, alpha)
        for Ny in (1, 4):
            P, (U, W, ub, wb) = _solve_3d(Nx, Ny, Nz, N, M, 1.0, 1.1, 1.5, alpha, 0.0, 'cosx')
            for A3, A2 in ((U, U2), (W, W2), (ub, ub2), (wb, wb2)):
                for j in range(Ny):
                    assert _rel(A3[:, j], A2) < 1e-11, (alpha, Ny, _rel(A3[:, j], A2))


def test_y_invariant_energy_equals_2d():
    Nx, Nz, N, M = 16, 16, 8, 8
    (U2, W2, ub2, wb2), (tau2, d, alpha, gu, gw) = _solve_2d(Nx, Nz, N, M, 1.0, 1.1, 1.5, 0.0)
    Eps = np.linspace(0, 0.1, 6)
    delta = np.linspace(-0.2, 0.2, 5)
    ub = np.transpose(ub2, (2, 1, 0)); wb = np.transpose(wb2, (2, 1, 0))
    ee2, ru2, rl2 = hops.energy_defect(tau2, ub, wb, d, alpha, gu, gw, Eps, delta, Nx, N, M, 6, 5, 1)
    P, (U, W, ubar, wbar) = _solve_3d(Nx, 1, Nz, N, M, 1.0, 1.1, 1.5, 0.0, 0.0, 'cosx')
    for dom in ('physical', 'fourier'):
        ee, ru, rl = h3.energy_defect_3d(P.tau2, ubar, wbar, P.kx, P.ky, 0.0, 0.0, P.gamma_u_bar, P.gamma_w_bar,
                                         Eps, delta, N, M, 1, sum_domain=dom)
        assert np.max(np.abs(ru - ru2)) < 1e-12 and np.max(np.abs(rl - rl2)) < 1e-12, dom


def test_solvers_agree():
    args = (8, 8, 12, 5, 5, 1.0, 1.3, 1.2, 0.07, 0.03, 'egg')
    _, ref = _solve_3d(*args, solver='coupled')
    for s in ('lean', 'operator'):
        _, out = _solve_3d(*args, solver=s)
        for A, B in zip(out, ref):
            assert _rel(A, B) < 1e-11, (s, _rel(A, B))


def test_mms_crossed_oblique():
    N = M = 8
    P = h3.make_problem(16, 16, 16, N, M, 1.0, 1.3, 1.5, 0.13, 0.21, profile='egg')
    A, B, rs = 3.0, 5.0, (1, 3)
    xi_u, nu_u, xi_w, nu_w, zeta, psi = manufactured_data(P, A, B, rs)
    Su, Sw = P.setups()
    G = h3.dno_tfe_helmholtz_3d(xi_u, Su, N, M)[0]
    J = h3.dno_tfe_helmholtz_3d(xi_w, Sw, N, M)[0]
    U, W, ub, wb = h3.two_layer_solve_3d_coupled(P, zeta, psi)
    eps, delta = 0.01, 0.01
    ex = exact_values(P, A, B, rs, eps, delta)
    K = P.Nx * P.Ny
    for k, c in dict(U=U, W=W, ubar=ub, wbar=wb, G=G, J=J, psi=psi, zeta=zeta).items():
        approx = fcn_sum(1, c.reshape(K, M + 1, N + 1), eps, delta, K, N, M, True).reshape(P.Nx, P.Ny)
        assert _rel(approx, ex[k]) < 1e-10, (k, _rel(approx, ex[k]))


def test_single_eps_delta_3d_fields():
    import test_single_eps_delta_3D as t
    err, info = t.run(make_plots=False, verbose=False)
    for k in t.KEYS:
        for st in (1, 2, 3):
            assert err[(k, st)][-1, -1] < 1e-11, (k, st, err[(k, st)][-1, -1])


def test_flat_interface_fresnel():
    for Mode in (1, 2):
        P, (U, W, ub, wb) = _solve_3d(8, 8, 16, 3, 3, 1.0, 1.5 + 0.1j, 1.3, 0.3, 0.4, 'cosxcosy', Mode=Mode)
        ee, ru, rl = h3.energy_defect_3d(P.tau2, ub, wb, P.kx, P.ky, 0.3, 0.4, P.gamma_u_bar, P.gamma_w_bar,
                                         np.array([0.0]), np.array([0.0]), 3, 3, 2)
        gu, gw = P.gamma_u_bar, P.gamma_w_bar
        r = (gu - P.tau2 * gw) / (gu + P.tau2 * gw)
        assert abs(ru[0, 0] - abs(r) ** 2) < 1e-13, (Mode, ru[0, 0], abs(r) ** 2)
        # lossless check of R + T = 1 with a real lower index
        P, (U, W, ub, wb) = _solve_3d(8, 8, 16, 3, 12, 1.0, 1.5, 1.3, 0.3, 0.4, 'cosxcosy', Mode=Mode)
        ee, ru, rl = h3.energy_defect_3d(P.tau2, ub, wb, P.kx, P.ky, 0.3, 0.4, P.gamma_u_bar, P.gamma_w_bar,
                                         np.array([0.0]), np.linspace(-0.1, 0.1, 5), 0, 12, 2)
        assert np.max(np.abs(ee)) < 1e-13, np.max(np.abs(ee))


def test_energy_conservation_crossed_dielectric():
    w = make_windows(1.0, 2.0, 1.0, 1.1, 0.0, 0.0, mode='joint')[0]          # [1, sqrt(2)/1.1]
    _, ob, dm, _ = w
    N = M = 12
    P, (U, W, ub, wb) = _solve_3d(32, 32, 16, N, M, 1.0, 1.1, ob, 0.01, 0.005, 'fs')
    Eps = np.linspace(0, 0.2, 21)
    delta = np.linspace(-dm, dm, 21)
    ee, ru, rl = h3.energy_defect_3d(P.tau2, ub, wb, P.kx, P.ky, 0.01, 0.005, P.gamma_u_bar, P.gamma_w_bar,
                                     Eps, delta, N, M, 2)
    assert np.max(np.abs(ee)) < 1e-11, np.max(np.abs(ee))


def test_xy_symmetry():
    P, (U, W, ub, wb) = _solve_3d(16, 16, 16, 6, 6, 1.0, 1.1, 1.2, 0.03, 0.03, 'egg')
    for A in (U, W, ub, wb):
        assert _rel(A.transpose(1, 0, 2, 3), A) < 1e-12


def test_T_dno_3d_expansion():
    Nx = Ny = 4
    alpha, beta, omega, n = 0.2, -0.3, 1.3, 1.7
    P = h3.make_problem(Nx, Ny, 8, 1, 12, 1.0, n, omega, alpha, beta, profile='cosxcosy')
    k2 = P.alphap[0] ** 2 + P.betap[0] ** 2 + P.gammap_w[0, 0] ** 2
    T = h3.T_dno_3d(alpha, beta, P.alphap, P.betap, P.gamma_w_bar, P.gammap_w, k2, 12)
    for dl in (0.01, -0.01):
        s = 1 + dl
        ap = alpha * s + P.kx[:, None]
        bq = beta * s + P.ky[None, :]
        exact = 1j * h3.outgoing_sqrt(k2 * s ** 2 - ap ** 2 - bq ** 2)
        approx = np.sum(T * dl ** np.arange(T.shape[-1]), axis=-1)
        ok = np.abs(P.gammap_w) > 0.8          # modes far from their Rayleigh point (radius of convergence)
        assert _rel(approx[ok], exact[ok]) < 1e-10, _rel(approx[ok], exact[ok])


def test_rayleigh_windows():
    w = rayleigh_3d(1.0, 0.0, 0.0, 0.5, 3.1)
    assert np.allclose(w, [1, np.sqrt(2), 2, np.sqrt(5), np.sqrt(8), 3])
    wins = make_windows(1.0, 7.0, 1.0, 1.1, 0.0, 0.0, Ny=1)
    assert [round(x[1], 12) for x in wins] == [1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
    assert np.allclose([x[2] for x in wins], [0.99 / (2 * q + 1) for q in range(1, 7)])
    # fixed-angle Rayleigh frequencies solve n w = |(a w + p, b w + q)|
    a_, b_ = 0.3, 0.2
    for wv in rayleigh_3d(1.0, 0, 0, 0.5, 3.0, angle=(a_, b_)):
        dmin = min(abs(wv - np.hypot(a_ * wv + p, b_ * wv + q)) for p in range(-5, 6) for q in range(-5, 6))
        assert dmin < 1e-10


def test_fourier_vs_physical_summation():
    N = M = 8
    _, ob, dm, _ = make_windows(1.0, 2.0, 1.0, 1.1, 0.0, 0.0, mode='joint')[0]
    P, (U, W, ub, wb) = _solve_3d(16, 16, 16, N, M, 1.0, 1.1, ob, 0.0, 0.0, 'cosxcosy')
    Eps, delta = np.linspace(0, 0.1, 6), np.linspace(-dm / 2, dm / 2, 5)
    args = (P.tau2, ub, wb, P.kx, P.ky, 0.0, 0.0, P.gamma_u_bar, P.gamma_w_bar, Eps, delta, N, M)
    for st, tol in ((1, 1e-14), (2, 1e-8)):
        r1 = h3.energy_defect_3d(*args, st, sum_domain='physical')
        r2 = h3.energy_defect_3d(*args, st, sum_domain='fourier')
        assert np.max(np.abs(r1[1] - r2[1])) < tol, (st, np.max(np.abs(r1[1] - r2[1])))


def test_conical_incidence_energy():
    # y-invariant grating, beta != 0 (conical mount): a genuinely 3D effect the 2D code cannot model
    wins = make_windows(1.0, 2.0, 1.0, 1.1, 0.0, 0.3, mode='joint', Ny=1)
    _, ob, dm, _ = wins[0]
    N = M = 10
    P, (U, W, ub, wb) = _solve_3d(32, 1, 24, N, M, 1.0, 1.1, ob, 0.0, 0.3, 'cosx')
    ee, ru, rl = h3.energy_defect_3d(P.tau2, ub, wb, P.kx, P.ky, 0.0, 0.3, P.gamma_u_bar, P.gamma_w_bar,
                                     np.linspace(0, 0.1, 6), np.linspace(-dm / 2, dm / 2, 5), N, M, 2)
    assert np.max(np.abs(ee)) < 1e-9, np.max(np.abs(ee))
