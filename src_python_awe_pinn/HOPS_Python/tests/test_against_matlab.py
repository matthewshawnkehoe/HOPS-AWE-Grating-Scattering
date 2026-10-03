"""Validation suite for the Python HOPS/AWE port.

Group A compares against reference data produced by running the ORIGINAL MATLAB
sources (src.zip) under GNU Octave 8.4 (scripts in reference_data/octave/).
Group B checks the mathematics directly (manufactured solutions, exactness of the
Chebyshev derivative, energy conservation, fast vs. slow solvers, Pade robustness).

Run with   pytest -q tests     or   python test_mms_error.py
"""
import os
import sys
import warnings

import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
REF = os.path.join(ROOT, 'reference_data')

from hops import *                                       # noqa: E402,F401,F403
from hops import bvp                                     # noqa: E402


def _load(name):
    path = os.path.join(REF, name)
    if not os.path.exists(path):
        import pytest
        pytest.skip(f'{name} not found')
    return sio.loadmat(path)


def _rel(a, b):
    a, b = np.asarray(a), np.asarray(b)
    b = b.reshape(a.shape)
    return np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300)


# ---------------------------------------------------------------- Group A
def test_units_vs_matlab():
    R = sio.loadmat(os.path.join(REF, 'ref_units.mat'), squeeze_me=True)
    D, x = cheb(8)
    assert _rel(D, R['D8']) < 1e-13 and _rel(x, R['x8']) < 1e-15
    xx, kk, ap, bp, eep, _ = setup_2d(16, 2 * np.pi, 0.1, 1.3)
    assert _rel(bp, R['betap']) < 1e-14 and _rel(eep, R['eep']) < 1e-14
    _, _, _, bp2, _, _ = setup_2d(16, 2 * np.pi, 0, csqrt((0.05 + 2.275j) ** 2 * 2.25))
    assert _rel(bp2, R['betap2']) < 1e-14
    assert _rel(dx(R['u'], R['p']), R['ux']) < 1e-14
    assert _rel(dz(R['u'], D, 1.7), R['uz']) < 1e-14
    gq = gamma_exp(0.1, ap[2], 1.3, bp[2], np.sqrt(1.7), 6)
    assert _rel(gq, R['gq']) < 1e-14
    f, fx = np.cos(xx), -np.sin(xx)
    assert _rel(E_exp(gq, f, 5, 6), R['E']) < 1e-14
    assert _rel(E_exp_lf(gq, f, 5, 6), R['Elf']) < 1e-14
    assert _rel(A_exp(gq, fx, 5, 6), R['Aa']) < 1e-14
    assert _rel(T_dno(0.1, ap, 1.3, bp, ap[0] ** 2 + bp[0] ** 2, 16, 6), R['T']) < 1e-13
    U, W = AInverse(R['Q'], R['R'], bp, bp2, 16, 0.7)
    assert _rel(U, R['Ui']) < 1e-14 and _rel(W, R['Wi']) < 1e-14
    C = R['C']
    assert _rel(taylorsum_2_coeff(C, 0.1, 0.05, 8, 8), R['ts1']) < 1e-14
    assert _rel(taylorsum2(C, 0.1, 0.05, 8, 8), R['ts2']) < 1e-14
    assert _rel(padesum2(C, 0.1, 0.05, 8, 8), R['ps1']) < 1e-13
    assert _rel(padesum2_safe(C, 0.1, 0.05, 8, 8), R['ps2']) < 1e-13
    for st, key in ((1, 'ts1'), (2, 'ps1'), (3, 'ps2')):     # vectorised engine
        assert _rel(sum_series(st, C, 0.1, 0.05, 8, 8), R[key]) < 1e-13
    p, a, b = padesum(R['c1'], 0.7, 5)
    assert _rel(p, R['pp1']) < 1e-14 and _rel(a, R['pa1']) < 1e-13
    assert _rel(padesum_safe(R['c2'], 0.3, 2)[0], R['pp2']) < 1e-14
    r, a, b, mu, nu, _, _ = padeapprox(R['c1'], 4, 4)
    assert _rel(r(0.7), R['pr']) < 1e-14 and (mu, nu) == (R['mu'], R['nu'])


def _single_setup():
    M, Nx, N = 8, 16, 9
    Nz = Nx
    d, n_u, n_w, r, a, b = 2 * np.pi, 1.0, 1.1, 2, 1.0, 1.0
    Dz, _ = cheb(Nz)
    identy = np.eye(Nz + 1)
    xx = (d / Nx) * np.arange(Nx)
    f, f_x = np.cos(xx), -np.sin(xx)
    gub = 1.5
    xx, pp, abp, gubp, _, _ = setup_2d(Nx, d, 0.0, gub)
    gwb = 1.1 * 1.5
    xx, pp, abp, gwbp, _, _ = setup_2d(Nx, d, 0.0, gwb)
    return locals()


def test_single_pipeline_vs_matlab():
    R = _load('ref_single.mat')
    s = _single_setup()
    g = lambda k: s[k]
    Nx, Nz, N, M, Dz, identy = g('Nx'), g('Nz'), g('N'), g('M'), g('Dz'), g('identy')
    xx, pp, abp, gubp, gwbp, f, f_x = g('xx'), g('pp'), g('abp'), g('gubp'), g('gwbp'), g('f'), g('f_x')
    xu, nu_u = setup_xi_u_nu_u_n_m(3.0, 2, xx, pp, abp, gubp, f, f_x, Nx, N, M)
    xw, nu_w = setup_xi_w_nu_w_n_m(5.0, 2, xx, pp, abp, gwbp, f, f_x, Nx, N, M)
    assert _rel(xu, R['xi_u_r_n_m']) < 1e-13 and _rel(nu_w, R['nu_w_r_n_m']) < 1e-13
    u = field_tfe_helmholtz_m_and_n(xu, f, pp, gubp, 0.0, 1.5, Dz, 1.0, Nx, Nz, N, M, identy, abp)
    w = field_tfe_helmholtz_m_and_n_lf(xw, f, pp, gwbp, 0.0, 1.65, Dz, 1.0, Nx, Nz, N, M, identy, abp)
    G = dno_tfe_helmholtz_m_and_n(u, f, pp, Dz, 1.0, Nx, Nz, N, M)
    J = dno_tfe_helmholtz_m_and_n_lf(w, f, pp, Dz, 1.0, Nx, Nz, N, M)
    for mine, key in ((u, 'u_n_m'), (w, 'w_n_m'), (G, 'G_n_m'), (J, 'J_n_m')):
        assert _rel(mine, R[key]) < 1e-11, key
    tau2 = 1 / 1.21
    args = (tau2, xu - xw, -nu_u - tau2 * nu_w, gubp, gwbp, N, Nx, f, f_x, pp, 0.0, 1.5, 1.65,
            Dz, 1.0, 1.0, Nz, M, identy, abp)
    U, W, ub, wb = two_layer_solve_fast(*args)
    for mine, key in ((U, 'U_n_m'), (W, 'W_n_m'), (ub, 'ubar_n_m'), (wb, 'wbar_n_m')):
        assert _rel(mine, R[key]) < 1e-10, key


def test_single_errors_vs_matlab():
    """Full test_single_eps_delta.m error tables (n,m,SumType) agree with MATLAB."""
    R = _load('ref_single.mat')
    import test_single_eps_delta as t
    err, _ = t.run(make_plots=False, verbose=False)
    for mine, key in (('U', 'err_U'), ('W', 'err_W'), ('G', 'err_G'), ('J', 'err_J'),
                      ('ubar', 'err_ubar'), ('u', 'err_u')):
        for st in (1, 2, 3):
            a = np.log10(err[(mine, st)] + 1e-300)
            b = np.log10(R[key][:, :, st - 1] + 1e-300)
            big = b > -10            # above round-off: must agree closely
            assert np.max(np.abs(a[big] - b[big])) < 1e-3, (mine, st)
            assert np.all(a[~big] < -9.5), (mine, st)    # both at round-off level


def test_mms_vs_matlab():
    """Reduced mms_error.m (same physics, 6x7 (eps,delta) samples)."""
    R = _load('ref_mms.mat')
    import mms_error as me
    out, _ = me.run(dict(N_Eps=6, N_delta=7), verbose=False)
    o = out[0]
    # low orders agree to round-off; very high orders are round-off noise in both codes.
    # ubar/wbar are traces at z = +-4 of an evanescent mode (|ubar| ~ 5e^{-3.7*4} ~ 1e-6),
    # so they are compared on the scale of the O(1) field they are read from.
    scale = np.max(np.abs(R['U_n_m'][:, :9, :9]))
    for k in ('U', 'W', 'ubar', 'wbar'):
        a, b = o['coefs'][k], R[k + '_n_m']
        # (the operator solver forms ubar/wbar as sums of operator products, whose
        #  round-off is a few times larger than one field solve: allow 5e-8)
        tol = 1e-9 if k in ('U', 'W') else 5e-8
        assert np.max(np.abs(a[:, :9, :9] - b[:, :9, :9])) / scale < tol, k
    for nm, key in (('U', 'errU'), ('W', 'errW'), ('ubar', 'errubar')):
        for st in (1, 2, 3):
            a, b = np.log10(o['err'][(nm, st)]), np.log10(R[key][:, :, st - 1])
            big = b > -11
            assert np.max(np.abs(a[big] - b[big])) < 2e-2, (nm, st)   # 5% of an error ~1e-9
            assert np.all(a[~big] < -10.5)


def _refl_check(case, solver='operator'):
    """Reduced refl_map.m run (q = 1, 3; 4 x 5 (eps,delta) samples).

    Notes
    * The highest-order HOPS/AWE coefficients carry amplified round-off (they grow
      to ~1e16 for silver); two *Python* variants (LU vs cached inverse) differ from
      each other exactly as much as from MATLAB, so only orders n,m <= 5 are
      compared tightly.  Taylor sums (which use orders <= min(N,M)/2) agree to
      ~1e-9; Pade sums sit on nearly singular Toeplitz systems (MATLAB warns
      "matrix singular to machine precision", rcond ~ 1e-20) and agree to ~1e-2
      at worst, ~1e-8 typically.
    * For the metal the transmitted part rl (hence D) is not compared: MATLAB's
      '<' compares real parts of complex numbers (-> no propagating modes, rl = 0),
      whereas Octave compares by modulus.  The Python port follows MATLAB.
    """
    R = _load(f'ref_refl_{case}.mat')['out'][0, 0]
    import refl_map as rm
    res, info = rm.run(case, qq=(1, 3), N_Eps=4, N_delta=5, verbose=False, solver=solver, workers=1)
    metal = np.iscomplexobj(info['n_w']) and np.imag(info['n_w']) != 0
    for r in res:
        ref = R[f'q{r["q"]}'][0, 0]
        assert _rel(r['U_n_m'][:, :6, :6], ref['U'][:, :6, :6]) < 1e-7
        assert _rel(r['ru_flat'], ref['ru_flat']) < 1e-12
        ub = np.transpose(r['ubar_n_m'], (2, 1, 0))
        wb = np.transpose(r['wbar_n_m'], (2, 1, 0))
        tau2 = (1 / info['n_w']) ** 2
        omega_bar = r['q'] + 0.5
        gub, gwb = omega_bar, csqrt((info['n_w'] * omega_bar) ** 2)
        for st, e, u in ((1, 'ee_t', 'ru_t'), (2, 'ee_p', 'ru_p'), (3, 'ee_s', 'ru_s')):
            ee, ru, rl = energy_defect(tau2, ub, wb, 2 * np.pi, 0.0, gub, gwb, r['Eps'], r['delta'],
                                       info['Nx'], info['N'], info['M'], 4, 5, st)
            RRp, RRo = ru / r['ru_flat'], ref[u] / ref['ru_flat']
            dif = np.abs(RRp - RRo) / np.abs(RRo)
            if st == 1:
                assert np.max(dif) < 1e-9, (case, st, np.max(dif))
            else:
                assert np.median(dif) < 1e-6 and np.max(dif) < 2e-2, (case, st, np.max(dif))
            if not metal:
                assert np.max(np.abs(ee - ref[e])) < 1e-6 * max(1.0, np.max(np.abs(ref[e]))), (case, st)


def test_refl_silver_vs_matlab():
    _refl_check('silver')


def test_refl_silver_vs_matlab_classic_solver():
    _refl_check('silver', solver='fast')


def test_refl_dielectric_vs_matlab():
    _refl_check('dielectric')


def test_refl_gold_vs_matlab():
    """A9: paper Fig. 10b (gold, cos 4x, Pade) against the original MATLAB code (reduced grid);
    for a metal MATLAB has rl = 0, so D = 1 - R (the absorptance) and follows from R."""
    _refl_check('gold')
    import refl_map as rm
    res, _ = rm.run('gold', qq=(1, 3), N_Eps=4, N_delta=5, verbose=False, workers=1)
    for r in res:
        assert np.max(np.abs(r['ee'] - (1 - r['ru']))) < 1e-14 and np.all(r['rl'] == 0)


# ---------------------------------------------------------------- Group B
def test_dz_exact_derivative():
    """dz differentiates polynomials exactly -> no sign flip is needed."""
    Nz, a = 16, 2.5
    D, t = cheb(Nz)
    z = (a / 2) * (t - 1) + a                    # upper-layer map, z in [0, a]
    u = np.tile(z ** 5 - 3 * z ** 2, (4, 1))
    assert np.max(np.abs(dz(u, D, a) - (5 * z ** 4 - 6 * z))) < 1e-10
    b = 1.7
    z = (b / 2) * (t - 1)                        # lower-layer map, z in [-b, 0]
    u = np.tile(np.sin(z), (3, 1))
    assert np.max(np.abs(dz(u, D, b) - np.cos(z))) < 1e-12
    assert np.max(np.abs(dz(u, -D, b) - np.cos(z))) > 1.0   # the "fix" is wrong


def test_fast_equals_slow_two_layer():
    s = _single_setup()
    Nx, Nz, Dz, identy = s['Nx'], s['Nz'], s['Dz'], s['identy']
    N = M = 3
    xx, pp, abp, gubp, gwbp, f, f_x = (s[k] for k in ('xx', 'pp', 'abp', 'gubp', 'gwbp', 'f', 'f_x'))
    xu, nu_u = setup_xi_u_nu_u_n_m(3.0, 2, xx, pp, abp, gubp, f, f_x, Nx, N, M)
    xw, nu_w = setup_xi_w_nu_w_n_m(5.0, 2, xx, pp, abp, gwbp, f, f_x, Nx, N, M)
    tau2 = 1 / 1.21
    args = (tau2, xu - xw, -nu_u - tau2 * nu_w, gubp, gwbp, N, Nx, f, f_x, pp, 0.0, 1.5, 1.65,
            Dz, 1.0, 1.0, Nz, M, identy, abp)
    U1, W1, _, _ = two_layer_solve(*args)
    U2, W2, _, _ = two_layer_solve_fast(*args)
    assert _rel(U1, U2) < 1e-12 and _rel(W1, W2) < 1e-12


def test_mms_recovers_exact_solution():
    """With the correct sign conventions U, W, ubar are recovered to ~1e-13
    (a sign error in psi or dz gives O(1) - O(1e6) errors)."""
    import test_single_eps_delta as t
    err, _ = t.run(make_plots=False, verbose=False)
    for k in ('U', 'W', 'G', 'J', 'ubar', 'wbar', 'u', 'w'):
        for st in (1, 2, 3):
            assert err[(k, st)][-1, -1] < 1e-9, (k, st, err[(k, st)][-1, -1])


def test_energy_conservation_dielectric():
    """Real indices: R + T = 1 (energy defect ~ 0) for small eps."""
    import refl_map as rm
    res, info = rm.run('dielectric', qq=(1,), N_Eps=5, N_delta=5, RunNumber=1, verbose=False, workers=1)
    ee = res[0]['ee']
    assert np.max(np.abs(ee)) < 1e-7, np.max(np.abs(ee))


def test_pade_handles_tiny_and_zero_coefficients():
    """Robustness: Pade-safe and robust Pade give finite answers for degenerate data."""
    c = np.zeros((7, 7), dtype=complex)
    c[0, 0] = 2.0
    c[1, 0] = 1e-17                                  # (near-)zero higher coefficients
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for st in (2, 3, 4):
            v = sum_series(st, c, 1e-3, 1e-4, 6, 6)
            assert np.isfinite(v) and abs(v - 2.0) < 1e-12, (st, v)
    # geometric series 1/(1-rho): Pade is exact where Taylor is not
    cc = np.array([1.0 ** k for k in range(9)], dtype=complex)
    assert abs(padesum(cc, 0.9, 4)[0] - 10.0) < 1e-10
    assert abs(padesum_robust(cc, 0.9, 4) - 10.0) < 1e-10


def test_operator_solver_equals_classic():
    """two_layer_solve_operator (new, fast) == two_layer_solve_fast (MATLAB) to round-off."""
    from hops.operators import two_layer_solve_operator
    s = _single_setup()
    Nx, Nz, Dz, identy, N, M = s['Nx'], s['Nz'], s['Dz'], s['identy'], s['N'], s['M']
    xx, pp, abp, gubp, gwbp, f, f_x = (s[k] for k in ('xx', 'pp', 'abp', 'gubp', 'gwbp', 'f', 'f_x'))
    xu, nu_u = setup_xi_u_nu_u_n_m(3.0, 2, xx, pp, abp, gubp, f, f_x, Nx, N, M)
    xw, nu_w = setup_xi_w_nu_w_n_m(5.0, 2, xx, pp, abp, gwbp, f, f_x, Nx, N, M)
    tau2 = 1 / 1.21
    args = (tau2, xu - xw, -nu_u - tau2 * nu_w, gubp, gwbp, N, Nx, f, f_x, pp, 0.0, 1.5, 1.65,
            Dz, 1.0, 1.0, Nz, M, identy, abp)
    A = two_layer_solve_fast(*args)
    B = two_layer_solve_operator(*args)
    for a, b in zip(A, B):
        assert _rel(b, a) < 1e-10
    import hops.operators as op                       # FFT (large-Nx) code path too
    old = op.DENSE_X_MAX
    op.DENSE_X_MAX = 0
    try:
        B = two_layer_solve_operator(*args)
    finally:
        op.DENSE_X_MAX = old
    for a, b in zip(A, B):
        assert _rel(b, a) < 1e-10
    # nonzero alpha exercises the 2 i alpha d_x terms
    _, pp2, abp2, gubp2, _, _ = setup_2d(Nx, 2 * np.pi, 0.1, 1.5)
    _, _, _, gwbp2, _, _ = setup_2d(Nx, 2 * np.pi, 0.1, 1.65)
    xu, nu_u = setup_xi_u_nu_u_n_m(3.0, 2, xx, pp2, abp2, gubp2, f, f_x, Nx, 3, 3)
    xw, nu_w = setup_xi_w_nu_w_n_m(5.0, 2, xx, pp2, abp2, gwbp2, f, f_x, Nx, 3, 3)
    args = (tau2, xu - xw, -nu_u - tau2 * nu_w, gubp2, gwbp2, 3, Nx, f, f_x, pp2, 0.1, 1.5, 1.65,
            Dz, 1.0, 1.0, Nz, 3, identy, abp2)
    for a, b in zip(two_layer_solve_fast(*args), two_layer_solve_operator(*args)):
        assert _rel(b, a) < 1e-10


# ---------------------------------------------------------------- Group C (materials)
def test_refractiveindex_reader():
    """refractiveindex.info reader: Sellmeier formula (N-BK7 n_d = 1.5168, fused silica 1.4585)
    and tabulated n,k (Johnson & Christy silver/gold)."""
    from hops.materials import refractive_index as n
    assert abs(n('BK7', 0.5876).real - 1.5168) < 2e-4
    assert abs(n('SiO2', 0.5876).real - 1.4585) < 2e-4
    ag = n('Ag', 0.6199)            # J&C table entry: 0.055, 4.01? (interpolated around 0.62 um)
    assert 0.03 < ag.real < 0.08 and 3.8 < ag.imag < 4.3
    assert n('0.05+2.275i') == 0.05 + 2.275j and n('20i') == 20j


def test_surface_plasmon_position():
    """Grating-coupled SPP: HOPS/AWE reflectivity dip vs. the flat-interface SPP condition
    n_spp(lambda) * P / lambda = 1 (silver, P = 0.5 um, first order, band q = 0)."""
    import warnings
    import refl_map as rm
    from hops.materials import load
    warnings.simplefilter('ignore')
    P = 0.5
    res, info = rm.run('Ag_disp', qq=(0,), N_Eps=12, N_delta=60, verbose=False, workers=1, M=8,
                       max_delta=0.03, relative=False, eps_max=0.06)
    lam = np.concatenate([r['lam'] for r in res]) * P / (2 * np.pi)
    R = np.concatenate([np.real(r['ru'][-1]) for r in res])
    m = (lam > 0.505) & (lam < 0.6)
    lam_dip = lam[np.argmin(np.where(m, R, 9))]
    L = np.linspace(0.505, 0.6, 20001)
    eps = load('Ag')(L) ** 2
    lam_spp = L[np.argmin(abs(np.real(np.sqrt(eps / (eps + 1))) * P / L - 1))]
    assert abs(lam_dip - lam_spp) < 0.004, (lam_dip, lam_spp)


# ---------------------------------------------------------------- Group D (oblique incidence)
def _mms_layer(alpha_bar, Eps, lower, delta=0.01, N=8, M=8, Nx=32, Nz=32, r=2, L=1.0):
    d = 2 * np.pi
    Dz, _ = cheb(Nz)
    identy = np.eye(Nz + 1)
    xx = (d / Nx) * np.arange(Nx)
    f, f_x = np.cos(xx), -np.sin(xx)
    k_bar = 1.5 * (1.1 if lower else 1.0)
    gb = csqrt(k_bar ** 2 - alpha_bar ** 2)
    xx, pp, abp, gbp, _, _ = setup_2d(Nx, d, alpha_bar, gb)
    k = (1 + delta) * k_bar
    gr = csqrt(k ** 2 - (abp[r] + delta * alpha_bar) ** 2)
    if not lower:
        xi, nu = setup_xi_u_nu_u_n_m(3.0, r, xx, pp, abp, gbp, f, f_x, Nx, N, M)
        u = field_tfe_helmholtz_m_and_n(xi, f, pp, gbp, alpha_bar, gb, Dz, L, Nx, Nz, N, M, identy, abp)
        G = dno_tfe_helmholtz_m_and_n(u, f, pp, Dz, L, Nx, Nz, N, M)
        xi_ex = 3 * np.exp(1j * pp[r] * xx) * np.exp(1j * gr * Eps * f)
        nu_ex = (-1j * gr + 1j * pp[r] * Eps * f_x) * xi_ex
    else:
        xi, nu = setup_xi_w_nu_w_n_m(3.0, r, xx, pp, abp, gbp, f, f_x, Nx, N, M)
        u = field_tfe_helmholtz_m_and_n_lf(xi, f, pp, gbp, alpha_bar, gb, Dz, L, Nx, Nz, N, M, identy, abp)
        G = dno_tfe_helmholtz_m_and_n_lf(u, f, pp, Dz, L, Nx, Nz, N, M)
        xi_ex = 3 * np.exp(1j * pp[r] * xx) * np.exp(-1j * gr * Eps * f)
        nu_ex = (-1j * gr - 1j * pp[r] * Eps * f_x) * xi_ex
    return np.max(np.abs(nu_ex - fcn_sum(2, G, Eps, delta, Nx, N, M)))


def test_mms_oblique_incidence():
    """alpha != 0: with config.ALPHA_FIX the upper and lower TFE fields/DNOs reproduce the
    manufactured solution to round-off (the MATLAB formulation is off by 1e-5..1e-2)."""
    from hops import config
    assert config.ALPHA_FIX
    for lower in (False, True):
        for Eps in (0.0, 1e-3, 0.05):
            assert _mms_layer(0.1, Eps, lower) < 1e-10, (lower, Eps)


def test_energy_defect_oblique_dielectric():
    """Paper Fig. 14 setting (alpha = 0.01, n_w = 1.1): energy defect as small as for alpha = 0."""
    import warnings
    import refl_map as rm
    warnings.simplefilter('ignore')
    med = {}
    for al in (0.0, 0.01):
        res, _ = rm.run('dielectric_alpha', qq=(1,), N_Eps=20, N_delta=20, verbose=False, workers=1, alpha=al)
        med[al] = np.median(np.abs(res[0]['ee']))
    assert med[0.01] < 3 * med[0.0] and med[0.01] < 1e-7, med


def test_coupled_solver_equals_classic():
    """B7: two_layer_solve_coupled (one interleaved recursion per layer) == two_layer_solve_fast (MATLAB)."""
    from hops.coupled import two_layer_solve_coupled
    s = _single_setup()
    Nx, Nz, Dz, identy, N, M = s['Nx'], s['Nz'], s['Dz'], s['identy'], s['N'], s['M']
    xx, pp, abp, gubp, gwbp, f, f_x = (s[k] for k in ('xx', 'pp', 'abp', 'gubp', 'gwbp', 'f', 'f_x'))
    for alpha in (0.0, 0.1):
        _, pp2, abp2, gu2, _, _ = setup_2d(Nx, 2 * np.pi, alpha, 1.5)
        _, _, _, gw2, _, _ = setup_2d(Nx, 2 * np.pi, alpha, 1.65)
        xu, nu_u = setup_xi_u_nu_u_n_m(3.0, 2, xx, pp2, abp2, gu2, f, f_x, Nx, N, M)
        xw, nu_w = setup_xi_w_nu_w_n_m(5.0, 2, xx, pp2, abp2, gw2, f, f_x, Nx, N, M)
        tau2 = 1 / 1.21
        args = (tau2, xu - xw, -nu_u - tau2 * nu_w, gu2, gw2, N, Nx, f, f_x, pp2, alpha, 1.5, 1.65,
                Dz, 1.0, 1.0, Nz, M, identy, abp2)
        A = two_layer_solve_fast(*args)
        B = two_layer_solve_coupled(*args)
        for a, b in zip(A, B):
            assert _rel(b, a) < 1e-10, (alpha, _rel(b, a))


def test_refl_dielectric_full_map_vs_matlab():
    """A8: the FULL paper Fig. 9 map (6 bands, 100 x 100, N = M = 16, Taylor) against the original MATLAB
    code run under Octave (reference_data/ref_refl_dielectric_full.mat): R and the energy defect D agree
    to round-off everywhere, so any visual difference to the paper is rendering, not numerics."""
    path = os.path.join(REF, 'ref_refl_dielectric_full.mat')
    if not os.path.exists(path):
        import pytest
        pytest.skip('ref_refl_dielectric_full.mat not found')
    import refl_map as rm
    ref = sio.loadmat(path, squeeze_me=True, struct_as_record=False)['out']
    res, _ = rm.run('dielectric', verbose=False, workers=2)
    py = {r['q']: r for r in res}
    for k in ref._fieldnames:
        m = getattr(ref, k)
        r = py[int(k[1:])]
        assert np.max(np.abs(r['ee'] - m.ee_t)) < 1e-11
        assert np.max(np.abs(r['ru'] - m.ru_t)) < 1e-11
