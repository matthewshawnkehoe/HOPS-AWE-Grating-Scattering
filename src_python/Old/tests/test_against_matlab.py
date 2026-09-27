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
        assert np.max(np.abs(a[:, :9, :9] - b[:, :9, :9])) / scale < 1e-9, k
    for nm, key in (('U', 'errU'), ('W', 'errW'), ('ubar', 'errubar')):
        for st in (1, 2, 3):
            a, b = np.log10(o['err'][(nm, st)]), np.log10(R[key][:, :, st - 1])
            big = b > -11
            assert np.max(np.abs(a[big] - b[big])) < 1e-3, (nm, st)
            assert np.all(a[~big] < -10.5)


def _refl_check(case):
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
    res, info = rm.run(case, qq=(1, 3), N_Eps=4, N_delta=5, verbose=False)
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


def test_refl_dielectric_vs_matlab():
    _refl_check('dielectric')


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
    res, info = rm.run('dielectric', qq=(1,), N_Eps=5, N_delta=5, RunNumber=1, verbose=False)
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
