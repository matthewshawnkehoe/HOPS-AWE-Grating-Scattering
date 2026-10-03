"""test_single_eps_delta.py -- Python port of test_single_eps_delta.m.

Checks every stage of the HOPS/AWE pipeline (interface data, fields u/w, DNOs G/J,
two-layer solution U/W, traces ubar/wbar) at a single (eps, delta), for all
truncation orders 0<=n<=N, 0<=m<=M, with Taylor, Pade and Pade-safe summation,
and draws the same plot_errors figures (1-6 and 11) as the MATLAB script.
"""
import os
import time

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from hops import (cheb, setup_2d, csqrt, setup_xi_u_nu_u_n_m, setup_xi_w_nu_w_n_m,
                  field_tfe_helmholtz_m_and_n, field_tfe_helmholtz_m_and_n_lf,
                  dno_tfe_helmholtz_m_and_n, dno_tfe_helmholtz_m_and_n_lf,
                  two_layer_solve, two_layer_solve_fast, fcn_sum, vol_fcn_sum, plot_errors)

DoTwoLayerTest = 3      # 3 = two_layer_solve_coupled (default, fastest), 1 = two_layer_solve (slow),
                        # 2 = two_layer_solve_fast (= MATLAB), else dummy
RunNumber = 1
SHOW = True
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'test_single_eps_delta')


def run(RunNumber=RunNumber, DoTwoLayerTest=DoTwoLayerTest, make_plots=True, verbose=True):
    if RunNumber == 1:
        M, Nx, Eps, sigma = 8, 16, 1e-7, 1e-2
    elif RunNumber == 2:
        M, Nx, Eps, sigma = 6, 24, 0.1, 0.1
    else:
        M, Nx, Eps, sigma = 12, 32, 0.1, 0.5
    N = M + 1
    Nz = Nx
    q = 1
    alpha_bar, d, c_0, n_u, n_w, Mode = 0.0, 2 * np.pi, 1.0, 1.0, 1.1, 2
    A_u, A_w, r = 3.0, 5.0, 2
    a, b = 1.0, 1.0
    identy = np.eye(Nz + 1)
    Dz, _ = cheb(Nz)
    xx = (d / Nx) * np.arange(Nx)
    f, f_x = np.cos(xx), -np.sin(xx)

    delta = sigma / (2 * q + 1)
    omega_bar = q + 0.5
    k_u_bar = n_u * omega_bar / c_0
    gamma_u_bar = csqrt(k_u_bar ** 2 - alpha_bar ** 2)
    xx, pp, abp, gubp, _, _ = setup_2d(Nx, d, alpha_bar, gamma_u_bar)
    k_w_bar = n_w * omega_bar / c_0
    gamma_w_bar = csqrt(k_w_bar ** 2 - alpha_bar ** 2)
    xx, pp, abp, gwbp, _, _ = setup_2d(Nx, d, alpha_bar, gamma_w_bar)
    pp_r, alpha_bar_r = pp[r], abp[r]

    xi_u, nu_u = setup_xi_u_nu_u_n_m(A_u, r, xx, pp, abp, gubp, f, f_x, Nx, N, M)
    u_n_m = field_tfe_helmholtz_m_and_n(xi_u, f, pp, gubp, alpha_bar, gamma_u_bar, Dz, a, Nx, Nz, N, M, identy, abp)
    G_n_m = dno_tfe_helmholtz_m_and_n(u_n_m, f, pp, Dz, a, Nx, Nz, N, M)
    xi_w, nu_w = setup_xi_w_nu_w_n_m(A_w, r, xx, pp, abp, gwbp, f, f_x, Nx, N, M)
    w_n_m = field_tfe_helmholtz_m_and_n_lf(xi_w, f, pp, gwbp, alpha_bar, gamma_w_bar, Dz, b, Nx, Nz, N, M, identy, abp)
    J_n_m = dno_tfe_helmholtz_m_and_n_lf(w_n_m, f, pp, Dz, b, Nx, Nz, N, M)
    tau2 = 1.0 if Mode == 1 else (n_u / n_w) ** 2
    zeta = xi_u - xi_w
    psi = -nu_u - tau2 * nu_w
    # physical interface data for oblique incidence (hops/config.py item 3)
    from hops.two_layer import phase_correction
    for _n in range(N + 1):
        for _m in range(M + 1):
            psi[:, _m, _n] = psi[:, _m, _n] - phase_correction(xi_u, xi_w, _m, _n, alpha_bar, f_x, tau2)

    t0 = time.time()
    args = (tau2, zeta, psi, gubp, gwbp, N, Nx, f, f_x, pp, alpha_bar, gamma_u_bar, gamma_w_bar,
            Dz, a, b, Nz, M, identy, abp)
    if DoTwoLayerTest == 1:
        U_n_m, W_n_m, ubar_n_m, wbar_n_m = two_layer_solve(*args)
    elif DoTwoLayerTest == 2:
        U_n_m, W_n_m, ubar_n_m, wbar_n_m = two_layer_solve_fast(*args)
    elif DoTwoLayerTest == 3:
        from hops.coupled import two_layer_solve_coupled
        U_n_m, W_n_m, ubar_n_m, wbar_n_m = two_layer_solve_coupled(*args)
    else:
        U_n_m = np.ones((Nx, M + 1, N + 1), dtype=complex)
        W_n_m = ubar_n_m = wbar_n_m = U_n_m
    if verbose:
        print(f'TwoLayerSolve time: {time.time() - t0:.2f} s')

    alpha_r = alpha_bar_r + delta * alpha_bar
    k_u, k_w = (1 + delta) * k_u_bar, (1 + delta) * k_w_bar
    g_ur, g_wr = csqrt(k_u ** 2 - alpha_r ** 2), csqrt(k_w ** 2 - alpha_r ** 2)
    ex = {}
    ex['xi_u'] = A_u * np.exp(1j * pp_r * xx) * np.exp(1j * g_ur * Eps * f)
    ex['nu_u'] = (-1j * g_ur + 1j * pp_r * Eps * f_x) * ex['xi_u']
    ex['xi_w'] = A_w * np.exp(1j * pp_r * xx) * np.exp(-1j * g_wr * Eps * f)
    ex['nu_w'] = (-1j * g_wr - 1j * pp_r * Eps * f_x) * ex['xi_w']
    ex['zeta'] = ex['xi_u'] - ex['xi_w']
    ex['psi'] = -ex['nu_u'] - tau2 * ex['nu_w']
    ex['ubar'] = A_u * np.exp(1j * pp_r * xx) * np.exp(1j * g_ur * a)
    ex['wbar'] = A_w * np.exp(1j * pp_r * xx) * np.exp(-1j * g_wr * (-b))
    # exact fields in transformed coordinates
    tz = np.cos(np.pi * np.arange(Nz + 1) / Nz)
    zp = (a / 2.0) * (tz - 1.0) + a
    Z = (a - Eps * f[:, None]) * zp[None, :] / a + Eps * f[:, None]
    u_ex = A_u * np.exp(1j * pp_r * xx)[:, None] * np.exp(1j * g_ur * Z)
    zp = (b / 2.0) * (tz - 1.0)
    Z = (b + Eps * f[:, None]) * zp[None, :] / b + Eps * f[:, None]
    w_ex = A_w * np.exp(1j * pp_r * xx)[:, None] * np.exp(-1j * g_wr * Z)

    coefs = dict(xi_u=xi_u, nu_u=nu_u, xi_w=xi_w, nu_w=nu_w, zeta=zeta, psi=psi,
                 G=G_n_m, J=J_n_m, U=U_n_m, W=W_n_m, ubar=ubar_n_m, wbar=wbar_n_m)
    target = dict(xi_u='xi_u', nu_u='nu_u', xi_w='xi_w', nu_w='nu_w', zeta='zeta', psi='psi',
                  G='nu_u', J='nu_w', U='xi_u', W='xi_w', ubar='ubar', wbar='wbar')
    err = {(k, st): np.zeros((N + 1, M + 1)) for k in list(coefs) + ['u', 'w'] for st in (1, 2, 3)}
    t0 = time.time()
    for n in range(N + 1):
        for m in range(M + 1):
            for st in (1, 2, 3):
                for k, cf in coefs.items():
                    err[(k, st)][n, m] = np.max(np.abs(ex[target[k]] - fcn_sum(st, cf, Eps, delta, Nx, n, m)))
                err[('u', st)][n, m] = np.max(np.abs(u_ex - vol_fcn_sum(st, u_n_m, Eps, delta, Nx, Nz, n, m)))
                err[('w', st)][n, m] = np.max(np.abs(w_ex - vol_fcn_sum(st, w_n_m, Eps, delta, Nx, Nz, n, m)))
    if verbose:
        print(f'Summation time: {time.time() - t0:.2f} s')

    if make_plots:
        os.makedirs(OUTDIR, exist_ok=True)
        E = lambda k, st: err[(k, st)]
        figs = [
            (1, 'Taylor', [E(k, 1) for k in ('xi_u', 'nu_u', 'zeta', 'xi_w', 'nu_w', 'psi')], True),
            (3, 'Pade', [E(k, 2) for k in ('xi_u', 'nu_u', 'zeta', 'xi_w', 'nu_w', 'psi')], True),
            (5, 'Pade', [E(k, 3) for k in ('xi_u', 'nu_u', 'zeta', 'xi_w', 'nu_w', 'psi')], True),
            (2, 'Taylor', [E(k, 1) for k in ('G', 'U', 'ubar', 'J', 'W', 'wbar')], False),
            (4, 'Pade', [E(k, 2) for k in ('G', 'U', 'ubar', 'J', 'W', 'wbar')], False),
            (6, 'Pade Safe', [E(k, 3) for k in ('G', 'U', 'ubar', 'J', 'W', 'wbar')], False),
            (11, 'Fields', [E('u', 1), E('u', 2), E('u', 3), E('w', 1), E('w', 2), E('w', 3)], None)]
        for num, lab, arrs, data in figs:
            # plot_errors.m always labels the panels G,U,ubar,J,W,wbar
            fig = plot_errors(num, lab, N, M, *arrs)
            fig.savefig(os.path.join(OUTDIR, f'test_single_eps_delta_fig{num}.png'), dpi=110)
    return err, dict(N=N, M=M, Eps=Eps, delta=delta)


def test_single_eps_delta():
    """Entry point when PyCharm/pytest runs this file as a test (files named test_*.py are
    collected by pytest).  Runs the full check, saves the 7 plot_errors figures to
    figures/test_single_eps_delta/, prints the error table and asserts the manufactured solution is recovered."""
    matplotlib.use('Agg')
    err, info = run(make_plots=True, verbose=True)
    _print_table(err, info)
    for k in ('xi_u', 'nu_u', 'zeta', 'psi', 'G', 'J', 'U', 'W', 'ubar', 'wbar', 'u', 'w'):
        for st in (1, 2, 3):
            assert err[(k, st)][-1, -1] < 1e-9, (k, st, err[(k, st)][-1, -1])


def _print_table(err, info):
    print(f"\nErrors at N={info['N']}, M={info['M']}, eps={info['Eps']:g}, delta={info['delta']:g}:")
    for k in ('xi_u', 'nu_u', 'zeta', 'psi', 'G', 'J', 'U', 'W', 'ubar', 'wbar', 'u', 'w'):
        print(f"  {k:5s}  Taylor {err[(k, 1)][-1, -1]:.2e}   Pade {err[(k, 2)][-1, -1]:.2e}"
              f"   Pade-safe {err[(k, 3)][-1, -1]:.2e}")


if __name__ == '__main__':
    import sys
    if '--no-show' in sys.argv:
        matplotlib.use('Agg')
    err, info = run()
    _print_table(err, info)
    if SHOW and '--no-show' not in sys.argv:
        plt.show()
