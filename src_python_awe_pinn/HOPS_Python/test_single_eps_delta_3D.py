"""test_single_eps_delta_3D.py -- 3D analogue of test_single_eps_delta.py / test_single_eps_delta.m.

Checks every stage of the 3D HOPS/AWE pipeline (interface data, volume fields u/w, DNOs G/J,
two-layer solution U/W, traces ubar/wbar) at a single (eps, delta) for all truncation orders
0 <= n <= N, 0 <= m <= M, with Taylor, Pade and Pade-safe summation, against the manufactured
doubly periodic solution  A exp(i p_r x + i q_s y +- i gamma_rs z), and draws the same
plot_errors figures (1-6 and 11) as the 2D script.

    python test_single_eps_delta_3D.py            # RunNumber 1 (as the MATLAB default)
    python test_single_eps_delta_3D.py --run 3 --alpha 0.1 --beta 0.2 --profile egg
    pytest test_single_eps_delta_3D.py            # PyCharm: right-click -> Run 'pytest ...'
"""
import argparse
import os
import time

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import hops3d as h3
from hops3d.mms import manufactured_data, exact_values
from hops3d.windows import make_windows
from hops.summation import fcn_sum, vol_fcn_sum
from hops.plotting import plot_errors

DoTwoLayerTest = 'coupled'   # 'coupled' | 'lean' (= two_layer_solve_fast ordering) | 'operator' | None (dummy)
RunNumber = 1
SHOW = True
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures_3d', 'test_single_eps_delta')


def run(RunNumber=RunNumber, DoTwoLayerTest=DoTwoLayerTest, make_plots=True, verbose=True,
        alpha=0.0, beta=0.0, profile='cosxcosy', rs=(2, 1)):
    if RunNumber == 1:
        M, Nx, Eps, sigma = 8, 16, 1e-7, 1e-2
    elif RunNumber == 2:
        M, Nx, Eps, sigma = 6, 16, 0.1, 0.1
    else:
        M, Nx, Eps, sigma = 10, 16, 0.1, 0.5
    N = M + 1
    Ny, Nz = Nx, Nx
    n_u, n_w, Mode = 1.0, 1.1, 2
    A_u, A_w = 3.0, 5.0
    a, b = 1.0, 1.0
    # first 3D frequency window (between the Rayleigh frequencies 1 and sqrt 2; 2D: omega_bar = 1.5)
    _, omega_bar, dmax, _ = make_windows(1.0, 2.0, n_u, n_w, alpha, beta)[0]
    delta = sigma * dmax / 0.99
    P = h3.make_problem(Nx, Ny, Nz, N, M, n_u, n_w, omega_bar, alpha, beta, profile=profile, a=a, b=b, Mode=Mode)
    xi_u, nu_u, xi_w, nu_w, zeta, psi = manufactured_data(P, A_u, A_w, rs)
    Su, Sw = P.setups()
    u_n_m, ou = h3.field_tfe_helmholtz_3d(xi_u, Su, N, M)
    w_n_m, ow = h3.field_tfe_helmholtz_3d(xi_w, Sw, N, M)
    std = lambda D: np.transpose(D[..., 0], (2, 3, 1, 0))
    G_n_m, J_n_m = std(ou['DNO']), std(ow['DNO'])

    t0 = time.time()
    if DoTwoLayerTest:
        U_n_m, W_n_m, ubar_n_m, wbar_n_m = h3.SOLVERS_3D[DoTwoLayerTest](P, zeta, psi)
    else:
        U_n_m = np.ones_like(zeta)
        W_n_m = ubar_n_m = wbar_n_m = U_n_m
    if verbose:
        print(f'TwoLayerSolve ({DoTwoLayerTest}) time: {time.time() - t0:.2f} s')

    ex = exact_values(P, A_u, A_w, rs, Eps, delta, volume=True)
    coefs = dict(xi_u=xi_u, nu_u=nu_u, xi_w=xi_w, nu_w=nu_w, zeta=zeta, psi=psi,
                 G=G_n_m, J=J_n_m, U=U_n_m, W=W_n_m, ubar=ubar_n_m, wbar=wbar_n_m)
    K = Nx * Ny
    flat = {k: v.reshape(K, M + 1, N + 1) for k, v in coefs.items()}
    exf = {k: v.reshape(-1) for k, v in ex.items() if k not in ('u', 'w')}
    uv = u_n_m.reshape(K, Nz + 1, M + 1, N + 1)
    wv = w_n_m.reshape(K, Nz + 1, M + 1, N + 1)
    uex, wex = ex['u'].reshape(K, Nz + 1), ex['w'].reshape(K, Nz + 1)
    err = {(k, st): np.zeros((N + 1, M + 1)) for k in list(coefs) + ['u', 'w'] for st in (1, 2, 3)}
    t0 = time.time()
    for n in range(N + 1):
        for m in range(M + 1):
            for st in (1, 2, 3):
                for k, cf in flat.items():
                    err[(k, st)][n, m] = np.max(np.abs(exf[k] - fcn_sum(st, cf, Eps, delta, K, n, m)))
                err[('u', st)][n, m] = np.max(np.abs(uex - vol_fcn_sum(st, uv, Eps, delta, K, Nz, n, m)))
                err[('w', st)][n, m] = np.max(np.abs(wex - vol_fcn_sum(st, wv, Eps, delta, K, Nz, n, m)))
    if verbose:
        print(f'Summation time: {time.time() - t0:.2f} s')

    if make_plots:
        os.makedirs(OUTDIR, exist_ok=True)
        E = lambda k, st: err[(k, st)]
        figs = [
            (1, 'Taylor', [E(k, 1) for k in ('xi_u', 'nu_u', 'zeta', 'xi_w', 'nu_w', 'psi')]),
            (3, 'Pade', [E(k, 2) for k in ('xi_u', 'nu_u', 'zeta', 'xi_w', 'nu_w', 'psi')]),
            (5, 'Pade', [E(k, 3) for k in ('xi_u', 'nu_u', 'zeta', 'xi_w', 'nu_w', 'psi')]),
            (2, 'Taylor', [E(k, 1) for k in ('G', 'U', 'ubar', 'J', 'W', 'wbar')]),
            (4, 'Pade', [E(k, 2) for k in ('G', 'U', 'ubar', 'J', 'W', 'wbar')]),
            (6, 'Pade Safe', [E(k, 3) for k in ('G', 'U', 'ubar', 'J', 'W', 'wbar')]),
            (11, 'Fields', [E('u', 1), E('u', 2), E('u', 3), E('w', 1), E('w', 2), E('w', 3)])]
        names = {1: (r'\xi_u', r'\nu_u', r'\zeta', r'\xi_w', r'\nu_w', r'\psi'),
                 11: ('u (Taylor)', 'u (Pade)', 'u (Pade safe)', 'w (Taylor)', 'w (Pade)', 'w (Pade safe)')}
        names[3] = names[5] = names[1]
        for num, lab, arrs in figs:
            kw = {'names': names[num]} if num in names else {}
            fig = plot_errors(num, lab + ' (3D)', N, M, *arrs, **kw)
            fig.savefig(os.path.join(OUTDIR, f'test_single_eps_delta_3D_fig{num}.png'), dpi=100)
    return err, dict(N=N, M=M, Eps=Eps, delta=delta, omega_bar=omega_bar)


KEYS = ('xi_u', 'nu_u', 'zeta', 'psi', 'G', 'J', 'U', 'W', 'ubar', 'wbar', 'u', 'w')


def _print_table(err, info):
    print(f"\nErrors at N={info['N']}, M={info['M']}, eps={info['Eps']:g}, delta={info['delta']:g} "
          f"(omega_bar = {info['omega_bar']:.4f}):")
    for k in KEYS:
        print(f"  {k:5s}  Taylor {err[(k, 1)][-1, -1]:.2e}   Pade {err[(k, 2)][-1, -1]:.2e}"
              f"   Pade-safe {err[(k, 3)][-1, -1]:.2e}")


def test_single_eps_delta_3D():
    """pytest entry point: RunNumber 1 at normal incidence, then an oblique crossed case."""
    matplotlib.use('Agg')
    err, info = run(make_plots=True, verbose=True)
    _print_table(err, info)
    for k in KEYS:
        for st in (1, 2, 3):
            assert err[(k, st)][-1, -1] < 1e-9, (k, st, err[(k, st)][-1, -1])
    err, info = run(RunNumber=2, make_plots=False, verbose=False, alpha=0.1, beta=0.2, profile='egg', rs=(1, 3))
    _print_table(err, info)
    for k in KEYS:      # eps = 0.1: Taylor is truncated at half order (MATLAB taylorsum_2_coeff), Pade converges
        assert err[(k, 2)][-1, -1] < 1e-4, (k, err[(k, 2)][-1, -1])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', type=int, default=RunNumber, help='RunNumber 1, 2 or 3 (as the MATLAB script)')
    ap.add_argument('--solver', default=DoTwoLayerTest, choices=list(h3.SOLVERS_3D))
    ap.add_argument('--alpha', type=float, default=0.0)
    ap.add_argument('--beta', type=float, default=0.0)
    ap.add_argument('--profile', default='cosxcosy')
    ap.add_argument('--no-show', action='store_true')
    a = ap.parse_args()
    if a.no_show:
        matplotlib.use('Agg')
    err, info = run(a.run, a.solver, alpha=a.alpha, beta=a.beta, profile=a.profile)
    _print_table(err, info)
    if SHOW and not a.no_show:
        plt.show()
