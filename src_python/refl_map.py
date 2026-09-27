"""refl_map.py -- Python port of refl_map.m (Reflectivity Map R and energy defect D).

Run directly (e.g. in PyCharm):   python refl_map.py            (default: SCENARIO below)
or from the command line:         python refl_map.py --scenario dielectric

Scenarios (paper = Kehoe & Nicholls, J. Sci. Comput. 100:9, 2024)
  'silver'      : refl_map.m as shipped in src.zip -> paper Fig. 10a
                  (n_w = 0.05+2.275i, f = cos(4x), N=M=15, Pade)
  'gold'        : paper Fig. 10b (n_w = 1.48+1.883i)
  'dielectric'  : plots/test_scenarios.m / GitHub README / paper Fig. 9
                  (n_w = 1.1, f = cos(x), N=M=16, Taylor)
  'dielectric_alpha' : paper Fig. 14 (alpha = 0.01; paper used N_eps=N_delta=1000)
"""
import argparse
import os
import time

# One BLAS thread per process: the bands run in parallel processes, and the matrices
# (33 x 33) are too small to benefit from threaded BLAS (avoids oversubscription).
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from hops import (cheb, setup_2d, csqrt, setup_zeta_psi_n_m, two_layer_solve_fast,
                  energy_defect, fourier_repr_lipschitz, fourier_repr_rough)
from hops.plotting import safe_log10
from hops.operators import two_layer_solve_operator, two_layer_solve_auto

SOLVERS = {'operator': two_layer_solve_operator, 'fast': two_layer_solve_fast,
           'auto': two_layer_solve_auto}

# ----------------------------------------------------------------------------
# User settings (mirror the top of refl_map.m)
# ----------------------------------------------------------------------------
SCENARIO = 'silver'
PlotLambda = 1
PlotRelative = 1
RunNumber = 100
Mode = 2                  # 1 = TE, 2 = TM
N_delta = 100
N_Eps = 100
SHOW = True               # plt.show() at the end (set False for batch runs)
SOLVER = 'operator'       # 'operator' (fast, new) | 'fast' (= two_layer_solve_fast.m) | 'auto'
WORKERS = os.cpu_count()  # frequency bands solved in parallel processes (1 = serial)
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')

SCENARIOS = {
    'silver':     dict(M=15, n_w=0.05 + 2.275j, profile='cos4x', Taylor=False, alpha_bar=0.0),
    'gold':       dict(M=15, n_w=1.48 + 1.883j, profile='cos4x', Taylor=False, alpha_bar=0.0),
    'dielectric': dict(M=16, n_w=1.1,           profile='cosx',  Taylor=True,  alpha_bar=0.0),
    'dielectric_alpha': dict(M=16, n_w=1.1,     profile='cosx',  Taylor=True,  alpha_bar=0.01),
}


def profile_fn(name, xx):
    if name == 'cos4x':
        return np.cos(4 * xx), -4 * np.sin(4 * xx)
    if name == 'cosx':
        return np.cos(xx), -np.sin(xx)
    if name == 'lipschitz':
        return fourier_repr_lipschitz(40, xx)
    if name == 'rough':
        return fourier_repr_rough(40, xx)
    raise ValueError(name)


def _band(q, c):
    """Everything for one frequency band q (module level so it can run in a worker process)."""
    Nx, Nz, N, M, d, c_0, alpha_bar = c['Nx'], c['Nz'], c['N'], c['M'], c['d'], c['c_0'], c['alpha_bar']
    n_u, n_w, f, f_x, Eps, N_Eps, N_delta, sigma = (c['n_u'], c['n_w'], c['f'], c['f_x'], c['Eps'],
                                                    c['N_Eps'], c['N_delta'], c['sigma'])
    identy = np.eye(Nz + 1)
    delta = np.array([0.0]) if N_delta == 1 else np.linspace(-sigma / (2 * q + 1), sigma / (2 * q + 1), N_delta)
    omega_bar = c_0 * (2 * np.pi / d) * (q + 0.5)
    omega = (1 + delta) * omega_bar
    lam = 2 * np.pi * c_0 / omega
    k_u_bar = n_u * omega_bar / c_0
    gamma_u_bar = csqrt(k_u_bar ** 2 - alpha_bar ** 2)
    xx, pp, alpha_bar_p, gamma_u_bar_p, _, _ = setup_2d(Nx, d, alpha_bar, gamma_u_bar)
    k_w_bar = n_w * omega_bar / c_0
    gamma_w_bar = csqrt(k_w_bar ** 2 - alpha_bar ** 2)
    xx, pp, alpha_bar_p, gamma_w_bar_p, _, _ = setup_2d(Nx, d, alpha_bar, gamma_w_bar)
    zeta_n_m, psi_n_m = setup_zeta_psi_n_m(xx, pp, alpha_bar, gamma_u_bar, f, f_x, Nx, N, M)
    tau2 = 1.0 if c['Mode'] == 1 else (n_u / n_w) ** 2

    t0 = time.time()
    solve = SOLVERS[c['solver']]
    U_n_m, W_n_m, ubar_n_m, wbar_n_m = solve(
        tau2, zeta_n_m, psi_n_m, gamma_u_bar_p, gamma_w_bar_p, N, Nx, f, f_x, pp, alpha_bar,
        gamma_u_bar, gamma_w_bar, c['Dz'], c['a'], c['b'], Nz, M, identy, alpha_bar_p)
    t1 = time.time()
    ub = np.transpose(ubar_n_m, (2, 1, 0))    # permute(.,[3 2 1]) -> (N+1, M+1, Nx)
    wb = np.transpose(wbar_n_m, (2, 1, 0))
    ee_flat, ru_flat, rl_flat = energy_defect(tau2, ub, wb, d, alpha_bar, gamma_u_bar, gamma_w_bar,
                                              Eps, delta, Nx, 0, 0, N_Eps, N_delta, 1)
    st = 1 if c['Taylor'] else 2
    ee, ru, rl = energy_defect(tau2, ub, wb, d, alpha_bar, gamma_u_bar, gamma_w_bar,
                               Eps, delta, Nx, N, M, N_Eps, N_delta, st)
    if c['verbose']:
        print(f'q = {q}: two_layer_solve ({c["solver"]}) {t1 - t0:.1f} s, '
              f'energy_defect {time.time() - t1:.1f} s', flush=True)
    RR = ru / ru_flat if PlotRelative else ru
    return dict(q=q, delta=delta, omega=omega, lam=lam, Eps=Eps, ee=ee, ru=ru, rl=rl,
                ru_flat=ru_flat, RR=RR, U_n_m=U_n_m, W_n_m=W_n_m,
                ubar_n_m=ubar_n_m, wbar_n_m=wbar_n_m)


def run(scenario=SCENARIO, qq=(1, 2, 3, 4, 5, 6), N_Eps=N_Eps, N_delta=N_delta,
        RunNumber=RunNumber, Mode=Mode, Taylor=None, verbose=True, solver=None, workers=None):
    solver = solver or SOLVER
    workers = WORKERS if workers is None else workers
    sc = SCENARIOS[scenario]
    if RunNumber == 1:
        M, Nx, Eps_Max, sigma = 4, 16, 1e-2, 1e-2
    elif RunNumber == 2:
        M, Nx, Eps_Max, sigma = 6, 24, 0.1, 0.1
    elif RunNumber == 3:
        M, Nx, Eps_Max, sigma = 8, 32, 0.1, 0.5
    else:                                  # 100: HOPS/AWE paper
        M, Nx, Eps_Max, sigma = sc['M'], 32, 0.2, 0.99
    N = M
    Nz = 32
    alpha_bar = sc['alpha_bar']
    d = 2 * np.pi
    c_0 = 1.0
    n_u = 1.0
    n_w = sc['n_w']
    Taylor = sc['Taylor'] if Taylor is None else Taylor
    a, b = 1.0, 1.0
    identy = np.eye(Nz + 1)
    Dz, _ = cheb(Nz)
    Eps = np.linspace(0, Eps_Max, N_Eps)
    xx = (d / Nx) * np.arange(Nx)
    f, f_x = profile_fn(sc['profile'], xx)

    common = dict(N_delta=N_delta, sigma=sigma, c_0=c_0, d=d, n_u=n_u, n_w=n_w, alpha_bar=alpha_bar,
                  Nx=Nx, Nz=Nz, N=N, M=M, f=f, f_x=f_x, Mode=Mode, Dz=Dz, a=a, b=b, Eps=Eps,
                  N_Eps=N_Eps, Taylor=Taylor, solver=solver, verbose=verbose)
    workers = min(workers or 1, len(qq))
    if workers > 1:
        # the frequency bands q are independent -> solve them in parallel processes
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(_band, qq, [common] * len(qq)))
    else:
        results = [_band(q, common) for q in qq]
    return results, dict(scenario=scenario, N=N, M=M, Nx=Nx, Nz=Nz, n_w=n_w, Taylor=Taylor)


def plot(results, info, outdir=OUTDIR, tag=None):
    """MATLAB-like rendering: each band q is its own contourf call with automatic
    levels, all sharing one colour scale ('colormap hot'), as in refl_map.m.
    For R, isolated Pade spikes above 1 (a handful of points out of 60,000 next to
    the plasmon resonance) are drawn in the top colour, so the colour scale ends at 1
    as in the paper."""
    from hops.plotting import matlab_contourf, shared_colorbar
    from matplotlib.colors import Normalize
    os.makedirs(outdir, exist_ok=True)
    tag = tag or info['scenario']
    xkey = 'lam' if PlotLambda else 'omega'
    xlabel = r'$\lambda$' if PlotLambda else r'$\omega$'
    figs = []
    for num, (title, getZ) in enumerate([('$D$', lambda r: safe_log10(r['ee'])),
                                         ('$R$', lambda r: np.real(r['RR']))], start=1):
        Zs = [getZ(r) for r in results]
        allz = np.concatenate([z[np.isfinite(z)].ravel() for z in Zs])
        vmin, vmax = allz.min(), allz.max()
        if num == 2 and PlotRelative and np.mean(allz > 1.0 + 1e-9) < 1e-3:
            vmax = 1.0
            Zs = [np.minimum(z, 1.0) for z in Zs]
        norm = Normalize(vmin=vmin, vmax=vmax)
        fig, ax = plt.subplots(num=num, figsize=(7, 5), clear=True)
        for r, Z in zip(results, Zs):
            matlab_contourf(ax, r[xkey], r['Eps'], Z, 'hot', norm)
        shared_colorbar(fig, ax, 'hot', norm)
        ax.set_xlabel(xlabel, fontsize=16)
        ax.set_ylabel(r'$\varepsilon$', fontsize=18)
        ax.set_title(title, fontsize=16)
        fig.tight_layout()
        fname = os.path.join(outdir, f'refl_map_{tag}_{"D" if num == 1 else "R"}.png')
        fig.savefig(fname, dpi=130)
        figs.append(fname)
    return figs


def load_results(scenario, outdir=OUTDIR):
    """Re-load a saved run (figures/refl_map_<scenario>.npz) for re-plotting."""
    d = np.load(os.path.join(outdir, f'refl_map_{scenario}.npz'))
    qs = sorted({int(k.split('_')[0][1:]) for k in d.files})
    return [dict(q=q, **{k: d[f'q{q}_{k}'] for k in ('lam', 'omega', 'Eps', 'ee', 'ru', 'RR')})
            for q in qs], dict(scenario=scenario)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--scenario', default=SCENARIO, choices=list(SCENARIOS))
    ap.add_argument('--q', type=int, nargs='*', default=[1, 2, 3, 4, 5, 6])
    ap.add_argument('--neps', type=int, default=N_Eps)
    ap.add_argument('--ndelta', type=int, default=N_delta)
    ap.add_argument('--solver', default=SOLVER, choices=['operator', 'fast', 'auto'])
    ap.add_argument('--workers', type=int, default=WORKERS)
    ap.add_argument('--no-show', action='store_true')
    args = ap.parse_args()
    if args.no_show:
        matplotlib.use('Agg')
    t0 = time.time()
    res, info = run(args.scenario, tuple(args.q), args.neps, args.ndelta,
                    solver=args.solver, workers=args.workers)
    for fn in plot(res, info):
        print('saved', fn)
    np.savez_compressed(os.path.join(OUTDIR, f'refl_map_{args.scenario}.npz'),
                        **{f'q{r["q"]}_{k}': r[k] for r in res for k in ('lam', 'omega', 'Eps', 'ee', 'ru', 'RR')})
    print(f'total {time.time() - t0:.1f} s')
    if SHOW and not args.no_show:
        plt.show()
