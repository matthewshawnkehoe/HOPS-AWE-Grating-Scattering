"""compare_refl_map.py -- parametric PINN vs HOPS/AWE reflectivity map R(eps, omega) and energy defect D.

HOPS/AWE (refl_map.py) makes ONE joint expansion in (eps, delta) per frequency band and sums it on the
whole (eps, omega) grid.  The PINN counterpart is a *parametric* PINN: the networks take (x, z, eps, omega)
as inputs and are trained once on the governing equations (6a)-(6h) sampled over the whole band
    eps in [0, eps_max],  omega in omega_bar (1 +- sigma/(2q+1))   (the band of refl_map.m),
after which R, T, D follow on any grid by evaluating the networks on z = a and z = -b.
Both use the same R, T, D formulas (paper Sect. 2.2 / energy_defect.m).

Outputs results/refl_map/<scenario>_q<q>.png/.npz:  R/R_flat and log10|D| of HOPS and PINN, |R_PINN - R_HOPS|,
and spectra R(lambda) at two eps values.

    python compare_refl_map.py                                  # dielectric (paper Fig. 9, band q = 1)
    python compare_refl_map.py --scenario gold --budget quick   # paper Fig. 10b band q = 1
"""
import argparse
import json
import os
import time

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from pinn_hops import Grating2D, GratingPINN
from pinn_hops.hops_reference import hops_refl_map

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'refl_map')

SCENARIOS = {                  # refl_map.py scenario name -> the same physics for the PINN
    'dielectric': dict(n_w=1.1, profile='cosx', K=6, eps_max=0.2),
    'silver': dict(n_w=0.05 + 2.275j, profile='cos4x', K=12, eps_max=0.2),
    'gold': dict(n_w=1.48 + 1.883j, profile='cos4x', K=12, eps_max=0.2),
}
BUDGETS = {'quick': dict(adam_iters=2000, lbfgs_iters=1000), 'standard': dict(adam_iters=4000, lbfgs_iters=4000),
           'long': dict(adam_iters=6000, lbfgs_iters=8000)}


def run(scenario='dielectric', q=1, budget='standard', n_grid=21, sigma=0.99, verbose=True, seed=0):
    s = SCENARIOS[scenario]
    omega_bar = q + 0.5
    dmax = sigma / (2 * q + 1)
    om_rng = (omega_bar * (1 - dmax), omega_bar * (1 + dmax))
    eps_rng = (0.0, s['eps_max'])
    # ---- HOPS/AWE (refl_map.py) on an n_grid x n_grid (eps, delta) grid
    t0 = time.time()
    res, info = hops_refl_map(scenario, q=(q,), N_Eps=n_grid, N_delta=n_grid)
    t_hops = time.time() - t0
    r = res[0]
    Eps, omega = r['Eps'], r['omega']
    # ---- parametric PINN over the same band
    G = Grating2D(n_w=s['n_w'], profile=s['profile'], eps=0.5 * s['eps_max'], omega=omega_bar)
    P = GratingPINN(G, K=s['K'], width=64, depth=4, eps_range=eps_rng, omega_range=om_rng, n_int=4000,
                    n_if=512, seed=seed)
    P.train(**BUDGETS[budget], verbose=verbose, log_every=1000)
    t1 = time.time()
    Rp = np.zeros((Eps.size, omega.size))
    Tp = np.zeros_like(Rp)
    Dp = np.zeros_like(Rp)
    for i, e in enumerate(Eps):
        for j, o in enumerate(omega):
            Rp[i, j], Tp[i, j], Dp[i, j] = P.energy(e, o)
    t_eval = time.time() - t1
    Rh, Th, Dh = np.real(r['ru']), np.real(r['rl']), np.real(r['ee'])
    Rflat = np.real(r['ru_flat'])
    summary = dict(scenario=scenario, q=q, budget=budget, grid=f'{n_grid}x{n_grid}',
                   hops_summation='Taylor' if info['Taylor'] else 'Pade',
                   t_hops_s=t_hops, t_pinn_train_s=P.train_time, t_pinn_eval_s=t_eval,
                   R_absdiff_max=float(np.abs(Rp - Rh).max()), R_absdiff_median=float(np.median(np.abs(Rp - Rh))),
                   R_reldiff_median=float(np.median(np.abs(Rp - Rh) / np.abs(Rh))),
                   D_hops_median_log10=float(np.log10(np.median(np.abs(Dh)) + 1e-300)),
                   D_hops_max_log10=float(np.log10(np.abs(Dh).max() + 1e-300)),
                   D_pinn_median_log10=float(np.log10(np.median(np.abs(Dp)) + 1e-300)),
                   D_pinn_max_log10=float(np.log10(np.abs(Dp).max() + 1e-300)),
                   lossless=bool(info['lossless']))
    if verbose:
        print(json.dumps(summary, indent=1))
    os.makedirs(OUT, exist_ok=True)
    tag = f'{scenario}_q{q}'
    np.savez_compressed(os.path.join(OUT, tag + '.npz'), Eps=Eps, omega=omega, R_hops=Rh, R_pinn=Rp, T_hops=Th,
                        T_pinn=Tp, D_hops=Dh, D_pinn=Dp, R_flat=Rflat)
    with open(os.path.join(OUT, tag + '.json'), 'w') as fh:
        json.dump(summary, fh, indent=1)
    P.save(os.path.join(OUT, tag + '_pinn.pt'))
    _figure(tag, Eps, omega, Rh, Rp, Dh, Dp, Rflat, summary, info)
    return summary


def _figure(tag, Eps, omega, Rh, Rp, Dh, Dp, Rflat, summary, info):
    lam = 2 * np.pi / omega
    fig, axs = plt.subplots(2, 4, figsize=(20, 8.5))
    RR = [Rh / Rflat, Rp / Rflat]
    lo, hi = min(z.min() for z in RR), max(z.max() for z in RR)
    for ax, Z, t in ((axs[0, 0], RR[0], 'HOPS/AWE (refl_map.py)'), (axs[0, 1], RR[1], 'parametric PINN')):
        m = ax.contourf(lam, Eps, Z, levels=np.linspace(lo, hi, 15), cmap='hot')
        ax.contour(lam, Eps, Z, levels=np.linspace(lo, hi, 15)[1:-1], colors='k', linewidths=0.4)
        fig.colorbar(m, ax=ax)
        ax.set_title(r'$R/R_{flat}$: ' + t)
    m = axs[0, 2].pcolormesh(lam, Eps, np.log10(np.abs(Rp - Rh) + 1e-17), cmap='viridis', shading='auto')
    fig.colorbar(m, ax=axs[0, 2])
    axs[0, 2].set_title(r'$\log_{10}|R_{PINN} - R_{HOPS}|$')
    for i, c in zip((len(Eps) // 2, len(Eps) - 1), ('C0', 'C3')):
        axs[0, 3].plot(lam, Rh[i], '-', color=c, lw=2, label=f'HOPS, eps={Eps[i]:.2g}')
        axs[0, 3].plot(lam, Rp[i], '--', color='k', lw=1, label=f'PINN, eps={Eps[i]:.2g}')
    axs[0, 3].plot(lam, Rflat[0], ':', color='0.5', label='flat')
    axs[0, 3].set_title('spectra R($\\lambda$)')
    axs[0, 3].legend(fontsize=7)
    Ds = [np.log10(np.abs(Dh) + 1e-17), np.log10(np.abs(Dp) + 1e-17)]
    lo, hi = min(z.min() for z in Ds), max(z.max() for z in Ds)
    lab = 'log10|D|' if info['lossless'] else 'log10 A (absorptance = D)'
    for ax, Z, t in ((axs[1, 0], Ds[0], 'HOPS/AWE'), (axs[1, 1], Ds[1], 'parametric PINN')):
        m = ax.contourf(lam, Eps, Z, levels=np.linspace(lo, hi, 15), cmap='hot')
        fig.colorbar(m, ax=ax)
        ax.set_title(f'{lab}: {t}')
    m = axs[1, 2].pcolormesh(lam, Eps, np.log10(np.abs(Dp - Dh) + 1e-17), cmap='viridis', shading='auto')
    fig.colorbar(m, ax=axs[1, 2])
    axs[1, 2].set_title(r'$\log_{10}|D_{PINN} - D_{HOPS}|$')
    axs[1, 3].axis('off')
    txt = '\n'.join(f'{k}: {v:.4g}' if isinstance(v, float) else f'{k}: {v}' for k, v in summary.items())
    axs[1, 3].text(0.0, 1.0, txt, va='top', family='monospace', fontsize=9)
    for ax in axs.ravel()[:7]:
        ax.set_xlabel(r'$\lambda$')
        ax.set_ylabel(r'$\varepsilon$')
    axs[0, 3].set_ylabel('R')
    fig.suptitle(f'{tag}: HOPS/AWE vs parametric PINN (same equations (6), same R/T/D formulas)', fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, tag + '.png'), dpi=100)
    plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scenario', default='dielectric', choices=list(SCENARIOS))
    ap.add_argument('--q', type=int, default=1)
    ap.add_argument('--budget', default='standard', choices=list(BUDGETS))
    ap.add_argument('--grid', type=int, default=21)
    a = ap.parse_args()
    matplotlib.use('Agg')
    run(a.scenario, a.q, a.budget, a.grid)
