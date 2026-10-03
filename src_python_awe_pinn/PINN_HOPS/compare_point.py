"""compare_point.py -- PINN vs HOPS/AWE at single (eps, omega) points.

For each case a PINN is trained on the governing equations (6a)-(6h) of Kehoe & Nicholls (2024) and
compared with the HOPS solution of the same problem (HOPS_Python, delta = 0, N = 16, Pade in eps):
  * R, T and the energy defect D = 1 - R - T  (the refl_map / energy_defect quantities)
  * interface data U = u(x, g), W = w(x, g) and the full fields on a uniform (x, z) grid
  * the PDE loss of both solutions (HOPS evaluated through its spectral interpolant, see
    hops_reference.HOPSFieldTorch) and the wall-clock times.
Outputs: results/point/<case>.png, results/point/summary.csv (and .json).

    python compare_point.py                     # all cases, 'standard' budget (~10 min per case)
    python compare_point.py --cases dielectric --budget quick
"""
import argparse
import csv
import json
import os
import time
from dataclasses import asdict

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from pinn_hops import Grating2D, GratingPINN
from pinn_hops.hops_reference import hops_point, HOPSFieldTorch

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'point')

CASES = {
    # paper Fig. 9 setting (vacuum / n = 1.1, f = cos x), centre of band q = 1 (omega = 3/2, lambda = 4.19)
    'dielectric': dict(G=dict(n_w=1.1, profile='cosx', eps=0.1, omega=1.5), net=dict(K=6)),
    'dielectric_eps0.2': dict(G=dict(n_w=1.1, profile='cosx', eps=0.2, omega=1.5), net=dict(K=6)),
    'dielectric_oblique': dict(G=dict(n_w=1.1, profile='cosx', eps=0.1, omega=1.5, alpha=0.1), net=dict(K=6)),
    # paper Fig. 10 materials (f = cos 4x), band q = 1 centre
    'silver': dict(G=dict(n_w=0.05 + 2.275j, profile='cos4x', eps=0.1, omega=1.5), net=dict(K=12)),
    'gold': dict(G=dict(n_w=1.48 + 1.883j, profile='cos4x', eps=0.1, omega=1.5), net=dict(K=12)),
}
BUDGETS = {'quick': dict(adam_iters=1500, lbfgs_iters=500), 'standard': dict(adam_iters=3000, lbfgs_iters=3000),
           'long': dict(adam_iters=5000, lbfgs_iters=6000)}


def run_case(name, budget='standard', verbose=True, seed=0):
    cfg = CASES[name]
    G = Grating2D(**cfg['G'])
    t0 = time.time()
    ref = hops_point(G, N=16, Nx=128, Nz=48, summation='pade')    # resolved: HOPS loss ~1e-15
    t_hops = time.time() - t0
    P = GratingPINN(G, width=48, depth=4, seed=seed, **cfg['net'])
    if verbose:
        print(f'== {name}: {asdict(G)}')
    P.train(**BUDGETS[budget], verbose=verbose, log_every=1000)
    R, T, D = P.energy()
    # interface data and fields
    x = ref['x']
    U = P.evaluate('u', x, G.g(x))
    W = P.evaluate('w', x, G.g(x))
    H = HOPSFieldTorch(ref, G)
    nx, nz = 96, 97
    X, Z = np.meshgrid(np.linspace(0, 2 * np.pi, nx), np.linspace(-G.b, G.a, nz), indexing='ij')
    up = Z > G.g(X)
    import torch
    with torch.no_grad():
        hu = H('u', torch.tensor(X[up]), torch.tensor(Z[up]), None, None)
        hw = H('w', torch.tensor(X[~up]), torch.tensor(Z[~up]), None, None)
    Fh = np.zeros(X.shape, complex)
    Fh[up] = hu[0].numpy() + 1j * hu[1].numpy()
    Fh[~up] = hw[0].numpy() + 1j * hw[1].numpy()
    Fp = np.zeros(X.shape, complex)
    Fp[up] = P.evaluate('u', X[up], Z[up])
    Fp[~up] = P.evaluate('w', X[~up], Z[~up])
    # PDE loss of the HOPS solution vs the PINN's own loss on the same collocation points
    P.n_int, P.n_if = 4000, 512
    P.sample()
    Lp, parts_p = P.loss()
    P.override = H
    Lh, parts_h = P.loss()
    P.override = None
    row = dict(case=name, n_w=str(G.n_w), profile=G.profile, eps=G.eps, omega=G.omega, alpha=G.alpha,
               R_hops=ref['R'], R_pinn=R, R_relerr=abs(R - ref['R']) / abs(ref['R']),
               T_hops=ref['T'], T_pinn=T, D_hops=ref['D'], D_pinn=D, D_abs_diff=abs(D - ref['D']),
               U_relerr=np.abs(U - ref['U']).max() / np.abs(ref['U']).max(),
               W_relerr=np.abs(W - ref['W']).max() / np.abs(ref['W']).max(),
               field_relerr_L2=np.linalg.norm(Fp - Fh) / np.linalg.norm(Fh),
               loss_pinn=float(Lp), loss_hops=float(Lh), t_pinn_s=P.train_time, t_hops_s=t_hops)
    if verbose:
        print('   ' + ', '.join(f'{k}={v:.4g}' if isinstance(v, float) else f'{k}={v}' for k, v in row.items()))
    _figure(name, G, X, Z, up, Fp, Fh, P, row)
    return row, P


def _figure(name, G, X, Z, up, Fp, Fh, P, row):
    os.makedirs(OUT, exist_ok=True)
    inc = np.where(up, np.exp(-1j * G.gamma_u * Z), 0)          # total field above = u + incident
    ph = np.exp(1j * G.alpha * X)
    Tp, Th = ph * (Fp + inc), ph * (Fh + inc)
    fig, axs = plt.subplots(1, 4, figsize=(19, 4.2))
    v = np.abs(Th).max()
    for ax, F, t in ((axs[0], Th, 'HOPS/AWE: Re $u_{tot}$, $w$'), (axs[1], Tp, 'PINN: Re $u_{tot}$, $w$')):
        m = ax.pcolormesh(X, Z, np.real(F), cmap='RdBu_r', vmin=-v, vmax=v, shading='gouraud')
        ax.plot(X[:, 0], G.g(X[:, 0]), 'k-', lw=1.2)
        ax.set_title(t)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$z$')
        fig.colorbar(m, ax=ax)
    E = np.log10(np.abs(Fp - Fh) + 1e-16)
    m = axs[2].pcolormesh(X, Z, E, cmap='viridis', shading='gouraud')
    axs[2].plot(X[:, 0], G.g(X[:, 0]), 'w-', lw=1)
    axs[2].set_title(r'$\log_{10}|$PINN $-$ HOPS$|$')
    fig.colorbar(m, ax=axs[2])
    h = P.history
    axs[3].semilogy([r[4] for r in h], [r[2] for r in h], 'k.-', label='total')
    for k, c in (('pde', 'C0'), ('if_dir', 'C1'), ('if_neu', 'C2'), ('tbc', 'C3')):
        axs[3].semilogy([r[4] for r in h], [r[3][k] for r in h], '.-', color=c, label=k, lw=0.8)
    axs[3].axhline(row['loss_hops'], color='m', ls='--', label='HOPS solution')
    axs[3].set_xlabel('training time (s)')
    axs[3].set_title('PINN loss')
    axs[3].legend(fontsize=7)
    fig.suptitle(f"{name}: $n^w$={G.n_w}, f={G.profile}, $\\varepsilon$={G.eps}, $\\omega$={G.omega}, "
                 f"$\\alpha$={G.alpha}  |  R: HOPS {row['R_hops']:.6g}, PINN {row['R_pinn']:.6g}  |  "
                 f"D: HOPS {row['D_hops']:.1e}, PINN {row['D_pinn']:.1e}  |  "
                 f"time: HOPS {row['t_hops_s']:.2f} s, PINN {row['t_pinn_s']:.0f} s", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f'{name}.png'), dpi=110)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cases', nargs='*', default=list(CASES), choices=list(CASES))
    ap.add_argument('--budget', default='standard', choices=list(BUDGETS))
    a = ap.parse_args(argv)
    matplotlib.use('Agg')
    rows = []
    for c in a.cases:
        rows.append(run_case(c, a.budget)[0])
        os.makedirs(OUT, exist_ok=True)
        with open(os.path.join(OUT, 'summary.csv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        with open(os.path.join(OUT, 'summary.json'), 'w') as fh:
            json.dump(rows, fh, indent=1, default=str)
    print(f"saved {os.path.join(OUT, 'summary.csv')}")


if __name__ == '__main__':
    main()
