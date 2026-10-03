"""summarize.py -- tables and figures across all compare_maps.py scenarios.

  results/summary_all.md   per scenario: max / median |R - R_true| for AWE (refl_map summation), AWE with
                           full-order Taylor, hybrid; energy defect; times
  results/adaptive.md      the adaptive hybrid for a sweep of indicator tolerances (computed offline from
                           the saved maps): fraction of points that need a least-squares solve, and max error
  results/summary.png      bar chart of the max error in R per scenario
"""
import json
import os

import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from compare_maps import SCENARIOS, OUT

TOLS = (1e-18, 1e-14, 1e-10, 1e-6)


def main():
    rows, adapt = [], []
    for name in SCENARIOS:
        f = os.path.join(OUT, name, 'maps.npz')
        if not os.path.exists(f):
            continue
        d = np.load(f)
        s = json.load(open(os.path.join(OUT, name, 'summary.json')))
        Rt, Dt = d['R_true'], d['D_true']
        e = lambda R: np.abs(R - Rt)
        r = dict(name=name, n_w=s['n_w'], mode=s['mode'], summation=s['summation'], NM=f"{s['N']}/{s['M']}", Nx=s['Nx'],
                 awe=e(d['AWE_R']), full=e(d['R_taylor']), pi=e(d['PI_R']),
                 d_awe=np.abs(d['AWE_D'] - Dt).max(), d_pi=np.abs(d['PI_D'] - Dt).max(),
                 t_awe=s['AWE']['time_s'], t_pi=s['PI']['time_s'], t_truth=s['t_truth_s'], spear=s['indicator_corr'],
                 rf=e(d['PI+RF_R']) if 'PI+RF_R' in d.files else None)
        rows.append(r)
        n = Rt.size
        per = s['PI']['time_s'] / n
        for tol in TOLS:
            solve = d['indicator'] > tol
            Ra = np.where(solve, d['PI_R'], d['R_taylor'])
            adapt.append(dict(name=name, tol=tol, frac=float(solve.mean()), max_err=float(e(Ra).max()),
                              time_est=float(solve.sum() * per + n * 0.006)))
    L = ['| scenario | AWE summation | max err R: AWE | AWE full-order Taylor | **hybrid** | median err R: AWE / hybrid | '
         'max err D: AWE / hybrid | time: AWE / hybrid / pointwise HOPS | indicator vs error (Spearman) |',
         '|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        extra = f" (+RF: {r['rf'].max():.1e})" if r['rf'] is not None else ''
        L.append(f"| {r['name']} (n_w = {r['n_w'].strip('()').replace('+0j', '')}, {r['mode']}, Nx = {r['Nx']}) | "
                 f"{r['summation']} {r['NM']} | {r['awe'].max():.1e} | {r['full'].max():.1e} | **{r['pi'].max():.1e}**{extra} | "
                 f"{np.median(r['awe']):.1e} / {np.median(r['pi']):.1e} | {r['d_awe']:.1e} / {r['d_pi']:.1e} | "
                 f"{r['t_awe']:.1f} / {r['t_pi']:.0f} / {r['t_truth']:.0f} s | {r['spear']:.2f} |")
    open(os.path.join(OUT, 'summary_all.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    A = ['| scenario | ' + ' | '.join(f'tol {t:.0e}: solved / max err' for t in TOLS) + ' |', '|---' * (len(TOLS) + 1) + '|']
    for name in [r['name'] for r in rows]:
        v = [a for a in adapt if a['name'] == name]
        A.append(f'| {name} | ' + ' | '.join(f"{a['frac']:.0%} / {a['max_err']:.1e}" for a in v) + ' |')
    open(os.path.join(OUT, 'adaptive.md'), 'w').write('\n'.join(A) + '\n')
    print('\n'.join(A))
    fig, ax = plt.subplots(figsize=(13, 5))
    y = np.arange(len(rows))
    ax.bar(y - 0.27, [r['awe'].max() for r in rows], 0.27, label='HOPS/AWE (refl_map.py)', color='C0', log=True)
    ax.bar(y, [r['full'].max() for r in rows], 0.27, label='AWE, full-order Taylor', color='C9', log=True)
    ax.bar(y + 0.27, [max(r['pi'].max(), 1e-16) for r in rows], 0.27, label='HOPS/AWE + PINN (hybrid)', color='C3', log=True)
    ax.set_xticks(y)
    ax.set_xticklabels([r['name'] for r in rows], rotation=30, ha='right', fontsize=8)
    ax.set_ylabel('max |R - R_true| over the (eps, lambda) map')
    ax.set_ylim(1e-16, 1)
    ax.grid(alpha=0.3, axis='y')
    ax.legend(fontsize=8)
    ax.set_title('Reflectivity-map error: HOPS/AWE vs the physics-informed summation of the same HOPS/AWE fields')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'summary.png'), dpi=110)


if __name__ == '__main__':
    main()
