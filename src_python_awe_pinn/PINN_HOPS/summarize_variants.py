"""summarize_variants.py -- one table and one figure for all PINN variants vs HOPS/AWE.

Reads results/variants/*.json (gradient-trained PINN variants, DeepXDE), results/lsq/points.json
(least-squares interface PINN) and results/variants/timing.json; writes results/variants/summary.md
and summary.png.
"""
import glob
import json
import os

import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
V = os.path.join(HERE, 'results', 'variants')
LABEL = {'tanh': 'PINN tanh (baseline)', 'sin': 'PINN sin', 'laaf': 'PINN adaptive tanh (LAAF)',
         'ipinn_tanh_sin': 'I-PINN (tanh above, sin below)', 'tanh+lsq': 'PINN tanh + LSQ output layer',
         'lsgd': 'PINN, LSGD (Adam + LSQ output layer)', 'deepxde_tanh': 'DeepXDE PFNN tanh',
         'deepxde_LAAF-10_tanh': 'DeepXDE PFNN LAAF-10 tanh'}
ORDER = list(LABEL)


def load():
    rows = []
    for f in glob.glob(os.path.join(V, '*.json')):
        if os.path.basename(f) in ('timing.json',):
            continue
        r = json.load(open(f))
        if isinstance(r, dict) and 'variant' in r:
            rows.append(r)
    lsq = json.load(open(os.path.join(HERE, 'results', 'lsq', 'points.json')))
    for r in lsq:
        if r['case'] in ('dielectric', 'gold'):
            rows.append(dict(case=r['case'], variant='lsq', R_relerr=r['R_relerr'], D_absdiff=r['D_absdiff'],
                             U_relerr=r['U_relerr'], loss=r['loss'], t_train=r['t_lsq_s']))
    return rows


def main():
    rows = load()
    tim = json.load(open(os.path.join(V, 'timing.json'))) if os.path.exists(os.path.join(V, 'timing.json')) else {}
    order = ORDER + ['lsq']
    lab = dict(LABEL, lsq='least-squares interface PINN (this work)')
    lines = ['| variant | case | rel. err R | abs. err D | max rel. err U | PINN loss | wall time (s) |',
             '|---|---|---|---|---|---|---|']
    for case in ('dielectric', 'gold'):
        for v in order:
            for r in rows:
                if r['case'] == case and r['variant'] == v:
                    lines.append(f"| {lab[v]} | {case} | {r['R_relerr']:.1e} | {r['D_absdiff']:.1e} | {r['U_relerr']:.1e} | "
                                 f"{r['loss']:.1e} | {r['t_train']:.0f} |")
        if case in tim:
            lines.append(f"| HOPS/AWE (reference) | {case} | - | - | - | - | {tim[case]['hops_s']:.2f} |")
    open(os.path.join(V, 'summary.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.5))
    for ax, case in zip(axs, ('dielectric', 'gold')):
        names, vals, dvals = [], [], []
        for v in order:
            for r in rows:
                if r['case'] == case and r['variant'] == v:
                    names.append(lab[v])
                    vals.append(max(r['R_relerr'], 1e-16))
                    dvals.append(max(r['D_absdiff'], 1e-16))
        y = np.arange(len(names))
        ax.barh(y - 0.2, vals, 0.4, log=True, color='C0', label='relative error in R')
        ax.barh(y + 0.2, dvals, 0.4, log=True, color='C3', label='absolute error in D')
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=8)
        ax.set_xscale('log')
        ax.set_xlim(1e-16, 1)
        ax.axvline(1e-15, color='k', ls=':', lw=0.8)
        ax.set_title(f'{case}: error vs HOPS/AWE (same equations (6))')
        ax.grid(alpha=0.3, axis='x')
        ax.legend(fontsize=8, loc='lower right')
    fig.tight_layout()
    fig.savefig(os.path.join(V, 'summary.png'), dpi=110)


if __name__ == '__main__':
    main()
