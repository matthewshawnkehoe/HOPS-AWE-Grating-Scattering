"""compare_orders.py -- accuracy vs expansion order: HOPS/AWE with N = M in {6, 8, 10, 12, 16} summed as in
refl_map.py, vs the physics-informed (hybrid) sum of the SAME coefficient fields.

Reference: the cached pointwise truth of compare_maps.py (results/<scenario>/truth.npz).
Writes results/orders.json, results/orders.md, results/orders.png.
"""
import json
import os
import time
import warnings

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from hybrid import AWEBand, PISum
from compare_maps import SCENARIOS, OUT

CASES = ('dielectric_q1', 'high_index_q1', 'gold_q1_nx64')
ORDERS = (6, 8, 10, 12, 16)


def main():
    warnings.filterwarnings('ignore')
    fn = os.path.join(OUT, 'orders.json')
    rows = json.load(open(fn)) if os.path.exists(fn) else []
    done = {(r['case'], r['M']) for r in rows}
    for case in CASES:
        sc, q, over, n, _, _ = SCENARIOS[case]
        Rt = np.load(os.path.join(OUT, case, 'truth.npz'))['R']
        for M in ORDERS:
            if (case, M) in done:
                continue
            B = AWEBand(sc, q=q, N_Eps=n, N_delta=n, M=M, **{k: v for k, v in over.items() if k != 'M'})
            S = PISum(B, compress=1e-13)
            t0 = time.time()
            pi = S.map()
            rows.append(dict(case=case, M=M, awe_max=float(np.abs(B.R_awe - Rt).max()),
                             awe_median=float(np.median(np.abs(B.R_awe - Rt))),
                             taylor_full_max=float(np.abs(pi['R_taylor'] - Rt).max()),
                             hybrid_max=float(np.abs(pi['R'] - Rt).max()),
                             hybrid_median=float(np.median(np.abs(pi['R'] - Rt))),
                             t_band=B.t_hops, t_hybrid=time.time() - t0, unknowns=int(S.nfeat['u'] + S.nfeat['w'])))
            print(rows[-1], flush=True)
            json.dump(rows, open(fn, 'w'), indent=1)
    L = ['| scenario | N = M | AWE (refl_map) max err R | AWE Taylor full order | hybrid max err R | hybrid median | '
         'unknowns | time band / hybrid |', '|---|---|---|---|---|---|---|---|']
    for r in rows:
        L.append(f"| {r['case']} | {r['M']} | {r['awe_max']:.1e} | {r['taylor_full_max']:.1e} | {r['hybrid_max']:.1e} | "
                 f"{r['hybrid_median']:.1e} | {r['unknowns']} | {r['t_band']:.1f} s / {r['t_hybrid']:.0f} s |")
    open(os.path.join(OUT, 'orders.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    matplotlib.use('Agg')
    fig, axs = plt.subplots(1, len(CASES), figsize=(5.5 * len(CASES), 4.3))
    for ax, case in zip(axs, CASES):
        v = sorted([r for r in rows if r['case'] == case], key=lambda r: r['M'])
        Ms = [r['M'] for r in v]
        ax.semilogy(Ms, [r['awe_max'] for r in v], 'o-', label='HOPS/AWE (refl_map summation)')
        ax.semilogy(Ms, [r['taylor_full_max'] for r in v], 's--', label='AWE, full-order Taylor')
        ax.semilogy(Ms, [max(r['hybrid_max'], 1e-16) for r in v], 'o-', color='C3', label='hybrid (physics-informed sum)')
        ax.set_xlabel('expansion order N = M')
        ax.set_ylabel('max |R - R_true| over the map')
        ax.set_title(case)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'orders.png'), dpi=110)


if __name__ == '__main__':
    main()
