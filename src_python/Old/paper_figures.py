"""paper_figures.py -- regenerate the MMS error figures of the HOPS/AWE paper
(J. Sci. Comput. 100:9, 2024) with the Python port.

  Fig. 2 : N=M=4 for eps_max = 1e-2, 1e-4, 1e-6, 1e-8   (Taylor, relative error)
  Fig. 3 : (N=M, eps_max) = (4,1e-2), (8,1e-4), (12,1e-6), (16,1e-8)
  Fig. 4/5 : as Fig. 2/3 with n_w = 10.1  (use --nw 10.1)
Physical parameters (33): d=2pi, alpha=0, n_u=1, n_w=1.1, r=4, A_r=5, B_r=3, TM.
Numerical parameters (34): Nx=Nz=32, a=1, b=1 (lower boundary at z=-1).

Reflectivity maps (Figs. 9, 10) come from refl_map.py --scenario dielectric|silver|gold.
"""
import argparse
import os

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import mms_error as me
from hops.plotting import parula, matlab_contourf, shared_colorbar, safe_log10

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')
CASES = {2: [(4, 1e-2), (4, 1e-4), (4, 1e-6), (4, 1e-8)],
         3: [(4, 1e-2), (8, 1e-4), (12, 1e-6), (16, 1e-8)]}


def make(fig_no, n_w=1.1, show=False):
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
    for ax, (nm, eps), lab in zip(axs.ravel(), CASES[fig_no], 'abcd'):
        out, c = me.run(dict(N=nm, M=nm, Eps_Max=eps, a=1.0, b=1.0, n_w=n_w, relative=True,
                             sum_types=(1,)), verbose=False)
        o = out[0]
        Z = safe_log10(o['err'][('U', 1)])
        fin = Z[np.isfinite(Z)]
        norm = Normalize(fin.min(), fin.max())
        matlab_contourf(ax, o['omega'], o['Eps'], Z, parula, norm)
        shared_colorbar(fig, ax, parula, norm)
        ax.set_title('Relative Error')
        ax.set_xlabel(r'$\omega=\omega_1(1+\delta)$' + f'\n({lab}) $N=M={nm},\\ \\varepsilon={eps:g}$')
        ax.set_ylabel(r'$\varepsilon$')
        print(f'Fig {fig_no}{lab}: N=M={nm}, eps={eps:g}: log10 rel. error in [{fin.min():.2f}, {fin.max():.2f}]')
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    fn = os.path.join(OUT, f'paper_fig{fig_no}_nw{n_w:g}.png')
    fig.savefig(fn, dpi=120)
    print('saved', fn)
    return fn


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--fig', type=int, nargs='*', default=[2, 3])
    ap.add_argument('--nw', type=float, default=1.1)
    ap.add_argument('--no-show', action='store_true')
    a = ap.parse_args()
    if a.no_show:
        matplotlib.use('Agg')
    for f in a.fig:
        make(f, a.nw)
    if not a.no_show:
        plt.show()
