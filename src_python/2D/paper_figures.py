"""paper_figures.py -- reproduce Figs. 2-8 of Kehoe & Nicholls, "A stable HOPS/AWE method
for the numerical solution of grating scattering problems", J. Sci. Comput. 100:9 (2024).

All figures plot the relative error  |U_exact - U_approx|_inf / |U_exact|_inf  of the
upper-layer Dirichlet data for the manufactured solutions of Sec. 6.5
(r = 4, A_r = 5, B_r = 3, d = 2 pi, alpha = 0, n_u = 1, TM, sigma = 0.99, omega_1 = 3/2).

  Fig. 2 : n_w = 1.1,  N=M=4, eps_max = 1e-2, 1e-4, 1e-6, 1e-8          Taylor  (33),(34)
  Fig. 3 : n_w = 1.1,  (N=M, eps_max) = (4,1e-2) (8,1e-4) (12,1e-6) (16,1e-8)  Taylor
  Fig. 4 : as Fig. 2 with n_w = 10.1                                     (35),(34)
  Fig. 5 : as Fig. 3 with n_w = 10.1
           (34): Nx = Nz = 32, a = 1, b = 1 (lower boundary z = -1), f = cos(4x)/4
  Fig. 6 : smooth   f_s = cos(4x)/4,  Nx = 256,  Nz = 128, N=M=20, eps_max = 2, a = b = 4, Pade (39)
  Fig. 7 : rough    f_{r,P}, P = 120, Nx = 1024, Nz = 128, N=M=20, eps_max = 2, a = b = 4, Pade (40)
  Fig. 8 : Lipschitz f_{L,P}, P = 120, Nx = 1024, Nz = 128, ...                          (40)
           Figs. 6-8: left panel n_w = 1.1, right panel n_w = 10.1

Resolution for Figs. 6-8 (--res):
  reduced (default)  Fig. 6: Nx=128, Nz=64 (n_w=1.1) / Nz=128 (n_w=10.1); Fig. 7a: Nx=512,
                     Nz=64.  Figs. 7b and 8 keep the paper resolution (see REDUCED below).
  paper              Nx=256 (Fig. 6) / 1024 (Figs. 7, 8), Nz=128, as in (39)/(40).
Run times (one core, this port):  Figs. 2-5: ~10-30 s each.
  reduced: Fig. 6 ~5 / 15 min per panel, Fig. 7a ~20 min; Figs. 7b, 8 as paper.
  paper:   Fig. 6 ~30 min per panel, Figs. 7, 8: ~1.5-2 h per panel (Nx = 1024, Nz = 128).
  The panels are independent: --workers 6 runs them in parallel processes
  (each Nx = 1024 panel needs ~1.5 GB of RAM).  --quick gives a fast low-resolution
  preview of Figs. 6-8 (Nx = 256, Nz = 64, N = M = 12, 50 x 50 grid) -- NOT paper resolution.
Every panel is cached in figures/paper_cache/*.npz; --replot redraws from the cache.

Usage:   python paper_figures.py                   (Figs. 2-5)
         python paper_figures.py --fig 6 7 8 --workers 6
         python paper_figures.py --fig 2 3 4 5 6 7 8 --quick
"""
import argparse
import os
import time

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import mms_error as me
from hops.plotting import parula, matlab_contourf, shared_colorbar, safe_log10

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')
CACHE = os.path.join(OUT, 'paper_cache')

SMALL = dict(Nx=32, Nz=32, a=1.0, b=1.0, profile='cos4x_over4', sum_types=(1,))       # (34), Taylor
LARGE = dict(Nz=128, N=20, M=20, Eps_Max=2.0, a=4.0, b=4.0, sum_types=(2,))           # (39)/(40), Pade
PROFILE_NAME = {'cos4x_over4': r'$f_s$', 'rough': r'$f_{r,P}$', 'lipschitz': r'$f_{L,P}$'}

# figure -> list of panel configurations
FIGS = {
    2: [dict(SMALL, N=4, M=4, Eps_Max=e, n_w=1.1) for e in (1e-2, 1e-4, 1e-6, 1e-8)],
    3: [dict(SMALL, N=nm, M=nm, Eps_Max=e, n_w=1.1)
        for nm, e in ((4, 1e-2), (8, 1e-4), (12, 1e-6), (16, 1e-8))],
    4: [dict(SMALL, N=4, M=4, Eps_Max=e, n_w=10.1) for e in (1e-2, 1e-4, 1e-6, 1e-8)],
    5: [dict(SMALL, N=nm, M=nm, Eps_Max=e, n_w=10.1)
        for nm, e in ((4, 1e-2), (8, 1e-4), (12, 1e-6), (16, 1e-8))],
    6: [dict(LARGE, Nx=256, profile='cos4x_over4', n_w=nw) for nw in (1.1, 10.1)],
    7: [dict(LARGE, Nx=1024, profile='rough', n_w=nw) for nw in (1.1, 10.1)],
    8: [dict(LARGE, Nx=1024, profile='lipschitz', n_w=nw) for nw in (1.1, 10.1)],
}
QUICK = dict(Nx=256, Nz=64, N=12, M=12, N_Eps=50, N_delta=50)

# Reduced resolution (default) for Figs. 6-8, chosen from a convergence study against the
# paper resolution (see README "Figs. 6-8: how much resolution is needed"):
#   Fig. 6 (f_s, smooth): Nx = 128 is enough (Nx = 64 is not).  For n_w = 1.1, Nz = 64 is
#     converged and even ~0.7 digit MORE accurate than Nz = 128 (less round-off); n_w = 10.1
#     (k_w b ~ 60 in the lower layer) needs Nz = 128.   -> 2x-6x faster.
#   Fig. 7a (f_r, C^4, n_w = 1.1): Nx = 512, Nz = 64 reproduces Nx = 1024 (86 % of the cells
#     within 0.5 digit, same 95th percentile) and the paper.   -> ~8x faster.
#   Fig. 7b (n_w = 10.1): Nx = 512, Nz = 128 is NOT converged -> paper resolution kept.
#   Fig. 8 (f_L, Lipschitz): needs Nx = 1024 AND Nz = 128 (the slowly decaying Fourier modes
#     of the kinked profile create z-boundary layers); 512 x 128 is 4 digits worse
#     -> paper resolution kept.
REDUCED = {
    (6, 1.1): dict(Nx=128, Nz=64), (6, 10.1): dict(Nx=128, Nz=128),
    (7, 1.1): dict(Nx=512, Nz=64), (7, 10.1): dict(Nx=1024, Nz=128),
    (8, 1.1): dict(Nx=1024, Nz=128), (8, 10.1): dict(Nx=1024, Nz=128),
}


def panel_config(fig_no, i, quick=False, res='reduced'):
    cfg = dict(FIGS[fig_no][i], relative=True, quantities=('U',), solver='auto')
    if fig_no >= 6:
        if quick:
            cfg.update(QUICK)
        elif res == 'reduced':
            cfg.update(REDUCED[(fig_no, cfg['n_w'])])
    return cfg


def _key(cfg):
    return (f"{cfg['profile']}_Nx{cfg['Nx']}_Nz{cfg['Nz']}_NM{cfg['N']}_eps{cfg['Eps_Max']:g}"
            f"_nw{cfg['n_w']:g}_a{cfg['a']:g}_st{cfg['sum_types'][0]}"
            f"_g{cfg.get('N_Eps', 100)}x{cfg.get('N_delta', 100)}")


def compute_panel(cfg, verbose=True):
    """Run (or load from cache) one panel; returns (omega, Eps, relative error)."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, _key(cfg) + '.npz')
    if os.path.exists(path):
        d = np.load(path)
        return d['omega'], d['Eps'], d['err']
    t0 = time.time()
    if verbose:
        print(f'computing {_key(cfg)} ...', flush=True)
    out, c = me.run(cfg, verbose=verbose)
    o = out[0]
    st = cfg['sum_types'][0]
    err = o['err'][('U', st)]
    np.savez_compressed(path, omega=o['omega'], Eps=o['Eps'], err=err)
    if verbose:
        print(f'  done {_key(cfg)} in {time.time() - t0:.0f} s', flush=True)
    return o['omega'], o['Eps'], err


def _compute_star(args):
    return compute_panel(*args)


def make(fig_no, quick=False, workers=1, replot=False, res='reduced'):
    cfgs = [panel_config(fig_no, i, quick, res) for i in range(len(FIGS[fig_no]))]
    if replot:
        missing = [c for c in cfgs if not os.path.exists(os.path.join(CACHE, _key(c) + '.npz'))]
        if missing:
            raise FileNotFoundError(f'no cached data for Fig. {fig_no}; run without --replot first')
    if workers > 1 and len(cfgs) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=min(workers, len(cfgs))) as ex:
            panels = list(ex.map(_compute_star, [(c, True) for c in cfgs]))
    else:
        panels = [compute_panel(c) for c in cfgs]
    return plot(fig_no, cfgs, panels, quick, res)


def plot(fig_no, cfgs, panels, quick=False, res='reduced'):
    two = len(cfgs) == 2
    fig, axs = plt.subplots(1 if two else 2, 2, figsize=(12, 4.8) if two else (11, 8.5), squeeze=False)
    sum_lab = 'Taylor' if cfgs[0]['sum_types'][0] == 1 else 'Pade'
    for ax, cfg, (omega, Eps, err), lab in zip(axs.ravel(), cfgs, panels, 'abcd'):
        Z = safe_log10(err)
        fin = Z[np.isfinite(Z)]
        norm = Normalize(fin.min(), fin.max())
        matlab_contourf(ax, omega, Eps, Z, parula, norm)
        shared_colorbar(fig, ax, parula, norm)
        ax.set_title('Relative Error')
        ax.set_ylabel(r'$\varepsilon$')
        if fig_no <= 5:
            sub = f'({lab}) $N=M={cfg["N"]},\\ \\varepsilon={cfg["Eps_Max"]:g}$'
        else:
            sub = (f'({lab}) $N=M={cfg["N"]},\\ \\varepsilon={cfg["Eps_Max"]:g}$, '
                   f'$n^w={cfg["n_w"]:g}$, {PROFILE_NAME[cfg["profile"]]}, '
                   f'$N_x={cfg["Nx"]},\\ N_z={cfg["Nz"]}$')
        ax.set_xlabel(r'$\omega=\omega_1(1+\delta)$' + '\n' + sub)
        print(f'Fig {fig_no}{lab}: N=M={cfg["N"]}, eps={cfg["Eps_Max"]:g}, n_w={cfg["n_w"]:g}, '
              f'Nx={cfg["Nx"]}: log10 rel. error in [{fin.min():.2f}, {fin.max():.2f}]')
    tag = ''
    if fig_no >= 6:
        tag = (' (quick preview)' if quick else
               ' (paper resolution)' if res == 'paper' else ' (reduced, converged resolution)')
    fig.suptitle(f'Fig. {fig_no}: relative error in $U$, {sum_lab} summation{tag}', fontsize=12)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    suffix = ('_quick' if quick else '_paperres' if res == 'paper' else '') if fig_no >= 6 else ''
    fn = os.path.join(OUT, f'paper_fig{fig_no}{suffix}.png')
    fig.savefig(fn, dpi=120)
    print('saved', fn)
    return fn


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fig', type=int, nargs='*', default=[2, 3, 4, 5], choices=sorted(FIGS))
    ap.add_argument('--quick', action='store_true', help='low-resolution preview of Figs. 6-8')
    ap.add_argument('--res', choices=['reduced', 'paper'], default='reduced',
                    help="Figs. 6-8: 'reduced' (default, converged, 2-6x faster) or 'paper' (Nx=256/1024, Nz=128)")
    ap.add_argument('--workers', type=int, default=1, help='panels computed in parallel processes')
    ap.add_argument('--replot', action='store_true', help='only redraw from figures/paper_cache')
    ap.add_argument('--no-show', action='store_true')
    a = ap.parse_args()
    if a.no_show:
        matplotlib.use('Agg')
    for f in a.fig:
        make(f, a.quick, a.workers, a.replot, a.res)
    if not a.no_show:
        plt.show()
