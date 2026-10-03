"""paper_figures_3D.py -- the 3D (doubly periodic) analogues of Figs. 2-8 of Kehoe & Nicholls,
"A stable HOPS/AWE method for the numerical solution of grating scattering problems",
J. Sci. Comput. 100:9 (2024)  (2D version: paper_figures.py).

Every panel is the relative error  |U_exact - U_approx|_inf / |U_exact|_inf  of the upper-layer
Dirichlet data for the manufactured solution  u_rs = A exp(i p_r x + i q_s y + i gamma_rs z)
(rs = (4, 0), A = 5, B = 3, alpha = beta = 0, n_u = 1, TM, sigma = 0.99), as a function of
(omega, eps) over the first 3D frequency window omega in [1, sqrt 2] (omega_bar = 1.2071).

  Fig. 2 : n_w = 1.1,  N=M=4, eps_max = 1e-2, 1e-4, 1e-6, 1e-8                         Taylor
  Fig. 3 : n_w = 1.1,  (N=M, eps_max) = (4,1e-2) (8,1e-4) (12,1e-6) (16,1e-8)           Taylor
  Fig. 4 : as Fig. 2 with n_w = 10.1
  Fig. 5 : as Fig. 3 with n_w = 10.1
           Nx = Ny = Nz = 32, a = b = 1, f = f_s = cos(4x)cos(4y)/4        (2D: Nx = Nz = 32, cos(4x)/4)
  Fig. 6 : smooth    f_s,                       eps_max = 1 (2D: 2), a = b = 4, Pade
  Fig. 7 : rough     f_{r,P}(x) + f_{r,P}(y) (/2), P = 20                   (2D: P = 120, Nx = 1024)
  Fig. 8 : Lipschitz f_{L,P}(x) + f_{L,P}(y) (/2), P = 20
           Figs. 6-8: left panel n_w = 1.1, right panel n_w = 10.1.
Resolution for Figs. 6-8 (3D: Nx*Ny unknowns per Chebyshev level, so the 2D Nx = 256/1024 is out of
reach):  default  Nx = Ny = 64, Nz = 48 (n_w = 1.1) / 64 (n_w = 10.1), N = M = 16, 40 x 40 grid,
eps_max = 1 (~2 min and ~1.7 GB RAM per panel).  Convergence study (f_s, n_w = 1.1, rel. error at
eps = 0.5 / 1 / 2):  Nx = 32: 6e-4 / 5e-3 / 2e-2;  Nx = 48: 2e-6 / 5e-5 / 8e-4;  Nx = 64: 3e-9 / 3e-7 / 2e-5
-- the large-eps error is spatial resolution, exactly as in 2D (where Nx = 64 was not enough and
Nx = 128 was).  --eps-max 2 --nx 96 approaches the paper setting (needs ~4 GB per panel).
--quick  Nx = Ny = 32, Nz = 32, N = M = 12, 30 x 30 grid.  Figs. 2-5: ~10 s per panel.
Every panel is cached in figures_3d/paper/paper_cache/*.npz; --replot redraws from the cache.

Usage:   python paper_figures_3D.py                      (Figs. 2-5)
         python paper_figures_3D.py --fig 6 7 8 --workers 2
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

import mms_error_3D as me
from hops.plotting import parula, matlab_contourf, shared_colorbar, safe_log10

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures_3d', 'paper')
CACHE = os.path.join(OUT, 'paper_cache')

SMALL = dict(Nx=32, Ny=32, Nz=32, a=1.0, b=1.0, profile='fs', sum_types=(1,), N_Eps=50, N_delta=50)
LARGE = dict(N=16, M=16, Eps_Max=1.0, a=4.0, b=4.0, sum_types=(2,), N_Eps=40, N_delta=40)
PROFILE_NAME = {'fs': r'$f_s$', 'rough:20': r'$f_{r,P}$', 'lipschitz:20': r'$f_{L,P}$'}

FIGS = {
    2: [dict(SMALL, N=4, M=4, Eps_Max=e, n_w=1.1) for e in (1e-2, 1e-4, 1e-6, 1e-8)],
    3: [dict(SMALL, N=nm, M=nm, Eps_Max=e, n_w=1.1) for nm, e in ((4, 1e-2), (8, 1e-4), (12, 1e-6), (16, 1e-8))],
    4: [dict(SMALL, N=4, M=4, Eps_Max=e, n_w=10.1) for e in (1e-2, 1e-4, 1e-6, 1e-8)],
    5: [dict(SMALL, N=nm, M=nm, Eps_Max=e, n_w=10.1) for nm, e in ((4, 1e-2), (8, 1e-4), (12, 1e-6), (16, 1e-8))],
    6: [dict(LARGE, profile='fs', n_w=nw) for nw in (1.1, 10.1)],
    7: [dict(LARGE, profile='rough:20', n_w=nw) for nw in (1.1, 10.1)],
    8: [dict(LARGE, profile='lipschitz:20', n_w=nw) for nw in (1.1, 10.1)],
}
DEFAULT_LARGE = {1.1: dict(Nx=64, Ny=64, Nz=48), 10.1: dict(Nx=64, Ny=64, Nz=64)}
QUICK = dict(Nx=32, Ny=32, Nz=32, N=12, M=12, N_Eps=30, N_delta=30)


def panel_config(fig_no, i, quick=False, eps_max=None, nx=None):
    cfg = dict(FIGS[fig_no][i], relative=True, quantities=('U',), solver='coupled', windows=(1,))
    if fig_no >= 6:
        cfg.update(QUICK if quick else DEFAULT_LARGE[cfg['n_w']])
        if eps_max:
            cfg['Eps_Max'] = eps_max
        if nx:
            cfg.update(Nx=nx, Ny=nx)
    return cfg


def _key(cfg):
    return (f"{cfg['profile'].replace(':', '')}_Nx{cfg['Nx']}_Nz{cfg['Nz']}_NM{cfg['N']}_eps{cfg['Eps_Max']:g}"
            f"_nw{cfg['n_w']:g}_a{cfg['a']:g}_st{cfg['sum_types'][0]}_g{cfg['N_Eps']}x{cfg['N_delta']}")


def compute_panel(cfg, verbose=True):
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
    err = o['err'][('U', cfg['sum_types'][0])]
    np.savez_compressed(path, omega=o['omega'], Eps=o['Eps'], err=err)
    if verbose:
        print(f'  done {_key(cfg)} in {time.time() - t0:.0f} s', flush=True)
    return o['omega'], o['Eps'], err


def _compute_star(args):
    return compute_panel(*args)


def make(fig_no, quick=False, workers=1, replot=False, eps_max=None, nx=None):
    cfgs = [panel_config(fig_no, i, quick, eps_max, nx) for i in range(len(FIGS[fig_no]))]
    if replot and any(not os.path.exists(os.path.join(CACHE, _key(c) + '.npz')) for c in cfgs):
        raise FileNotFoundError(f'no cached data for Fig. {fig_no}; run without --replot first')
    if workers > 1 and len(cfgs) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=min(workers, len(cfgs))) as ex:
            panels = list(ex.map(_compute_star, [(c, True) for c in cfgs]))
    else:
        panels = [compute_panel(c) for c in cfgs]
    return plot(fig_no, cfgs, panels, quick)


def plot(fig_no, cfgs, panels, quick=False):
    two = len(cfgs) == 2
    fig, axs = plt.subplots(1 if two else 2, 2, figsize=(12, 4.8) if two else (11, 8.5), squeeze=False)
    sum_lab = 'Taylor' if cfgs[0]['sum_types'][0] == 1 else 'Pade'
    for ax, cfg, (omega, Eps, err), lab in zip(axs.ravel(), cfgs, panels, 'abcd'):
        Z = safe_log10(err)
        fin = Z[np.isfinite(Z)]
        norm = Normalize(fin.min(), fin.max())
        matlab_contourf(ax, omega, Eps, Z, parula, norm)
        shared_colorbar(fig, ax, parula, norm)
        ax.set_title('Relative Error (3D)')
        ax.set_ylabel(r'$\varepsilon$')
        if fig_no <= 5:
            sub = f'({lab}) $N=M={cfg["N"]},\\ \\varepsilon={cfg["Eps_Max"]:g}$'
        else:
            sub = (f'({lab}) $N=M={cfg["N"]},\\ \\varepsilon={cfg["Eps_Max"]:g}$, $n^w={cfg["n_w"]:g}$, '
                   f'{PROFILE_NAME[cfg["profile"]]}, $N_x=N_y={cfg["Nx"]},\\ N_z={cfg["Nz"]}$')
        ax.set_xlabel(r'$\omega=\bar\omega(1+\delta)$' + '\n' + sub)
        print(f'3D Fig {fig_no}{lab}: N=M={cfg["N"]}, eps={cfg["Eps_Max"]:g}, n_w={cfg["n_w"]:g}, '
              f'Nx=Ny={cfg["Nx"]}: log10 rel. error in [{fin.min():.2f}, {fin.max():.2f}]')
    tag = (' (quick preview)' if quick else '') if fig_no >= 6 else ''
    fig.suptitle(f'3D analogue of Fig. {fig_no}: relative error in $U$, {sum_lab} summation{tag}', fontsize=12)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    fn = os.path.join(OUT, f'paper_fig{fig_no}_3D{"_quick" if quick and fig_no >= 6 else ""}.png')
    fig.savefig(fn, dpi=120)
    print('saved', fn)
    return fn


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fig', type=int, nargs='*', default=[2, 3, 4, 5], choices=sorted(FIGS))
    ap.add_argument('--quick', action='store_true', help='low-resolution preview of Figs. 6-8')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--eps-max', type=float, help='Figs. 6-8: eps_max (default 1; the 2D paper uses 2)')
    ap.add_argument('--nx', type=int, help='Figs. 6-8: Nx = Ny (default 64)')
    ap.add_argument('--replot', action='store_true')
    ap.add_argument('--no-show', action='store_true')
    a = ap.parse_args()
    if a.no_show:
        matplotlib.use('Agg')
    for f in a.fig:
        make(f, a.quick, a.workers, a.replot, a.eps_max, a.nx)
    if not a.no_show:
        plt.show()
