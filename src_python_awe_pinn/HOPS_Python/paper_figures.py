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
Every panel is cached in figures/paper/paper_cache/*.npz; --replot redraws from the cache.

Figs. 9, 10, 14 (reflectivity maps): python paper_figures.py --fig 9 10 14   (~1 min each)
  Fig. 9 is also compared with the ORIGINAL MATLAB code run under Octave (full 6-band, 100 x 100 map,
  reference_data/ref_refl_dielectric_full.mat) -> figures/paper/paper_fig9_matlab_vs_python.png.

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

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'paper')
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
    cfg = dict(FIGS[fig_no][i], relative=True, quantities=('U',), solver='coupled')
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


# ----------------------------------------------------------------------------
# Figs. 9, 10, 14: reflectivity maps R and energy defect D (refl_map.py scenarios)
# ----------------------------------------------------------------------------
REFL_FIGS = {
    9: [('dielectric', 'R', '(a) Reflectivity Map'), ('dielectric', 'D', '(b) Energy Defect')],
    10: [('silver', 'R', '(a) Silver'), ('gold', 'R', '(b) Gold')],
    14: [('dielectric_alpha', 'R', '(a) Reflectivity Map'), ('dielectric_alpha', 'D', '(b) Energy Defect')],
}
MATLAB_FULL = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference_data',
                           'ref_refl_dielectric_full.mat')


def _refl_panel(ax, fig, results, what, title, relative=True):
    from hops.plotting import matlab_contourf, shared_colorbar
    getZ = (lambda r: safe_log10(r['ee'])) if what == 'D' else \
        (lambda r: np.real(r['RR']) if relative else np.real(r['ru']))
    Zs = [getZ(r) for r in results]
    allz = np.concatenate([z[np.isfinite(z)].ravel() for z in Zs])
    vmin, vmax = allz.min(), allz.max()
    if what == 'R' and vmax > 1.0 and np.mean(allz > 1.0 + 1e-9) < 1e-3:
        vmax, Zs = 1.0, [np.minimum(z, 1.0) for z in Zs]
    norm = Normalize(vmin, vmax)
    for r, Z in zip(results, Zs):
        matlab_contourf(ax, r['lam'], r['Eps'], Z, 'hot', norm, step=1.0 if (what == 'D' and vmax - vmin >= 4) else None)
    shared_colorbar(fig, ax, 'hot', norm)
    ax.set_xlabel(r'$\lambda$' + '\n' + title)
    ax.set_ylabel(r'$\varepsilon$')
    ax.set_title('$D$' if what == 'D' else '$R$')
    return vmin, vmax


def make_refl(fig_no, workers=2):
    """Paper Figs. 9, 10, 14 (N = M = 16 Taylor / 15 Pade, Nx = Nz = 32, 100 x 100 per band).
    Fig. 9 additionally: the same maps computed by the ORIGINAL MATLAB code (run under Octave,
    reference_data/ref_refl_dielectric_full.mat) drawn with the same renderer, and their difference."""
    import refl_map as rm
    cache = {}
    for sc, _, _ in REFL_FIGS[fig_no]:
        if sc not in cache:
            cache[sc] = rm.run(sc, verbose=False, workers=workers)[0]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, (sc, what, title) in zip(axs, REFL_FIGS[fig_no]):
        lo, hi = _refl_panel(ax, fig, cache[sc], what, title)
        print(f'Fig {fig_no} {title}: {what} in [{lo:.3g}, {hi:.3g}]')
    fig.suptitle(f'Fig. {fig_no} (Python)', fontsize=12)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    fn = os.path.join(OUT, f'paper_fig{fig_no}.png')
    fig.savefig(fn, dpi=120)
    print('saved', fn)
    if fig_no == 9 and os.path.exists(MATLAB_FULL):
        fn2 = compare_fig9_with_matlab(cache['dielectric'])
        return fn, fn2
    return fn


def load_matlab_fig9():
    """Octave run of the unmodified MATLAB code (test_scenarios.m settings, 6 bands, 100 x 100)."""
    import scipy.io as sio
    o = sio.loadmat(MATLAB_FULL, squeeze_me=True, struct_as_record=False)['out']
    res = []
    for k in sorted(o._fieldnames, key=lambda s: int(s[1:])):
        m = getattr(o, k)
        q = int(k[1:])
        omega = (q + 0.5) * (1 + np.asarray(m.delta))
        res.append(dict(q=q, lam=2 * np.pi / omega, Eps=np.asarray(m.Eps), ee=m.ee_t, ru=m.ru_t,
                        RR=np.real(m.ru_t) / np.real(m.ru_flat)))
    return res


def compare_fig9_with_matlab(py_results):
    mat = load_matlab_fig9()
    py = {r['q']: r for r in py_results}
    fig, axs = plt.subplots(2, 3, figsize=(17, 8.5))
    for row, (res, lab) in enumerate(((mat, 'MATLAB code (Octave)'), ([py[r['q']] for r in mat], 'Python'))):
        _refl_panel(axs[row, 0], fig, res, 'R', f'{lab}: Reflectivity Map')
        _refl_panel(axs[row, 1], fig, res, 'D', f'{lab}: Energy Defect')
    # differences
    from hops.plotting import matlab_contourf, shared_colorbar
    dR = [np.abs(np.real(py[m['q']]['ru']) - np.real(m['ru'])) for m in mat]
    dD = [np.abs(py[m['q']]['ee'] - m['ee']) for m in mat]
    for ax, dl, t in ((axs[0, 2], dR, r'$|R_{Python} - R_{MATLAB}|$'), (axs[1, 2], dD, r'$|D_{Python} - D_{MATLAB}|$')):
        Zs = [safe_log10(np.maximum(d, 1e-17)) for d in dl]
        norm = Normalize(-17, max(-12, max(np.nanmax(z) for z in Zs)))
        for m, Z in zip(mat, Zs):
            ax.pcolormesh(m['lam'], m['Eps'], Z, cmap=parula, norm=norm, shading='auto')
        shared_colorbar(fig, ax, parula, norm)
        ax.set_title(t + '  (log10)')
        ax.set_xlabel(r'$\lambda$')
        ax.set_ylabel(r'$\varepsilon$')
    mx = max(np.max(d) for d in dD)
    fig.suptitle(f'Paper Fig. 9: original MATLAB code vs Python port (same data -> same picture); '
                 f'max |D_py - D_mat| = {mx:.1e}', fontsize=12)
    fig.tight_layout()
    fn = os.path.join(OUT, 'paper_fig9_matlab_vs_python.png')
    fig.savefig(fn, dpi=110)
    print('saved', fn, f'(max |D_py - D_mat| = {mx:.2e})')
    return fn


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fig', type=int, nargs='*', default=[2, 3, 4, 5], choices=sorted(FIGS) + sorted(REFL_FIGS),
                    help='2-8: MMS error figures; 9, 10, 14: reflectivity map / energy defect figures')
    ap.add_argument('--quick', action='store_true', help='low-resolution preview of Figs. 6-8')
    ap.add_argument('--res', choices=['reduced', 'paper'], default='reduced',
                    help="Figs. 6-8: 'reduced' (default, converged, 2-6x faster) or 'paper' (Nx=256/1024, Nz=128)")
    ap.add_argument('--workers', type=int, default=1, help='panels computed in parallel processes')
    ap.add_argument('--replot', action='store_true', help='only redraw from figures/paper/paper_cache')
    ap.add_argument('--no-show', action='store_true')
    a = ap.parse_args()
    if a.no_show:
        matplotlib.use('Agg')
    for f in a.fig:
        if f in REFL_FIGS:
            make_refl(f, workers=max(2, a.workers))
        else:
            make(f, a.quick, a.workers, a.replot, a.res)
    if not a.no_show:
        plt.show()
