"""mms_error_3D.py -- 3D analogue of mms_error.py / mms_error.m (Method of Manufactured Solutions).

Manufactured doubly periodic outgoing solutions
    u_rs = A exp(i p_r x + i q_s y + i gamma^u_rs z),   w_rs = B exp(i p_r x + i q_s y - i gamma^w_rs z)
are fed to the 3D two-layer HOPS/AWE solver; the recovered interface data (U, W, G, J, ubar, wbar)
are compared with the exact values over a grid of (eps, delta) for every frequency window.

Presets (3D resolutions are necessarily smaller than the 2D ones: Nx*Ny unknowns per level)
  'matlab'  : as mms_error.m: N = M = 12, Nx = Ny = Nz = 32, a = b = 4, eps_max = 0.2,
              f = cos(4x)cos(4y)/4 (3D f_s), Pade, absolute errors
  'fig2'    : paper Fig. 2 analogue: N = M = 4, eps_max = 1e-2, a = b = 1, Taylor, relative
  'fig3'    : paper Fig. 3 analogue: N = M = 8, eps_max = 1e-4
  'fig6'    : paper Fig. 6 analogue: N = M = 16, eps_max = 2, Nx = Ny = 64, Nz = 48, a = b = 4, Pade
  'oblique' : alpha = 0.1, beta = 0.2, egg-crate profile, mode (1, 3)
Frequency windows: the 3D bands between consecutive Rayleigh frequencies (hops3d/windows.py); window 1
is omega in [1, sqrt 2] (2D band q = 1 is [1, 2]).  --windows 1 2 3 selects several.

Usage:  python mms_error_3D.py [--preset fig2] [--nm 4] [--epsmax 1e-2] [--nw 1.1] [--windows 1 2]
"""
import argparse
import os
import time

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import hops3d as h3
from hops3d.mms import manufactured_data, exact_values
from hops3d.windows import make_windows
from hops3d.energy import sum_interface_3d
from hops.plotting import safe_log10, parula, matlab_contourf, shared_colorbar

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures_3d', 'mms_error')
SHOW = True

BASE = dict(N=12, M=12, Nx=32, Ny=32, Nz=32, N_delta=40, N_Eps=40, windows=(1,), Eps_Max=0.2, sigma=0.99,
            alpha_bar=0.0, beta_bar=0.0, n_u=1.0, n_w=1.1, Mode=2, A=5.0, B=3.0, rs=(4, 0),
            a=4.0, b=4.0, profile='fs', plot_sum='pade', relative=False, sum_types=(1, 2, 3),
            taylor_full_order=False, solver='coupled', quantities=('G', 'J', 'U', 'W', 'ubar', 'wbar'))
PRESETS = {
    'matlab': dict(),
    'fig2': dict(N=4, M=4, Eps_Max=1e-2, a=1.0, b=1.0, plot_sum='taylor', relative=True, sum_types=(1,)),
    'fig3': dict(N=8, M=8, Eps_Max=1e-4, a=1.0, b=1.0, plot_sum='taylor', relative=True, sum_types=(1,)),
    'fig6': dict(N=16, M=16, Nx=64, Ny=64, Nz=48, Eps_Max=2.0, a=4.0, b=4.0, plot_sum='pade',
                 relative=True, sum_types=(2,), quantities=('U', 'ubar', 'W', 'wbar')),
    'oblique': dict(N=10, M=10, Nx=16, Ny=16, Nz=24, alpha_bar=0.1, beta_bar=0.2, profile='egg', rs=(1, 3),
                    a=1.0, b=1.0, relative=True),
}


def run(cfg, verbose=True):
    c = dict(BASE)
    c.update(cfg)
    N, M = c['N'], c['M']
    Eps = np.linspace(0, c['Eps_Max'], c['N_Eps'])
    wins = make_windows(1.0, 8.0, c['n_u'], c['n_w'], c['alpha_bar'], c['beta_bar'], c['sigma'],
                        Nx=c['Nx'], Ny=c['Ny'])
    names = tuple(c['quantities'])
    want_dno = 'G' in names or 'J' in names
    out = []
    for iw in c['windows']:
        key, omega_bar, dmax, edges = wins[iw - 1]
        delta = np.linspace(-dmax, dmax, c['N_delta'])
        t0 = time.time()
        P = h3.make_problem(c['Nx'], c['Ny'], c['Nz'], N, M, c['n_u'], c['n_w'], omega_bar, c['alpha_bar'],
                            c['beta_bar'], profile=c['profile'], a=c['a'], b=c['b'], Mode=c['Mode'])
        xi_u, nu_u, xi_w, nu_w, zeta, psi = manufactured_data(P, c['A'], c['B'], c['rs'])
        coefs = {}
        if want_dno:
            Su, Sw = P.setups()
            coefs['G'] = h3.dno_tfe_helmholtz_3d(xi_u, Su, N, M)[0]
            coefs['J'] = h3.dno_tfe_helmholtz_3d(xi_w, Sw, N, M)[0]
        coefs['U'], coefs['W'], coefs['ubar'], coefs['wbar'] = h3.SOLVERS_3D[c['solver']](P, zeta, psi)
        if verbose:
            print(f'window {iw} (omega in [{edges[0]:.4f}, {edges[1]:.4f}], omega_bar = {omega_bar:.4f}): '
                  f'HOPS/AWE recursions {time.time() - t0:.1f} s', flush=True)
        t1 = time.time()
        err = {}
        exact_all = [exact_values(P, c['A'], c['B'], c['rs'], Eps, dl) for dl in delta]
        for st in c['sum_types']:
            for nm in names:
                approx = sum_interface_3d(st, coefs[nm], Eps, delta, N, M, c['taylor_full_order'])  # (E, D, Nx, Ny)
                ex = np.stack([e[nm] for e in exact_all], axis=1)                                    # (E, D, Nx, Ny)
                e = np.max(np.abs(ex - approx), axis=(-2, -1))
                if c['relative']:
                    e = e / np.max(np.abs(ex), axis=(-2, -1))
                err[(nm, st)] = e
        if verbose:
            print(f'   summation {time.time() - t1:.1f} s')
        out.append(dict(window=iw, omega=omega_bar * (1 + delta), delta=delta, Eps=Eps, err=err, coefs=coefs,
                        omega_bar=omega_bar))
    return out, c


def plot(out, c, tag='matlab', outdir=OUTDIR):
    from matplotlib.colors import Normalize
    os.makedirs(outdir, exist_ok=True)
    st = 1 if c['plot_sum'] == 'taylor' else 2
    lab = 'Taylor' if st == 1 else 'Pade'
    files = []
    for num, nm in ((1, 'U'), (2, 'ubar')):
        if (nm, st) not in out[0]['err']:
            continue
        fig, ax = plt.subplots(num=num, figsize=(6.5, 5), clear=True)
        Zs = [safe_log10(res['err'][(nm, st)]) for res in out]
        fin = np.concatenate([z[np.isfinite(z)] for z in Zs])
        norm = Normalize(fin.min(), fin.max())
        for res, Z in zip(out, Zs):
            matlab_contourf(ax, res['omega'], res['Eps'], Z, parula, norm)
        shared_colorbar(fig, ax, parula, norm)
        ax.set_xlabel(r'$\omega=\bar\omega(1+\delta)$', fontsize=15)
        ax.set_ylabel(r'$\varepsilon$', fontsize=17)
        sym = 'U' if nm == 'U' else r'\bar{u}'
        ax.set_title(('Relative Error' if c['relative'] else 'Error') +
                     f' in ${sym}$ (3D, {lab}, N=M={c["N"]}, $N_x$=$N_y$={c["Nx"]})', fontsize=12)
        fig.tight_layout()
        fn = os.path.join(outdir, f'mms3D_{tag}_{nm}_{lab.lower()}.png')
        fig.savefig(fn, dpi=130)
        files.append(fn)
    return files


def summary(out, c):
    print('max error over the (eps, delta) grid:')
    for res in out:
        for (nm, st), e in sorted(res['err'].items(), key=lambda t: (t[0][1], t[0][0])):
            lab = {1: 'Taylor', 2: 'Pade', 3: 'Pade-safe'}[st]
            print(f'  window {res["window"]}  {lab:9s} {nm:5s}  max={np.nanmax(e):.3e}  median={np.nanmedian(e):.3e}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--preset', default='matlab', choices=list(PRESETS))
    ap.add_argument('--nm', type=int, help='override N = M')
    ap.add_argument('--epsmax', type=float)
    ap.add_argument('--nw', type=float)
    ap.add_argument('--nx', type=int, help='override Nx = Ny')
    ap.add_argument('--nz', type=int)
    ap.add_argument('--neps', type=int)
    ap.add_argument('--ndelta', type=int)
    ap.add_argument('--windows', type=int, nargs='*', help='frequency windows (1 = [1, sqrt 2], 2 = [sqrt 2, 2], ...)')
    ap.add_argument('--profile')
    ap.add_argument('--alpha', type=float)
    ap.add_argument('--beta', type=float)
    ap.add_argument('--solver', choices=list(h3.SOLVERS_3D))
    ap.add_argument('--no-show', action='store_true')
    args = ap.parse_args()
    if args.no_show:
        matplotlib.use('Agg')
    cfg = dict(PRESETS[args.preset])
    if args.nm: cfg.update(N=args.nm, M=args.nm)
    if args.epsmax: cfg['Eps_Max'] = args.epsmax
    if args.nw: cfg['n_w'] = args.nw
    if args.nx: cfg.update(Nx=args.nx, Ny=args.nx)
    if args.nz: cfg['Nz'] = args.nz
    if args.neps: cfg['N_Eps'] = args.neps
    if args.ndelta: cfg['N_delta'] = args.ndelta
    if args.windows: cfg['windows'] = tuple(args.windows)
    if args.profile: cfg['profile'] = args.profile
    if args.alpha is not None: cfg['alpha_bar'] = args.alpha
    if args.beta is not None: cfg['beta_bar'] = args.beta
    if args.solver: cfg['solver'] = args.solver
    t0 = time.time()
    out, c = run(cfg)
    summary(out, c)
    tag = f'{args.preset}_NM{c["N"]}_eps{c["Eps_Max"]:g}_nw{c["n_w"]:g}'
    os.makedirs(OUTDIR, exist_ok=True)
    np.savez_compressed(os.path.join(OUTDIR, f'mms3D_{tag}.npz'),
                        **{f'w{o["window"]}_{nm}_st{st}': e for o in out for (nm, st), e in o['err'].items()},
                        **{f'w{o["window"]}_omega': o['omega'] for o in out}, Eps=out[0]['Eps'])
    for fn in plot(out, c, tag):
        print('saved', fn)
    print(f'total {time.time() - t0:.1f} s')
    if SHOW and not args.no_show:
        plt.show()
