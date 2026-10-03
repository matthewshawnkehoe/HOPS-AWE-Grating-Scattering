"""mms_error.py -- Python port of mms_error.m (Method of Manufactured Solutions).

Manufactured outgoing solutions  u_r = A e^{i p_r x + i gamma^u_r z},
w_r = B e^{i p_r x - i gamma^w_r z}  (paper Sec. 6.5) are fed to the HOPS/AWE
two-layer solver and the recovered interface data (U, W, G, J, ubar, wbar) are
compared with the exact values over a grid of (eps, delta).

Presets
  'matlab'   : exactly the parameters of mms_error.m in src.zip (default)
               N=M=16, Nx=Nz=32, a=b=4, eps_max=0.2, f=cos(4x)/4, Pade plots.
  'fig2'     : paper Fig. 2 (a)  N=M=4,  eps_max=1e-2, a=b=1, Taylor, relative error
  'fig3'     : paper Fig. 3 (b)  N=M=8,  eps_max=1e-4           (use --nm/--epsmax to vary)
  'fig6'     : paper Fig. 6 (a)  N=M=20, eps_max=2, Nx=256, Nz=128, Pade  (slow!)
Usage:  python mms_error.py [--preset fig2] [--nm 4] [--epsmax 1e-2] [--nw 1.1]
"""
import argparse
import os
import time

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from hops import (cheb, setup_2d, csqrt, setup_xi_u_nu_u_n_m, setup_xi_w_nu_w_n_m,
                  field_tfe_helmholtz_m_and_n, field_tfe_helmholtz_m_and_n_lf,
                  dno_tfe_helmholtz_m_and_n, dno_tfe_helmholtz_m_and_n_lf,
                  two_layer_solve_fast, fcn_sum, fourier_repr_lipschitz, fourier_repr_rough)
from hops.plotting import safe_log10

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'mms_error')
SHOW = True

BASE = dict(N=16, M=16, Nx=32, Nz=32, N_delta=100, N_Eps=100, qq=(1,), Eps_Max=0.2, sigma=0.99,
            alpha_bar=0.0, d=2 * np.pi, c_0=1.0, n_u=1.0, n_w=1.1, Mode=2, A=5.0, B=3.0, r=4,
            a=4.0, b=4.0, profile='cos4x_over4', plot_sum='pade', relative=False,
            sum_types=(1, 2, 3), taylor_full_order=False,
            solver='coupled', quantities=('G', 'J', 'U', 'W', 'ubar', 'wbar'))
PRESETS = {
    'matlab': dict(),
    'fig2': dict(N=4, M=4, Eps_Max=1e-2, a=1.0, b=1.0, plot_sum='taylor', relative=True, sum_types=(1,)),
    'fig3': dict(N=8, M=8, Eps_Max=1e-4, a=1.0, b=1.0, plot_sum='taylor', relative=True, sum_types=(1,)),
    'fig6': dict(N=20, M=20, Nx=256, Nz=128, Eps_Max=2.0, a=4.0, b=4.0, plot_sum='pade',
                 relative=True, sum_types=(2,)),
}


def profile_fn(name, xx):
    if name == 'cos4x_over4':
        return 0.25 * np.cos(4 * xx), -np.sin(4 * xx)
    if name == 'lipschitz':
        return fourier_repr_lipschitz(120, xx)
    if name == 'rough':
        return fourier_repr_rough(120, xx)
    raise ValueError(name)


def run(cfg, verbose=True):
    c = dict(BASE); c.update(cfg)
    N, M, Nx, Nz = c['N'], c['M'], c['Nx'], c['Nz']
    A, B, r, a, b = c['A'], c['B'], c['r'], c['a'], c['b']
    alpha_bar, d, c_0, n_u, n_w = c['alpha_bar'], c['d'], c['c_0'], c['n_u'], c['n_w']
    N_Eps, N_delta, sigma = c['N_Eps'], c['N_delta'], c['sigma']
    Eps = np.linspace(0, c['Eps_Max'], N_Eps)
    identy = np.eye(Nz + 1)
    Dz, _ = cheb(Nz)
    xx = (d / Nx) * np.arange(Nx)
    f, f_x = profile_fn(c['profile'], xx)
    names = tuple(c['quantities'])
    want_dno = 'G' in names or 'J' in names
    out = []
    for q in c['qq']:
        delta = np.linspace(-sigma / (2 * q + 1), sigma / (2 * q + 1), N_delta)
        omega_bar = q + 0.5
        omega = (1 + delta) * omega_bar
        k_u_bar = n_u * omega_bar / c_0
        gamma_u_bar = csqrt(k_u_bar ** 2 - alpha_bar ** 2)
        xx, pp, abp, gubp, _, _ = setup_2d(Nx, d, alpha_bar, gamma_u_bar)
        k_w_bar = n_w * omega_bar / c_0
        gamma_w_bar = csqrt(k_w_bar ** 2 - alpha_bar ** 2)
        xx, pp, abp, gwbp, _, _ = setup_2d(Nx, d, alpha_bar, gamma_w_bar)
        pp_r = pp[r]
        alpha_bar_r = abp[r]

        t0 = time.time()
        xi_u, nu_u = setup_xi_u_nu_u_n_m(A, r, xx, pp, abp, gubp, f, f_x, Nx, N, M)
        xi_w, nu_w = setup_xi_w_nu_w_n_m(B, r, xx, pp, abp, gwbp, f, f_x, Nx, N, M)
        G_n_m = J_n_m = None
        if want_dno:     # (skipped for the large Figs. 6-8 runs: needs the full volume fields)
            u_n_m = field_tfe_helmholtz_m_and_n(xi_u, f, pp, gubp, alpha_bar, gamma_u_bar, Dz, a, Nx, Nz, N, M, identy, abp)
            G_n_m = dno_tfe_helmholtz_m_and_n(u_n_m, f, pp, Dz, a, Nx, Nz, N, M)
            w_n_m = field_tfe_helmholtz_m_and_n_lf(xi_w, f, pp, gwbp, alpha_bar, gamma_w_bar, Dz, b, Nx, Nz, N, M, identy, abp)
            J_n_m = dno_tfe_helmholtz_m_and_n_lf(w_n_m, f, pp, Dz, b, Nx, Nz, N, M)
            del u_n_m, w_n_m
        tau2 = 1.0 if c['Mode'] == 1 else (n_u / n_w) ** 2
        zeta = xi_u - xi_w
        psi = -nu_u - tau2 * nu_w
        # physical interface data for oblique incidence (hops/config.py item 3): the manufactured
        # nu's use the phase-removed d_x, the true normal derivatives carry + i alpha g_x.
        from hops.two_layer import phase_correction
        for _n in range(N + 1):
            for _m in range(M + 1):
                psi[:, _m, _n] = psi[:, _m, _n] - phase_correction(xi_u, xi_w, _m, _n, alpha_bar, f_x, tau2)
        from hops.operators import two_layer_solve_operator, two_layer_solve_auto, two_layer_solve_lean
        from hops.coupled import two_layer_solve_coupled
        solve = {'coupled': two_layer_solve_coupled, 'operator': two_layer_solve_operator, 'fast': two_layer_solve_fast,
                 'auto': two_layer_solve_auto, 'lean': two_layer_solve_lean}[c['solver']]
        U_n_m, W_n_m, ubar_n_m, wbar_n_m = solve(
            tau2, zeta, psi, gubp, gwbp, N, Nx, f, f_x, pp, alpha_bar, gamma_u_bar, gamma_w_bar,
            Dz, a, b, Nz, M, identy, abp)
        if verbose:
            print(f'q = {q}: HOPS/AWE recursions {time.time() - t0:.1f} s', flush=True)
        coefs = dict(G=G_n_m, J=J_n_m, U=U_n_m, W=W_n_m, ubar=ubar_n_m, wbar=wbar_n_m)

        err = {(nm, st): np.zeros((N_Eps, N_delta)) for nm in names for st in c['sum_types']}
        E = Eps[:, None]
        for ell in range(N_delta):
            dl = delta[ell]
            alpha_r = alpha_bar_r + dl * alpha_bar
            k_u = (1 + dl) * k_u_bar
            k_w = (1 + dl) * k_w_bar
            g_ur = csqrt(k_u ** 2 - alpha_r ** 2)
            g_wr = csqrt(k_w ** 2 - alpha_r ** 2)
            ph = np.exp(1j * pp_r * xx)[None, :]
            xi_u_r = A * ph * np.exp(1j * g_ur * E * f)
            nu_u_r = (-1j * g_ur + 1j * pp_r * E * f_x) * xi_u_r
            xi_w_r = B * ph * np.exp(-1j * g_wr * E * f)
            nu_w_r = (-1j * g_wr - 1j * pp_r * E * f_x) * xi_w_r      # -d_N w sign, as MATLAB
            ubar = A * ph * np.exp(1j * g_ur * a) * np.ones_like(E)
            wbar = B * ph * np.exp(-1j * g_wr * (-b)) * np.ones_like(E)
            exact = dict(G=nu_u_r, J=nu_w_r, U=xi_u_r, W=xi_w_r, ubar=ubar, wbar=wbar)
            for st in c['sum_types']:
                for nm in names:
                    approx = fcn_sum(st, coefs[nm], Eps, dl, Nx, N, M, c['taylor_full_order'])       # (N_Eps, Nx)
                    e = np.max(np.abs(exact[nm] - approx), axis=-1)
                    if c['relative']:
                        e = e / np.max(np.abs(exact[nm]), axis=-1)
                    err[(nm, st)][:, ell] = e
        out.append(dict(q=q, omega=omega, delta=delta, Eps=Eps, err=err, coefs=coefs))
    return out, c


def plot(out, c, tag='matlab', outdir=OUTDIR):
    os.makedirs(outdir, exist_ok=True)
    st = 1 if c['plot_sum'] == 'taylor' else 2
    lab = 'Taylor' if st == 1 else 'Pade'
    files = []
    for num, nm in ((1, 'U'), (2, 'ubar')):
        fig, ax = plt.subplots(num=num, figsize=(6.5, 5), clear=True)
        from hops.plotting import parula, matlab_contourf, shared_colorbar
        from matplotlib.colors import Normalize
        Zs = [safe_log10(res['err'][(nm, st)]) for res in out]
        fin = np.concatenate([z[np.isfinite(z)] for z in Zs])
        norm = Normalize(fin.min(), fin.max())
        for res, Z in zip(out, Zs):
            matlab_contourf(ax, res['omega'], res['Eps'], Z, parula, norm)
        shared_colorbar(fig, ax, parula, norm)
        ax.set_xlabel(r'$\omega=\omega_1(1+\delta)$', fontsize=15)
        ax.set_ylabel(r'$\varepsilon$', fontsize=17)
        sym = 'U' if nm == 'U' else r'\bar{u}'
        ax.set_title(('Relative Error' if c['relative'] else 'Error') +
                     f' in ${sym}$ ({lab}, N=M={c["N"]})', fontsize=13)
        fig.tight_layout()
        fn = os.path.join(outdir, f'mms_{tag}_{nm}_{lab.lower()}.png')
        fig.savefig(fn, dpi=130)
        files.append(fn)
    return files


def summary(out, c):
    print('max error over the (eps, delta) grid:')
    for res in out:
        for (nm, st), e in sorted(res['err'].items(), key=lambda t: (t[0][1], t[0][0])):
            lab = {1: 'Taylor', 2: 'Pade', 3: 'Pade-safe'}[st]
            print(f'  q={res["q"]}  {lab:9s} {nm:5s}  max={np.nanmax(e):.3e}  median={np.nanmedian(e):.3e}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--preset', default='matlab', choices=list(PRESETS))
    ap.add_argument('--nm', type=int, help='override N = M')
    ap.add_argument('--epsmax', type=float)
    ap.add_argument('--nw', type=float)
    ap.add_argument('--neps', type=int)
    ap.add_argument('--ndelta', type=int)
    ap.add_argument('--solver', choices=['coupled', 'operator', 'fast', 'auto', 'lean'])
    ap.add_argument('--no-show', action='store_true')
    args = ap.parse_args()
    if args.no_show:
        matplotlib.use('Agg')
    cfg = dict(PRESETS[args.preset])
    if args.nm: cfg.update(N=args.nm, M=args.nm)
    if args.epsmax: cfg['Eps_Max'] = args.epsmax
    if args.nw: cfg['n_w'] = args.nw
    if args.neps: cfg['N_Eps'] = args.neps
    if args.ndelta: cfg['N_delta'] = args.ndelta
    if args.solver: cfg['solver'] = args.solver
    t0 = time.time()
    out, c = run(cfg)
    summary(out, c)
    tag = f'{args.preset}_NM{c["N"]}_eps{c["Eps_Max"]:g}_nw{c["n_w"]:g}'
    os.makedirs(OUTDIR, exist_ok=True)
    np.savez_compressed(os.path.join(OUTDIR, f'mms_{tag}.npz'),
                        **{f'q{o["q"]}_{nm}_st{st}': e for o in out for (nm, st), e in o['err'].items()},
                        **{f'q{o["q"]}_omega': o['omega'] for o in out}, Eps=out[0]['Eps'])
    for fn in plot(out, c, tag):
        print('saved', fn)
    print(f'total {time.time() - t0:.1f} s')
    if SHOW and not args.no_show:
        plt.show()
