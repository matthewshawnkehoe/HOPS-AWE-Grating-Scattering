"""compare_maps.py -- HOPS/AWE vs the HOPS/AWE + least-squares-PINN hybrid (physics-informed summation)
on reflectivity maps R(eps, lambda) and energy defects D, for several refl_map.py scenarios.

For every scenario (one refl_map.py band):
  AWE        the map of refl_map.py (joint (eps, delta) expansion, Taylor or Pade as in the paper)
  PI         the hybrid: the same HOPS/AWE coefficient fields, weights from the PINN least-squares
             problem at every (eps, omega)                                         (hybrid.PISum.map)
  PI-adapt   the hybrid only where the a-posteriori indicator (PINN loss of the AWE Taylor sum) > tol
  PI+RF      (some scenarios) PI enriched with random Fourier x sin features (LSQ-PINN of PINN_HOPS)
  truth      HOPS at every grid point separately (delta = 0: no frequency expansion), N = 24, Pade,
             finer grid (Nx = 64 or 128, Nz = 48) -- cached in results/<name>/truth.npz
Outputs results/<name>/{summary.json, maps.png, indicator.png, maps.npz}, results/summary.md.

    python compare_maps.py                           # all scenarios (~1 h on 2 cores)
    python compare_maps.py --only dielectric_q1 gold_q1
"""
import argparse
import json
import os
import time
import warnings

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from hybrid import AWEBand, PISum, RandomFeatures
from pinn_hops import Grating2D
from pinn_hops.hops_reference import hops_point

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')

SCENARIOS = {
    # name: (refl_map scenario, band q, overrides, grid, truth Nx, enrich)
    'dielectric_q1': ('dielectric', 1, {}, 21, 64, False),        # paper Fig. 9 (Taylor, N = M = 16)
    'dielectric_q2': ('dielectric', 2, {}, 21, 64, False),
    'dielectric_q3': ('dielectric', 3, {}, 21, 64, False),
    'te_dielectric_q1': ('thesis29', 1, {}, 21, 64, False),       # thesis Fig. 29 (TE, Taylor, N = M = 15)
    'silver_q1': ('silver', 1, {}, 15, 128, False),               # paper Fig. 10a (Pade, cos 4x)
    'gold_q1': ('gold', 1, {}, 15, 128, False),                   # paper Fig. 10b
    'gold_q3': ('gold', 3, {}, 15, 128, False),
    'silver_q1_nx64': ('silver', 1, dict(Nx=64), 15, 128, False),    # as the paper, but a resolved x-grid
    'gold_q1_nx64': ('gold', 1, dict(Nx=64), 15, 128, False),
    'gold_q3_nx64': ('gold', 3, dict(Nx=64), 15, 128, False),
    'high_index_q1': ('dielectric', 1, dict(n_w=1.8), 21, 64, True),    # 2 lower-layer Wood anomalies in the band
    'large_eps_q1': ('dielectric', 1, dict(eps_max=0.6), 21, 64, False),  # 3x the paper's eps range
}
ADAPT_TOL = 1e-18      # indicator threshold for the adaptive hybrid
COMPRESS = 1e-13       # POD of the HOPS/AWE basis (drops numerically dependent directions)


def truth_map(B, name, Nx):
    fn = os.path.join(OUT, name, 'truth.npz')
    if os.path.exists(fn):
        d = np.load(fn)
        return d['R'], d['T'], d['D'], float(d['time'])
    R, T, D = (np.zeros((len(B.Eps), len(B.omega))) for _ in range(3))
    t0 = time.time()
    prof = B.info['profile']
    for i, e in enumerate(B.Eps):
        for j, o in enumerate(B.omega):
            G = Grating2D(n_u=B.n_u, n_w=B.n_w, omega=float(o), eps=float(e), profile=prof,
                          alpha=B.alpha_bar * o / B.omega_bar, a=B.a, b=B.b, mode=B.mode)
            r = hops_point(G, N=24, Nx=Nx, Nz=48, summation='pade', fields=False)
            R[i, j], T[i, j], D[i, j] = r['R'], r['T'], r['D']
    t = time.time() - t0
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    np.savez(fn, R=R, T=T, D=D, time=t)
    return R, T, D, t


def stats(err):
    err = np.abs(err)
    return dict(median=float(np.median(err)), max=float(err.max()), p90=float(np.percentile(err, 90)))


def run(name, verbose=True):
    sc, q, over, n, tNx, enrich = SCENARIOS[name]
    B = AWEBand(sc, q=q, N_Eps=n, N_delta=n, **over)
    if verbose:
        print(f'== {name}: {sc} q={q} {over}  n_w={B.n_w:.4g} {B.mode} N={B.N} M={B.M} Nx={B.Nx} '
              f'{"Taylor" if B.info["Taylor"] else "Pade"}; AWE band {B.t_hops:.1f} s', flush=True)
    Rt, Tt, Dt, t_truth = truth_map(B, name, tNx)
    S = PISum(B, compress=COMPRESS)
    pi = S.map()
    ad = S.map(adaptive_tol=ADAPT_TOL)
    res = dict(AWE=dict(R=B.R_awe, D=B.D_awe, time=B.t_hops), PI=pi, PI_adapt=ad)
    if enrich:
        Se = PISum(B, enrich=RandomFeatures(B, K=min(12, B.Nx // 2 - 2), nz=8), compress=COMPRESS)
        res['PI+RF'] = Se.map()
    lossless = abs(B.n_w.imag) < 1e-9 and abs(B.n_u.imag) < 1e-9
    summ = dict(name=name, scenario=sc, q=q, overrides={k: str(v) for k, v in over.items()}, n_w=str(B.n_w),
                mode=B.mode, N=B.N, M=B.M, Nx=B.Nx, Nz=B.Nz, summation='Taylor' if B.info['Taylor'] else 'Pade',
                grid=f'{n}x{n}', lossless=lossless, t_truth_s=t_truth, t_awe_s=B.t_hops,
                truth_D_median=float(np.median(np.abs(Dt))))
    for k, v in res.items():
        summ[k] = dict(R_err=stats(v['R'] - Rt), D_err=stats(v['D'] - Dt),
                       D_median_log10=float(np.log10(np.median(np.abs(v['D'])) + 1e-300)),
                       time_s=float(v['time']))
        if 'n_solved' in v:
            summ[k]['n_solved'] = int(v['n_solved'])
    # indicator quality: correlation of log indicator with log AWE-Taylor error
    from scipy.stats import spearmanr
    summ['indicator_corr'] = float(spearmanr(pi['loss_awe_taylor'].ravel(), np.abs(pi['R_taylor'] - Rt).ravel())[0])
    d = os.path.join(OUT, name)
    os.makedirs(d, exist_ok=True)
    json.dump(summ, open(os.path.join(d, 'summary.json'), 'w'), indent=1)
    np.savez_compressed(os.path.join(d, 'maps.npz'), Eps=B.Eps, omega=B.omega, R_true=Rt, D_true=Dt,
                        R_flat=B.R_flat, **{f'{k}_{q_}': np.asarray(v[q_]) for k, v in res.items()
                                            for q_ in ('R', 'D')},
                        indicator=pi['loss_awe_taylor'], R_taylor=pi['R_taylor'])
    figure(name, B, Rt, Dt, res, pi, summ)
    if verbose:
        for k in res:
            s = summ[k]
            print(f"   {k:9s} R err median {s['R_err']['median']:.1e} max {s['R_err']['max']:.1e} | "
                  f"D err median {s['D_err']['median']:.1e} max {s['D_err']['max']:.1e} | time {s['time_s']:.1f} s"
                  + (f" ({s['n_solved']} LSQ solves)" if 'n_solved' in s else ''), flush=True)
    return summ


def figure(name, B, Rt, Dt, res, pi, summ):
    lam = 2 * np.pi / B.omega
    Eps = B.Eps
    fig, axs = plt.subplots(2, 4, figsize=(21, 8.8))
    RR = [B.R_awe / B.R_flat, pi['R'] / B.R_flat]
    lo, hi = min(np.nanmin(z) for z in RR), max(np.nanmax(z) for z in RR)
    lv = np.linspace(lo, hi, 15)
    for ax, Z, t in ((axs[0, 0], RR[0], 'HOPS/AWE (refl_map.py)'), (axs[0, 1], RR[1], 'hybrid (physics-informed sum)')):
        m = ax.contourf(lam, Eps, Z, levels=lv, cmap='hot')
        ax.contour(lam, Eps, Z, levels=lv[1:-1], colors='k', linewidths=0.4)
        fig.colorbar(m, ax=ax)
        ax.set_title(r'$R/R_{flat}$: ' + t)
    E = [np.log10(np.abs(B.R_awe - Rt) + 1e-17), np.log10(np.abs(pi['R'] - Rt) + 1e-17)]
    for ax, Z, t in ((axs[0, 2], E[0], 'AWE'), (axs[0, 3], E[1], 'hybrid')):
        m = ax.pcolormesh(lam, Eps, Z, cmap='viridis', vmin=-17, vmax=max(-2, E[0].max()), shading='auto')
        fig.colorbar(m, ax=ax)
        ax.set_title(rf'$\log_{{10}}|R - R_{{true}}|$: {t}')
    lab = 'log10|D|' if summ['lossless'] else 'log10|D - D_true|'
    Ds = ([np.log10(np.abs(B.D_awe) + 1e-17), np.log10(np.abs(pi['D']) + 1e-17)] if summ['lossless'] else
          [np.log10(np.abs(B.D_awe - Dt) + 1e-17), np.log10(np.abs(pi['D'] - Dt) + 1e-17)])
    for ax, Z, t in ((axs[1, 0], Ds[0], 'AWE'), (axs[1, 1], Ds[1], 'hybrid')):
        m = ax.pcolormesh(lam, Eps, Z, cmap='hot', vmin=-17, vmax=max(-2, Ds[0].max()), shading='auto')
        fig.colorbar(m, ax=ax)
        ax.set_title(f'{lab}: {t}')
    m = axs[1, 2].pcolormesh(lam, Eps, np.log10(pi['loss_awe_taylor'] + 1e-30), cmap='magma', shading='auto')
    fig.colorbar(m, ax=axs[1, 2])
    axs[1, 2].set_title('indicator: log10 PINN loss of the AWE Taylor sum')
    ax = axs[1, 3]
    ax.loglog(pi['loss_awe_taylor'].ravel() + 1e-30, np.abs(pi['R_taylor'] - Rt).ravel() + 1e-17, '.', ms=3,
              label='AWE Taylor sum (full order)')
    ax.loglog(pi['loss'].ravel() + 1e-30, np.abs(pi['R'] - Rt).ravel() + 1e-17, '.', ms=3, label='hybrid')
    ax.axvline(ADAPT_TOL, color='k', ls=':', lw=0.8)
    ax.set_xlabel('PINN loss (computable, no reference)')
    ax.set_ylabel('|R - R_true|')
    ax.set_title(f"indicator vs true error (Spearman {summ['indicator_corr']:.2f})")
    ax.legend(fontsize=8)
    for a_ in axs.ravel()[:7]:
        a_.set_xlabel(r'$\lambda$')
        a_.set_ylabel(r'$\varepsilon$')
    s = summ
    fig.suptitle(f"{name}: n_w = {B.n_w:.4g}, {B.mode}, {s['summation']} N={B.N} M={B.M}, Nx={B.Nx} | "
                 f"R err median/max: AWE {s['AWE']['R_err']['median']:.1e}/{s['AWE']['R_err']['max']:.1e}, "
                 f"hybrid {s['PI']['R_err']['median']:.1e}/{s['PI']['R_err']['max']:.1e} | time AWE {s['AWE']['time_s']:.1f} s, "
                 f"hybrid {s['PI']['time_s']:.0f} s, adaptive {s['PI_adapt']['time_s']:.0f} s", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name, 'maps.png'), dpi=95)
    plt.close(fig)


def table():
    rows = ['| scenario | AWE sum | R err AWE (median / max) | R err hybrid (median / max) | D err AWE (max) | '
            'D err hybrid (max) | time AWE / hybrid / adaptive (LSQ solves) |', '|---|---|---|---|---|---|---|']
    for name in SCENARIOS:
        fn = os.path.join(OUT, name, 'summary.json')
        if not os.path.exists(fn):
            continue
        s = json.load(open(fn))
        a, p, ad = s['AWE'], s['PI'], s['PI_adapt']
        extra = ''
        if 'PI+RF' in s:
            extra = f"<br>hybrid + random features: {s['PI+RF']['R_err']['median']:.1e} / {s['PI+RF']['R_err']['max']:.1e}"
        rows.append(f"| {name} (n_w = {s['n_w']}, {s['mode']}) | {s['summation']} N={s['N']},M={s['M']} | "
                    f"{a['R_err']['median']:.1e} / {a['R_err']['max']:.1e} | {p['R_err']['median']:.1e} / "
                    f"{p['R_err']['max']:.1e}{extra} | {a['D_err']['max']:.1e} | {p['D_err']['max']:.1e} | "
                    f"{a['time_s']:.1f} s / {p['time_s']:.0f} s / {ad['time_s']:.0f} s ({ad.get('n_solved', 0)}) |")
    open(os.path.join(OUT, 'summary.md'), 'w').write('\n'.join(rows) + '\n')
    print('\n'.join(rows))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--only', nargs='*', default=list(SCENARIOS))
    a = ap.parse_args()
    matplotlib.use('Agg')
    warnings.filterwarnings('ignore')
    for nm in a.only:
        if os.path.exists(os.path.join(OUT, nm, 'summary.json')):
            print('skip', nm)
            continue
        run(nm)
    table()
