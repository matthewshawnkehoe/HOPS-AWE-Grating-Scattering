"""run_scenarios.py -- every refl_map.py scenario of Hops_2D_AWE with HOPS/AWE and with HOPS/AWE + PINN.

For each scenario: the AWE maps exactly as refl_map.py computes them, the hybrid maps (physics-informed
summation of the same coefficient fields, adaptive), the indicator map, and figures:
  figures/<scenario>/refl_map_<scenario>_awe_{R,D}.png      refl_map.py rendering of HOPS/AWE
  figures/<scenario>/refl_map_<scenario>_hybrid_{R,D}.png   the same rendering of HOPS/AWE + PINN
  figures/<scenario>/compare.png                            |R_hyb - R_awe|, indicator, D (both), solved points
and results/scenarios/<scenario>.npz, results/scenarios/<scenario>.json (statistics).

Grid policy (the paper uses 100 x 100 per band; the hybrid costs ~0.06 s per corrected point):
  31 x 31 per window; 21 x 21 if the scenario has more than 6 windows (named dispersive materials);
  9 x 9 and basis order 8 for the Nx = 256 scenarios (thesis26/27); thesis28a/c (Nx = 1024, rough and
  Lipschitz profiles with 120 Fourier modes): AWE + indicator only.

    python run_scenarios.py                    # all 43 scenarios, 2 processes (~2-3 h)
    python run_scenarios.py --only dielectric gold
"""
import argparse
import json
import os
import sys
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from awepinn import core                 # noqa: E402
from awepinn.core import rm              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, 'figures')
RES = os.path.join(HERE, 'results', 'scenarios')
INDICATOR_ONLY = ('thesis28a', 'thesis28c')
# extra scenarios (not in refl_map.py): the paper's metal maps on a resolved x-grid
EXTRA = {'silver_nx64': ('silver', dict(Nx=64)), 'gold_nx64': ('gold', dict(Nx=64))}


def resolve(name):
    return EXTRA.get(name, (name, {}))
COARSE = ('thesis26', 'thesis27')


def grid_for(name):
    name = resolve(name)[0]
    sc = rm.SCENARIOS[name]
    if name in COARSE:
        return 9
    if sc.get('max_delta') or sc.get('windows') == 'joint':
        return 21
    return 31


def stats(res, info):
    lossless = bool(info['lossless'])
    cat = lambda k: np.concatenate([np.real(r[k]).ravel() for r in res])
    Ra, Rh, Da, Dh, ind = cat('ru_awe'), cat('ru'), cat('ee_awe'), cat('ee'), cat('indicator')
    sol = np.concatenate([r['solved'].ravel() for r in res])
    dR = np.abs(Ra - Rh)
    s = dict(windows=len(res), points=int(Ra.size), solved_frac=float(sol.mean()), lossless=lossless,
             t_awe=info['t_awe'], t_hybrid=info['t_hybrid'],
             dR_max=float(np.nanmax(dR)), dR_median=float(np.nanmedian(dR)), dR_p99=float(np.nanpercentile(dR, 99)),
             frac_awe_unreliable=float(np.mean(ind > 1e-10)), indicator_max=float(np.nanmax(ind)),
             awe_nonphysical=float(np.mean((Ra < -1e-6) | (Ra > 1 + 1e-6))),
             hyb_nonphysical=float(np.mean((Rh < -1e-6) | (Rh > 1 + 1e-6))))
    if lossless:
        s.update(D_awe_median_log10=float(np.log10(np.nanmedian(np.abs(Da)) + 1e-300)),
                 D_awe_max_log10=float(np.log10(np.nanmax(np.abs(Da)) + 1e-300)),
                 D_hyb_median_log10=float(np.log10(np.nanmedian(np.abs(Dh)) + 1e-300)),
                 D_hyb_max_log10=float(np.log10(np.nanmax(np.abs(Dh)) + 1e-300)))
    else:
        s.update(A_median=float(np.nanmedian(Dh)), dA_max=float(np.nanmax(np.abs(Da - Dh))))
    return s


def compare_figure(name, res, info, path):
    period = info.get('period')
    sc = period / (2 * np.pi) if period else 1.0
    unit = r' ($\mu$m)' if period else ''
    fig, axs = plt.subplots(2, 3, figsize=(19, 8.5))
    lossless = info['lossless']

    lam_all = np.concatenate([r['lam'] for r in res]) * sc
    logx = lam_all.max() / lam_all.min() > 5          # q = 0 bands reach the static limit: log axis

    def panel(ax, key_fn, title, cmap, vmin=None, vmax=None, clip=False):
        Zs = [key_fn(r) for r in res]
        z = np.concatenate([q[np.isfinite(q)].ravel() for q in Zs])
        lo, hi = (np.percentile(z, 0.5), np.percentile(z, 99.5)) if clip else (np.min(z), np.max(z))
        vmin = lo if vmin is None else vmin
        vmax = hi if vmax is None else vmax
        if logx:
            ax.set_xscale('log')
        for r, Z in zip(res, Zs):
            m = ax.pcolormesh(r['lam'] * sc, r['Eps'], Z, cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
        fig.colorbar(m, ax=ax)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(r'$\lambda$' + unit)
        ax.set_ylabel(r'$\varepsilon$')

    L = lambda v: np.log10(np.abs(v) + 1e-17)
    panel(axs[0, 0], lambda r: np.real(r['RR_awe']), 'HOPS/AWE: ' + ('R/R_flat' if info['relative'] else 'R') +
          ' (colours clipped to the 0.5-99.5 % range)', 'hot', clip=True)
    panel(axs[0, 1], lambda r: np.real(r['RR']), 'HOPS/AWE + PINN: ' + ('R/R_flat' if info['relative'] else 'R'), 'hot',
          clip=True)
    panel(axs[0, 2], lambda r: L(np.real(r['ru']) - np.real(r['ru_awe'])), r'$\log_{10}|R_{hybrid} - R_{AWE}|$',
          'viridis', -17, None)
    if lossless:
        panel(axs[1, 0], lambda r: L(r['ee_awe']), r'HOPS/AWE: $\log_{10}|D|$', 'hot', -17, 0)
        panel(axs[1, 1], lambda r: L(r['ee']), r'HOPS/AWE + PINN: $\log_{10}|D|$', 'hot', -17, 0)
    else:
        panel(axs[1, 0], lambda r: np.real(r['ee_awe']), 'HOPS/AWE: absorptance D (clipped 0.5-99.5 %)', 'magma', clip=True)
        panel(axs[1, 1], lambda r: L(np.real(r['ee']) - np.real(r['ee_awe'])), r'$\log_{10}|D_{hybrid} - D_{AWE}|$',
              'viridis', -17, None)
    panel(axs[1, 2], lambda r: L(r['indicator']), r'indicator: $\log_{10}$ PINN loss of the AWE sum', 'magma')
    for r in res:
        ii, jj = np.nonzero(r['solved'])
        axs[1, 2].plot(r['lam'][jj] * sc, r['Eps'][ii], 'c.', ms=1.2)
    s = info['stats']
    fig.suptitle(f"{name}: {info.get('desc', '')[:110]}\nmax |R_hyb - R_AWE| = {s['dR_max']:.1e},  "
                 f"{100 * s['solved_frac']:.0f} % of points corrected (cyan),  AWE {s['t_awe']:.1f} s, hybrid "
                 f"{s['t_hybrid']:.0f} s" + (f",  max log10|D|: AWE {s['D_awe_max_log10']:.1f}, hybrid "
                                             f"{s['D_hyb_max_log10']:.1f}" if lossless else ''), fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


def run_one(name, force=False):
    warnings.filterwarnings('ignore')
    out_json = os.path.join(RES, f'{name}.json')
    if os.path.exists(out_json) and not force:
        return json.load(open(out_json))
    t0 = time.time()
    n = grid_for(name)
    tol = np.inf if name in INDICATOR_ONLY else core.DEFAULT_TOL
    try:
        base, over = resolve(name)
        res, info = core.run(base, N_Eps=n, N_delta=n, tol=tol, workers=1, verbose=False,
                             basis_order=8 if name in COARSE else core.DEFAULT_ORDER, **over)
    except Exception as e:                                          # noqa: BLE001
        traceback.print_exc()
        return dict(name=name, error=repr(e))
    s = stats(res, info)
    s.update(name=name, grid=n, desc=info.get('desc', ''), n_w=str(info['n_w_spec']), n_u=str(info['n_u_spec']),
             mode=info['mode'], summation='Taylor' if info['Taylor'] else 'Pade', N=info['N'], M=info['M'],
             Nx=info['Nx'], indicator_only=name in INDICATOR_ONLY, wall=time.time() - t0)
    info['stats'] = s
    d = os.path.join(FIG, name)
    os.makedirs(d, exist_ok=True)
    awe = [dict(r, ru=r['ru_awe'], rl=r['rl_awe'], ee=r['ee_awe'], RR=r['RR_awe']) for r in res]
    rm.plot(awe, info, outdir=d, tag=f'{name}_awe')
    rm.plot(res, info, outdir=d, tag=f'{name}_hybrid')
    plt.close('all')
    compare_figure(name, res, info, os.path.join(d, 'compare.png'))
    os.makedirs(RES, exist_ok=True)
    arr = {}
    for k, r in enumerate(res):
        for key in ('lam', 'omega', 'Eps', 'ru', 'rl', 'ee', 'ru_awe', 'rl_awe', 'ee_awe', 'ru_flat', 'indicator', 'solved'):
            arr[f'w{k}_{key}'] = np.asarray(r[key])
        arr[f'w{k}_n_w'] = np.asarray(complex(r['n_w']))
        arr[f'w{k}_omega_bar'] = np.asarray(r['omega_bar'])
    np.savez_compressed(os.path.join(RES, f'{name}.npz'), **arr)
    json.dump(s, open(out_json, 'w'), indent=1)
    print(f"{name:22s} grid {n:2d} windows {s['windows']:2d}  max|dR| {s['dR_max']:.1e}  corrected "
          f"{100 * s['solved_frac']:4.0f} %  AWE {s['t_awe']:6.1f} s  hybrid {s['t_hybrid']:6.0f} s" +
          (f"  log10|D| max AWE {s['D_awe_max_log10']:6.1f} -> hybrid {s['D_hyb_max_log10']:6.1f}" if s['lossless'] else ''),
          flush=True)
    return s


def table(rows):
    rows = [r for r in rows if 'error' not in r]
    L = ['| scenario | n_w | mode, summation | windows x grid | corrected points | max \\|R_hyb - R_AWE\\| | '
         'max log10\\|D\\| AWE -> hybrid (lossless) | AWE / hybrid time |', '|---|---|---|---|---|---|---|---|']
    for s in rows:
        dd = (f"{s['D_awe_max_log10']:.1f} -> {s['D_hyb_max_log10']:.1f}" if s['lossless'] else
              f"absorbing (A ~ {s['A_median']:.2f})")
        L.append(f"| {s['name']} | {s['n_w']} | {s['mode']}, {s['summation']} {s['N']}/{s['M']} | {s['windows']} x "
                 f"{s['grid']}^2 | {100 * s['solved_frac']:.0f} %{' (indicator only)' if s['indicator_only'] else ''} | "
                 f"{s['dR_max']:.1e} | {dd} | {s['t_awe']:.1f} / {s['t_hybrid']:.0f} s |")
    open(os.path.join(HERE, 'results', 'scenarios.md'), 'w').write('\n'.join(L) + '\n')
    return L


def replot(name):
    """rebuild compare.png from the saved arrays"""
    s = json.load(open(os.path.join(RES, f'{name}.json')))
    if 'error' in s:
        return
    d = np.load(os.path.join(RES, f'{name}.npz'))
    nwin = len({k.split('_')[0] for k in d.files})
    sc = rm.SCENARIOS[resolve(name)[0]]
    relative = sc.get('relative', rm.PlotRelative)
    res = []
    for k in range(nwin):
        r = {q: d[f'w{k}_{q}'] for q in ('lam', 'Eps', 'ru', 'rl', 'ee', 'ru_awe', 'ee_awe', 'ru_flat', 'indicator', 'solved')}
        r['RR'] = r['ru'] / r['ru_flat'] if relative else r['ru']
        r['RR_awe'] = r['ru_awe'] / r['ru_flat'] if relative else r['ru_awe']
        res.append(r)
    period = sc.get('period')
    if period is None and rm.mat.parse_index(sc['n_w']) is None:
        period = rm.mat.MATERIALS.get(sc['n_w'], {}).get('period')
    info = dict(lossless=s['lossless'], relative=relative, period=period, desc=s['desc'], stats=s)
    compare_figure(name, res, info, os.path.join(FIG, name, 'compare.png'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--replot', action='store_true', help='only rebuild the compare figures from saved results')
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--procs', type=int, default=2)
    ap.add_argument('--force', action='store_true', help='recompute even if results/scenarios/<name>.json exists')
    a = ap.parse_args()
    names = a.only or list(rm.SCENARIOS)
    if a.replot:
        for k in names:
            if os.path.exists(os.path.join(RES, f'{k}.json')):
                replot(k)
        sys.exit(0)
    # cheap first, heavy last
    names = sorted(names, key=lambda k: (k in INDICATOR_ONLY, k in COARSE, len(rm.SCENARIOS[resolve(k)[0]]['q']) * grid_for(k) ** 2))
    if a.procs > 1:
        with ProcessPoolExecutor(max_workers=a.procs) as ex:
            rows = list(ex.map(run_one, names, [a.force] * len(names)))
    else:
        rows = [run_one(k, a.force) for k in names]
    rows = [json.load(open(os.path.join(RES, f'{k}.json'))) for k in list(rm.SCENARIOS) + list(EXTRA)
            if os.path.exists(os.path.join(RES, f'{k}.json'))]
    print('\n'.join(table(rows)))
