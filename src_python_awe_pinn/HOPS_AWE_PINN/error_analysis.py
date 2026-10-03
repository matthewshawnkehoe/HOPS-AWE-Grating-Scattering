"""error_analysis.py -- rigorous error analysis: standard HOPS/AWE vs HOPS/AWE + PINN against an independent
reference, for every scenario whose profile the reference solver supports.

Reference ("truth"): HOPS solved separately at each sampled (eps, omega) -- delta = 0, i.e. NO frequency
expansion -- with N = 24 orders, Pade summation, twice the scenario's Nx (>= 64) and Nz + 16, and the
window's own refractive indices (so piecewise-constant dispersion is the same as in AWE).

Samples: per window (at most 4 windows per scenario, spread over the bands) 5 eps values x 7 delta values
(including both band edges) of the run_scenarios.py grid -> up to 140 points per scenario.
At every sample:  R, D of refl_map.py (AWE as published: Taylor half order / Pade), of the AWE full-order
Taylor sum, of the hybrid with the least-squares weights (basis orders 6, 8, 10, 12) and of the adaptive
hybrid (tolerance sweep), plus the indicator.  Outputs:
  results/error_analysis/<scenario>.json    all samples
  results/error_analysis.md                  per-scenario table, paired statistics (Wilcoxon), calibration
  figures/error_analysis/*.png               error vs indicator, vs band position, vs basis order, tolerance curve
"""
import argparse
import json
import os
import sys
import time
import warnings

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats as sps

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from awepinn import core                                       # noqa: E402
from awepinn.core import rm, WindowBand, HybridSum              # noqa: E402
from run_scenarios import grid_for, resolve, EXTRA              # noqa: E402

PINN = os.path.join(os.path.dirname(core.HERE), 'PINN_HOPS')
if PINN not in sys.path:
    sys.path.insert(0, PINN)
from pinn_hops.hops_reference import hops_point                 # noqa: E402
from pinn_hops.problem import Grating2D, PROFILES               # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results', 'error_analysis')
FIG = os.path.join(HERE, 'figures', 'error_analysis')
PROFILE_MAP = {'fs1': 'cos4x_over4'}
ORDERS = (6, 8, 10, 12)
TOLS = (1e-20, 1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6)


def supported(name):
    sc = rm.SCENARIOS[resolve(name)[0]]
    p = PROFILE_MAP.get(sc['profile'], sc['profile'])
    return p in PROFILES and sc['eps_max'] <= 0.4 and sc['Nx'] <= 64


def analyse(name, verbose=True):
    fn = os.path.join(OUT, f'{name}.json')
    if os.path.exists(fn):
        return json.load(open(fn))
    warnings.filterwarnings('ignore')
    t0 = time.time()
    n = grid_for(name)
    base, over = resolve(name)
    res, info = rm.run(base, N_Eps=n, N_delta=n, verbose=False, workers=1, keep_fields=True, **over)
    res = sorted(res, key=lambda r: r['omega_bar'])
    if len(res) > 4:
        res = [res[i] for i in np.linspace(0, len(res) - 1, 4).round().astype(int)]
    prof = PROFILE_MAP.get(info['profile'], info['profile'])
    rows = []
    ie = np.linspace(1, n - 1, 5).round().astype(int)
    jd = np.linspace(0, n - 1, 7).round().astype(int)
    for w, r in enumerate(res):
        Bs = {o: WindowBand(r, info, basis_order=o) for o in ORDERS}
        Ss = {o: HybridSum(B, compress=core.COMPRESS) for o, B in Bs.items()}
        dmax = np.max(np.abs(r['delta']))
        for i in ie:
            for j in jd:
                e, dl = float(r['Eps'][i]), float(r['delta'][j])
                om = float(r['omega_bar'] * (1 + dl))
                G = Grating2D(n_u=complex(r['n_u']), n_w=complex(r['n_w']), omega=om, eps=e, profile=prof,
                              alpha=info['alpha'] * (1 + dl), a=info['a'], b=info['b'], mode=info['mode'])
                ref = hops_point(G, N=24, Nx=max(64, 2 * info['Nx']), Nz=info['Nz'] + 16, summation='pade',
                                 fields=False)
                same = hops_point(G, N=24, Nx=info['Nx'], Nz=info['Nz'], summation='pade', fields=False)
                row = dict(window=w, eps=e, delta=dl, delta_rel=dl / dmax if dmax else 0.0, omega=om,
                           R_true=float(ref['R']), D_true=float(ref['D']),
                           R_same=float(same['R']), D_same=float(same['D']),
                           R_awe=float(np.real(r['ru'][i, j])), D_awe=float(np.real(r['ee'][i, j])))
                for o, S in Ss.items():
                    s = S.point(e, om, tol=-1.0)                    # always least squares
                    row[f'R_hyb{o}'], row[f'D_hyb{o}'], row[f'loss_hyb{o}'] = float(s['R']), float(s['D']), s['loss']
                row.update(R_taylor=float(s['R_taylor']), D_taylor=float(s['D_taylor']),
                           indicator=float(s['loss_awe_taylor']))
                rows.append(row)
    out = dict(name=name, desc=info.get('desc', ''), lossless=bool(info['lossless']), n_w=str(info['n_w_spec']),
               mode=info['mode'], summation='Taylor' if info['Taylor'] else 'Pade', N=info['N'], M=info['M'],
               Nx=info['Nx'], samples=rows, time=time.time() - t0)
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(fn, 'w'))
    if verbose:
        e = lambda k: np.array([abs(x[k] - x['R_true']) for x in rows])
        print(f"{name:22s} {len(rows):4d} samples  max|R err|: AWE {e('R_awe').max():.1e}  hybrid(12) "
              f"{e('R_hyb12').max():.1e}   ({out['time']:.0f} s)", flush=True)
    return out


def adaptive(rows, tol, o=12):
    """R of the adaptive hybrid: Taylor where indicator <= tol, least squares elsewhere"""
    return np.array([x['R_taylor'] if x['indicator'] <= tol else x[f'R_hyb{o}'] for x in rows])


def report(all_out):
    os.makedirs(FIG, exist_ok=True)
    L = ['Errors against the converged reference (2 x Nx, Nz + 16, N = 24) -- total error -- and against HOPS on the',
         "scenario's own grid at delta = 0 (same Nx, Nz) -- the SUMMATION error, which is what the PINN step addresses.", '',
         '| scenario | n_w | samples | total max err R: AWE (refl_map) | AWE full Taylor | **hybrid** | summation max err R: AWE / **hybrid** | '
         'median total err R: AWE / hybrid | max err D: AWE / hybrid | hybrid better (> 2x) at | Wilcoxon p |',
         '|---|---|---|---|---|---|---|---|---|---|---|']
    pooled = dict(awe=[], hyb=[], ind=[], tay=[], drel=[])
    for o in all_out:
        rows = o['samples']
        Rt = np.array([x['R_true'] for x in rows])
        Dt = np.array([x['D_true'] for x in rows])
        ea = np.abs(np.array([x['R_awe'] for x in rows]) - Rt)
        et = np.abs(np.array([x['R_taylor'] for x in rows]) - Rt)
        eh = np.abs(adaptive(rows, core.DEFAULT_TOL) - Rt)
        da = np.abs(np.array([x['D_awe'] for x in rows]) - Dt)
        dh = np.abs(np.array([x['D_hyb12'] if x['indicator'] > core.DEFAULT_TOL else x['D_taylor'] for x in rows]) - Dt)
        Rs = np.array([x['R_same'] for x in rows])
        sa = np.abs(np.array([x['R_awe'] for x in rows]) - Rs)
        sh = np.abs(adaptive(rows, core.DEFAULT_TOL) - Rs)
        la, lh = np.log10(ea + 1e-17), np.log10(eh + 1e-17)
        better = np.mean(lh < la - 0.3)          # better by more than a factor 2
        try:
            p = sps.wilcoxon(la, lh, alternative='greater').pvalue
        except ValueError:
            p = np.nan
        L.append(f"| {o['name']} | {o['n_w']} | {len(rows)} | {ea.max():.1e} | {et.max():.1e} | **{eh.max():.1e}** | "
                 f"{sa.max():.1e} / **{sh.max():.1e}** | {np.median(ea):.1e} / {np.median(eh):.1e} | {da.max():.1e} / {dh.max():.1e} | "
                 f"{100 * better:.0f} % | {p:.0e} |")
        pooled.setdefault('sawe', []).extend(list(sa))
        pooled.setdefault('shyb', []).extend(list(sh))
        o['_e'] = dict(awe=ea, hyb=eh, tay=et)
        pooled['awe'] += list(ea)
        pooled['hyb'] += list(eh)
        pooled['tay'] += list(et)
        pooled['ind'] += [x['indicator'] for x in rows]
        pooled['drel'] += [abs(x['delta_rel']) for x in rows]
    P = {k: np.array(v) for k, v in pooled.items()}
    rho = sps.spearmanr(P['ind'], P['tay'])[0]
    m = (P['ind'] > 1e-24) & (P['tay'] > 1e-15)
    slope, icpt = np.polyfit(np.log10(P['ind'][m]), np.log10(P['tay'][m]), 1)
    L += ['', f"Pooled over {len(P['awe'])} samples: max error AWE {P['awe'].max():.1e}, hybrid {P['hyb'].max():.1e}; "
              f"median AWE {np.median(P['awe']):.1e}, hybrid {np.median(P['hyb']):.1e}; hybrid better by > 2x at "
              f"{100 * np.mean(np.log10(P['hyb'] + 1e-17) < np.log10(P['awe'] + 1e-17) - 0.3):.0f} % of the samples, "
              f"worse by > 2x at {100 * np.mean(np.log10(P['hyb'] + 1e-17) > np.log10(P['awe'] + 1e-17) + 0.3):.0f} %.",
          f"Summation error alone (against HOPS on the same grid): max AWE {P['sawe'].max():.1e}, hybrid {P['shyb'].max():.1e}; "
          f"median AWE {np.median(P['sawe']):.1e}, hybrid {np.median(P['shyb']):.1e}.",
          f"Indicator calibration (AWE full-order Taylor sum): Spearman rho = {rho:.2f}; "
          f"log10|R error| ~ {slope:.2f} log10(indicator) {icpt:+.2f}."]
    # tolerance sweep
    L += ['', '| indicator tolerance | ' + ' | '.join(f'{t:.0e}' for t in TOLS) + ' |', '|---' * (len(TOLS) + 1) + '|']
    fr, mx = [], []
    for t in TOLS:
        errs, solved = [], []
        for o in all_out:
            rows = o['samples']
            Rt = np.array([x['R_true'] for x in rows])
            errs += list(np.abs(adaptive(rows, t) - Rt))
            solved += [x['indicator'] > t for x in rows]
        fr.append(np.mean(solved))
        mx.append(np.max(errs))
    L.append('| fraction of samples solved | ' + ' | '.join(f'{100 * f:.0f} %' for f in fr) + ' |')
    L.append('| max error over all samples | ' + ' | '.join(f'{v:.1e}' for v in mx) + ' |')
    # basis order
    L += ['', '| scenario | ' + ' | '.join(f'basis order {o}' for o in ORDERS) + ' |', '|---' * (len(ORDERS) + 1) + '|']
    for o in all_out:
        rows = o['samples']
        Rt = np.array([x['R_true'] for x in rows])
        L.append(f"| {o['name']} | " + ' | '.join(f"{np.abs(np.array([x[f'R_hyb{q}'] for x in rows]) - Rt).max():.1e}"
                                                  for q in ORDERS) + ' |')
    open(os.path.join(HERE, 'results', 'error_analysis.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    # figures
    fig, axs = plt.subplots(1, 3, figsize=(19, 5.3))
    ax = axs[0]
    ax.loglog(P['ind'] + 1e-30, P['tay'] + 1e-17, '.', ms=3, alpha=0.5, label='AWE full-order Taylor sum')
    xx = np.logspace(-25, 10, 50)
    ax.loglog(xx, 10 ** (slope * np.log10(xx) + icpt), 'k--', lw=1, label=f'fit: slope {slope:.2f}')
    ax.axvline(core.DEFAULT_TOL, color='r', ls=':', label=f'tol = {core.DEFAULT_TOL:.0e}')
    ax.set_xlabel('indicator (PINN loss of the AWE sum; no reference needed)')
    ax.set_ylabel('true |R error|')
    ax.set_title(f'indicator calibration, {len(P["ind"])} samples, Spearman {rho:.2f}')
    ax.set_ylim(1e-17, 10)
    ax.set_xlim(1e-26, 1e16)
    ax.legend(fontsize=8)
    ax = axs[1]
    bins = np.array([-0.01, 0.2, 0.5, 0.8, 1.01])
    for k, lab, c in (('awe', 'HOPS/AWE (refl_map)', 'C0'), ('hyb', 'HOPS/AWE + PINN', 'C3')):
        sel = [P[k][(P['drel'] > lo) & (P['drel'] <= hi)] for lo, hi in zip(bins[:-1], bins[1:])]
        cen = np.array([0.5 * (lo + hi) for (lo, hi), s_ in zip(zip(bins[:-1], bins[1:]), sel) if s_.size])
        med = [np.median(s_) for s_ in sel if s_.size]
        mx_ = [np.max(s_) for s_ in sel if s_.size]
        ax.semilogy(cen, med, 'o-', color=c, label=f'{lab}, median')
        ax.semilogy(cen, mx_, 's--', color=c, label=f'{lab}, max')
    ax.set_xlabel(r'$|\delta| / \delta_{max}$  (0 = band centre, 1 = band edge)')
    ax.set_ylabel('|R error|')
    ax.set_title('error vs position in the frequency window (all scenarios)')
    ax.legend(fontsize=8)
    ax = axs[2]
    names = [o['name'] for o in all_out]
    y = np.arange(len(names))
    ax.barh(y - 0.2, [o['_e']['awe'].max() for o in all_out], 0.4, log=True, label='HOPS/AWE (refl_map)', color='C0')
    ax.barh(y + 0.2, [max(o['_e']['hyb'].max(), 1e-16) for o in all_out], 0.4, log=True, label='HOPS/AWE + PINN',
            color='C3')
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xscale('log')
    ax.set_xlim(1e-16, 10)
    ax.invert_yaxis()
    ax.set_xlabel('max |R - R_true| over the samples')
    ax.legend(fontsize=8)
    ax.set_title('per scenario')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'error_analysis.png'), dpi=105)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.loglog(np.array(fr) * 100 + 1e-3, mx, 'o-')
    for t, f_, m_ in zip(TOLS, fr, mx):
        ax.annotate(f'{t:.0e}', (f_ * 100 + 1e-3, m_), fontsize=7)
    ax.set_xlabel('% of points needing the least-squares solve')
    ax.set_ylabel('max |R error| over all samples')
    ax.set_title('adaptive hybrid: tolerance trade-off')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'tolerance.png'), dpi=110)
    plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--report-only', action='store_true')
    a = ap.parse_args()
    names = a.only or [k for k in list(rm.SCENARIOS) + list(EXTRA) if supported(k)]
    outs = []
    for k in names:
        if a.report_only and not os.path.exists(os.path.join(OUT, f'{k}.json')):
            continue
        try:
            outs.append(analyse(k))
        except Exception as e:      # noqa: BLE001
            print('FAILED', k, repr(e), flush=True)
    report(outs)
