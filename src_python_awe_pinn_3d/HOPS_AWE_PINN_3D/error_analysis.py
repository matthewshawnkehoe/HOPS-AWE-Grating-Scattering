"""error_analysis.py -- rigorous comparison of 3D HOPS/AWE and 3D HOPS/AWE + PINN against converged references.

Samples: in every frequency window (every k-th window for scenarios with many windows) two frequencies,
delta = -0.8 delta_max and +0.45 delta_max, times four heights eps = (0.15, 0.4, 0.65, 1.0) eps_max.
At each sample:
  R_true, D_true   HOPS solved POINTWISE in frequency (delta = 0: no frequency expansion), N = 24, Pade in eps,
                   on the refined grid (2 Nx, 2 Ny, Nz + 16)                        -> total error
  R_same, D_same   the same on the scenario's own grid (Nx, Ny, Nz)                  -> summation error
  R_awe            refl_map_3D.py's summation (MATLAB-style Taylor or Pade), exactly as in the maps
  R_taylor         the full-order Taylor double sum (what the indicator measures)
  R_hyb{6,8,10,12} PI-sum least-squares weights with basis orders n, m <= 6 ... 12 (discrete form)
  R_chain          PI-sum, order 12, with the 2D hybrid's chain-rule rows (comparison of the two forms)
  indicator        PINN loss of the full-order AWE Taylor sum
The adaptive hybrid (the method of the maps) = R_taylor where indicator <= tol, R_hyb12 elsewhere.

    python error_analysis.py                    # all scenarios (~1-2 h, 2 processes)
    python error_analysis.py --only gold        # subset;  --report-only rebuilds the report from saved json
"""
import argparse
import json
import os
import sys
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor
from functools import partial

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from awepinn3d import core, reference as ref     # noqa: E402
from awepinn3d.core import rm3                   # noqa: E402
import hops3d as h3                              # noqa: E402
from run_scenarios import EXTRA, resolve, n_windows   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results', 'error_analysis')
FIGD = os.path.join(HERE, 'figures', 'error_analysis')
DFRAC = (-0.8, 0.45)
EFRAC = (0.15, 0.4, 0.65, 1.0)
ORDERS = (6, 8, 10, 12)
TOL = core.DEFAULT_TOL


def analyse_window(win, c, opts):
    """refl_map_3D window -> sample rows (runs in place of refl_map_3D._window)"""
    key, omega_bar, dmax, edges = win
    base = dict(key=key, omega_bar=omega_bar, delta=np.zeros(1), omega=np.array([omega_bar]),
                lam=np.array([2 * np.pi / omega_bar]), Eps=c['Eps'], n_u=c['n_u_w'][key], n_w=c['n_w_w'][key])
    k = int(key[1:]) if key[1:].isdigit() else 0
    if k % opts['every']:
        return dict(base, rows=[])
    warnings.filterwarnings('ignore')
    t0 = time.time()
    N, M = c['N'], c['M']
    P, n_u, n_w, alpha, beta = core.window_problem(win, c)
    zeta, psi = h3.setup_zeta_psi_n_m_3d(alpha, beta, P.gamma_u_bar, P.f, P.f_x, P.f_y, N, M)
    U, W, ubar, wbar, vu, vw = h3.two_layer_solve_3d_coupled(P, zeta, psi, keep_volume=True)
    eps_s = np.array(EFRAC) * opts['eps_max']
    dl_s = np.array(DFRAC) * dmax
    args = (P.tau2, ubar, wbar, P.kx, P.ky, alpha, beta, P.gamma_u_bar, P.gamma_w_bar, eps_s, dl_s)
    ee, ru, rl = h3.energy_defect_3d(*args, N, M, 1 if c['Taylor'] else 2,
                                     taylor_full_order=c.get('taylor_full_order', False), sum_domain=c['sum_domain'])
    rows = {(i, j): dict(scenario=opts['name'], window=key, eps=float(e), delta=float(d), omega=float(omega_bar * (1 + d)),
                         R_awe=float(np.real(ru[i, j])), D_awe=float(np.real(ee[i, j])))
            for i, e in enumerate(eps_s) for j, d in enumerate(dl_s)}
    # truth: pointwise HOPS at each sample frequency (refined grid, and the scenario grid)
    info = dict(profile=opts['profile'], mode='TE' if c['Mode'] == 1 else 'TM', a=c['a'], b=c['b'],
                Nx=c['Nx'], Ny=c['Ny'], Nz=c['Nz'])
    for j, d in enumerate(dl_s):
        s = 1 + d
        for tag, fine in (('true', True), ('same', False)):
            Nx, Ny, Nz = ref.refined(info, fine)
            R, T, D = ref.hops_pointwise(info, n_u, n_w, omega_bar * s, alpha * s, beta * s, eps_s, Nx=Nx, Ny=Ny, Nz=Nz)
            for i in range(len(eps_s)):
                rows[i, j].update({f'R_{tag}': float(R[i]), f'D_{tag}': float(D[i])})
    # hybrids
    for form, orders in (('discrete', ORDERS), ('chain', (12,))):
        for order in orders:
            Wn = core.Window3D(vu, vw, P, n_u, n_w, omega_bar, alpha, beta, N, M, basis_order=order, form=form)
            tag = f'hyb{order}' if form == 'discrete' else 'chain'
            for (i, j), row in rows.items():
                if form == 'discrete' and order == 12:
                    sI = Wn.indicator(row['eps'], row['omega'])
                    row.update(indicator=sI['loss_awe_taylor'], R_taylor=sI['R_taylor'], D_taylor=sI['D_taylor'])
                if form == 'chain':
                    row['indicator_chain'] = Wn.indicator(row['eps'], row['omega'])['loss_awe_taylor']
                sol = Wn.solve(row['eps'], row['omega'])
                row.update({f'R_{tag}': sol['R'], f'D_{tag}': sol['D'], f'loss_{tag}': sol['loss']})
            if order == 12 and form == 'discrete':
                nunk = Wn.n_unknowns
            del Wn
    for row in rows.values():
        row.update(n_unknowns=nunk, lossless=bool(abs(np.imag(n_w)) < 1e-6 and abs(np.imag(n_u)) < 1e-6))
    if c['verbose']:
        print(f"  {opts['name']} {key}: {len(rows)} samples, {time.time() - t0:.0f} s", flush=True)
    return dict(base, rows=list(rows.values()))


def analyse(name):
    warnings.filterwarnings('ignore')
    out = os.path.join(OUT, f'{name}.json')
    if os.path.exists(out):
        return json.load(open(out))
    t0 = time.time()
    base, over = resolve(name)
    sc = dict(rm3.SCENARIOS[base]); sc.update(over)
    nw = n_windows(name)
    opts = dict(name=name, profile=sc['profile'], eps_max=sc['eps_max'], every=max(1, nw // 12))
    orig = rm3._window
    rm3._window = partial(analyse_window, opts=opts)
    try:
        res, info = rm3.run(base, N_Eps=2, N_delta=2, verbose=True, workers=1, **over)
    except Exception as e:                                          # noqa: BLE001
        traceback.print_exc()
        return dict(name=name, error=repr(e))
    finally:
        rm3._window = orig
    rows = [row for r in res for row in r['rows']]
    os.makedirs(OUT, exist_ok=True)
    d = dict(name=name, rows=rows, wall=time.time() - t0, summation='Taylor' if info['Taylor'] else 'Pade',
             n_w=str(info['n_w_spec']), Nx=info['Nx'], Ny=info['Ny'], Nz=info['Nz'], N=info['N'], M=info['M'])
    json.dump(d, open(out, 'w'))
    print(f'{name}: {len(rows)} samples in {time.time() - t0:.0f} s', flush=True)
    return d


# ---------------------------------------------------------------------------------------------------
def adaptive(rows, tol=TOL, key='R'):
    return np.array([r[f'{key}_taylor'] if r['indicator'] <= tol else r[f'{key}_hyb12'] for r in rows])


def wilcoxon(a, b):
    from scipy.stats import wilcoxon as w
    try:
        return float(w(a, b).pvalue)
    except ValueError:
        return float('nan')


def report():
    from scipy.stats import spearmanr
    names = list(rm3.SCENARIOS) + list(EXTRA)
    data = {}
    for n in names:
        p = os.path.join(OUT, f'{n}.json')
        if os.path.exists(p):
            d = json.load(open(p))
            if 'error' not in d and d['rows']:
                data[n] = d
    if not data:
        print('no results yet')
        return
    A = lambda rows, k: np.array([r[k] for r in rows])
    L = ['Errors against the converged reference (pointwise HOPS at delta = 0, N = 24, Pade, 2 Nx, 2 Ny, Nz + 16) '
         '-- total error -- and against the same pointwise HOPS on the scenario\'s own grid (same Nx, Ny, Nz) -- the '
         'SUMMATION error, which is what the PINN step addresses.', '',
         '| scenario | n_w | grid | samples | total max err R: AWE (refl_map_3D) | AWE full Taylor | **hybrid** | '
         'summation max err R: AWE / **hybrid** | median total err R: AWE / hybrid | max err D: AWE / hybrid | '
         'hybrid better (> 2x) at | Wilcoxon p |', '|---|---|---|---|---|---|---|---|---|---|---|---|']
    pool = dict(ea=[], eh=[], sa=[], sh=[], ind=[], et=[], st=[])
    order_rows, chain_rows = [], []
    for n, d in data.items():
        rows = d['rows']
        Rt, Rs = A(rows, 'R_true'), A(rows, 'R_same')
        Ra, Rtay, Rh = A(rows, 'R_awe'), A(rows, 'R_taylor'), adaptive(rows)
        Dt, Da, Dh = A(rows, 'D_true'), A(rows, 'D_awe'), adaptive(rows, key='D')
        ea, eh, et = np.abs(Ra - Rt), np.abs(Rh - Rt), np.abs(Rtay - Rt)
        sa, sh, st = np.abs(Ra - Rs), np.abs(Rh - Rs), np.abs(Rtay - Rs)
        for k, v in (('ea', ea), ('eh', eh), ('sa', sa), ('sh', sh), ('et', et), ('st', st)):
            pool[k].append(v)
        pool['ind'].append(A(rows, 'indicator'))
        better = np.mean(eh < 0.5 * ea)
        L.append(f"| {n} | {d['n_w']} | {d['Nx']}x{d['Ny']}x{d['Nz']} | {len(rows)} | {ea.max():.1e} | {et.max():.1e} | "
                 f"**{eh.max():.1e}** | {sa.max():.1e} / **{sh.max():.1e}** | {np.median(ea):.1e} / {np.median(eh):.1e} | "
                 f"{np.abs(Da - Dt).max():.1e} / {np.abs(Dh - Dt).max():.1e} | {100 * better:.0f} % | {wilcoxon(ea, eh):.0e} |")
        order_rows.append((n, [np.max(np.abs(A(rows, f'R_hyb{o}') - Rs)) for o in ORDERS]))
        chain_rows.append((n, np.max(np.abs(A(rows, 'R_hyb12') - Rs)), np.max(np.abs(A(rows, 'R_chain') - Rs)),
                           np.median(A(rows, 'indicator')), np.median(A(rows, 'indicator_chain'))))
    P = {k: np.concatenate(v) for k, v in pool.items()}
    nS = P['ea'].size
    L += ['', f"Pooled over {nS} samples ({len(data)} scenarios): max error AWE {P['ea'].max():.1e}, hybrid {P['eh'].max():.1e}; "
          f"median AWE {np.median(P['ea']):.1e}, hybrid {np.median(P['eh']):.1e}; hybrid better by > 2x at "
          f"{100 * np.mean(P['eh'] < 0.5 * P['ea']):.0f} % of the samples, worse by > 2x at {100 * np.mean(P['eh'] > 2 * P['ea']):.0f} %.",
          f"Summation error alone (against HOPS on the same grid): max AWE {P['sa'].max():.1e}, hybrid {P['sh'].max():.1e}; "
          f"median AWE {np.median(P['sa']):.1e}, hybrid {np.median(P['sh']):.1e}; hybrid better by > 2x at "
          f"{100 * np.mean(P['sh'] < 0.5 * P['sa']):.0f} %, worse by > 2x at {100 * np.mean(P['sh'] > 2 * P['sa']):.0f} %.",
          f"Excluding round-off (samples where the larger of the two errors is below 1e-12): hybrid better by > 2x at "
          f"{100 * np.mean((P['eh'] < 0.5 * P['ea']) & (np.maximum(P['ea'], P['eh']) > 1e-12)):.0f} %, worse by > 2x at "
          f"{100 * np.mean((P['eh'] > 2 * P['ea']) & (np.maximum(P['ea'], P['eh']) > 1e-12)):.0f} % of all samples "
          f"(total error); errors below 1e-12 for both at {100 * np.mean(np.maximum(P['ea'], P['eh']) <= 1e-12):.0f} %."]
    m = (P['st'] > 1e-15) & (P['ind'] > 0)
    rho = spearmanr(np.log10(P['ind'][m]), np.log10(P['st'][m])).correlation
    fit = np.polyfit(np.log10(P['ind'][m]), np.log10(P['st'][m]), 1)
    L.append(f'Indicator calibration (full-order AWE Taylor sum vs its summation error): Spearman rho = {rho:.2f}; '
             f'log10|R error| ~ {fit[0]:.2f} log10(indicator) {fit[1]:+.2f}.')
    tols = (1e-20, 1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6)
    allrows = [r for d in data.values() for r in d['rows']]
    Rs_all = A(allrows, 'R_same')
    L += ['', '| indicator tolerance | ' + ' | '.join(f'{t:.0e}' for t in tols) + ' |', '|---' * (len(tols) + 1) + '|',
          '| fraction of samples solved | ' + ' | '.join(f"{100 * np.mean(A(allrows, 'indicator') > t):.0f} %" for t in tols) + ' |',
          '| max summation error | ' + ' | '.join(f"{np.max(np.abs(adaptive(allrows, t) - Rs_all)):.1e}" for t in tols) + ' |',
          '| median summation error | ' + ' | '.join(f"{np.median(np.abs(adaptive(allrows, t) - Rs_all)):.1e}" for t in tols) + ' |']
    L += ['', 'Basis order (max summation error of the least-squares weights, every sample solved):', '',
          '| scenario | ' + ' | '.join(f'order {o}' for o in ORDERS) + ' |', '|---' * (len(ORDERS) + 1) + '|']
    L += [f'| {n} | ' + ' | '.join(f'{v:.1e}' for v in vals) + ' |' for n, vals in order_rows]
    L += ['', 'Residual form (order 12, every sample solved): hops3d\'s discrete flux form (default) vs the chain-rule '
          'form of the 2D hybrid:', '', '| scenario | max summation err R: discrete | chain rule | median indicator: discrete | chain rule |',
          '|---|---|---|---|---|']
    L += [f'| {n} | {a:.1e} | {b:.1e} | {c:.1e} | {e:.1e} |' for n, a, b, c, e in chain_rows]
    open(os.path.join(HERE, 'results', 'error_analysis.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    # ---- figures
    os.makedirs(FIGD, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(20, 6))
    ax[0].loglog(P['ind'][m], P['st'][m], '.', ms=3, alpha=0.5, label='full-order AWE Taylor sum')
    xx = np.logspace(np.log10(P['ind'][m].min()), np.log10(P['ind'][m].max()), 10)
    ax[0].loglog(xx, 10 ** np.polyval(fit, np.log10(xx)), 'k--', lw=1, label=f'fit: slope {fit[0]:.2f}')
    ax[0].axvline(TOL, color='r', ls=':', label=f'tol = {TOL:.0e}')
    ax[0].set_xlabel('indicator (PINN loss of the AWE sum; no reference needed)')
    ax[0].set_ylabel('|R error| (summation)')
    ax[0].set_title(f'indicator calibration, {m.sum()} samples, Spearman {rho:.2f}')
    ax[0].legend(fontsize=8)
    names_ = list(data)
    y = np.arange(len(names_))
    mx = lambda k: [np.max(v) for v in pool[k]]
    ax[1].barh(y + 0.2, mx('ea'), 0.4, label='3D HOPS/AWE (refl_map_3D)', color='C0')
    ax[1].barh(y - 0.2, mx('eh'), 0.4, label='3D HOPS/AWE + PINN', color='C3')
    ax[1].set_xscale('log'); ax[1].set_yticks(y); ax[1].set_yticklabels(names_, fontsize=7); ax[1].invert_yaxis()
    ax[1].set_xlabel('max |R - R_true| (total error, refined-grid reference)'); ax[1].legend(fontsize=8, loc='lower left')
    ax[1].set_title('total error per scenario')
    ax[2].barh(y + 0.2, mx('sa'), 0.4, label='3D HOPS/AWE', color='C0')
    ax[2].barh(y - 0.2, mx('sh'), 0.4, label='3D HOPS/AWE + PINN', color='C3')
    ax[2].set_xscale('log'); ax[2].set_yticks(y); ax[2].set_yticklabels(names_, fontsize=7); ax[2].invert_yaxis()
    ax[2].set_xlabel('max |R - R_same| (summation error, same-grid reference)'); ax[2].legend(fontsize=8, loc='lower left')
    ax[2].set_title('summation error per scenario')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGD, 'error_analysis.png'), dpi=90)
    plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--procs', type=int, default=2)
    ap.add_argument('--report-only', action='store_true')
    a = ap.parse_args()
    if not a.report_only:
        names = a.only or (list(rm3.SCENARIOS) + list(EXTRA))
        big = lambda k: rm3.SCENARIOS[resolve(k)[0]]['Nx'] >= 32 and rm3.SCENARIOS[resolve(k)[0]].get('Ny', 16) >= 32
        names = sorted(names, key=lambda k: (big(k), n_windows(k)))      # cheap first
        if a.procs > 1:
            with ProcessPoolExecutor(max_workers=a.procs) as ex:
                list(ex.map(analyse, names))
        else:
            for n in names:
                analyse(n)
    report()
