"""materials_hybrid_3D.py -- what do the 3D (crossed-grating) reflectivity maps and energy defects of the curated
material library (hops/materials.py, refractiveindex.info) teach us, once 3D HOPS/AWE + PINN is available as an
accurate reference and the PINN indicator as a reliability map?

Setup exactly as material_survey_3D.py: vacuum over the material, f = cos(x)cos(y), eps_max = 0.2, Pade,
N = M = 10, Nx = Ny = Nz = 16, 'joint' windows (also split at the substrate's Rayleigh frequencies when it is
lossless), omega in [1, 3] (bands q = 1, 2), the material's suggested period, the index evaluated per window;
11 x 11 points per window.  For every material:
  features   Re n, Im n, Re eps = Re n^2, lossless?, number of 3D Rayleigh anomalies of the substrate in the range,
             distance of Re eps to the surface-plasmon condition eps = -1
  AWE        fraction of the map where |R_AWE - R_hybrid| > 1e-6, non-physical R, indicator > 1e-10, max log10|D|
  hybrid     max |R_hyb - R_AWE|, max log10|D|, non-physical R
  truth      pointwise HOPS (delta = 0, N = 24, same grid) on the two edge columns of every window ->
             the TRUE summation errors of AWE and of the hybrid there
  physics    hybrid min R/R_flat (resonance depth), absorptance range
Outputs results/materials.csv, results/materials.md (Spearman correlations, decision tree, ranking),
figures/materials/*.png.

    python materials_hybrid_3D.py                # all curated materials, 2 processes
    python materials_hybrid_3D.py --analyse-only --gallery
"""
import argparse
import csv
import json
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats as sps

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from awepinn3d import core, reference as ref     # noqa: E402
from awepinn3d.core import rm3                   # noqa: E402
from hops import materials as mat                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'results', 'materials')
FIG = os.path.join(HERE, 'figures', 'materials')
SETUP = dict(profile='cosxcosy', M=10, Nx=16, Ny=16, Nz=16, relative=True, windows='joint', taylor=False)
Q = (1, 2)
SUPERSTRATES = ('air', 'vacuum', 'water', 'ethanol', 'glass_superstrate')


def rayleigh_count(res):
    """number of 3D Rayleigh frequencies |(p, q)| / Re n_w of the substrate inside the omega range (lossless only)"""
    n = 0
    for r in res:
        nw = complex(r['n_w'])
        if abs(nw.imag) > 1e-3 * abs(nw) or nw.real <= 0:
            continue
        lo, hi = r['omega'].min() * nw.real, r['omega'].max() * nw.real
        p = np.arange(-12, 13)
        rad = np.sqrt(p[:, None] ** 2 + p[None, :] ** 2)
        n += int(np.sum((rad > lo) & (rad <= hi)))
    return n


def one(key, n=11):
    fn = os.path.join(OUT, f'{key}.json')
    if os.path.exists(fn):
        return json.load(open(fn))
    warnings.filterwarnings('ignore')
    t0 = time.time()
    try:
        res, info = core.run('dielectric', qq=Q, N_Eps=n, N_delta=n, workers=1, verbose=False, n_w=key,
                             period=mat.MATERIALS[key]['period'], **SETUP)
    except Exception as e:                    # noqa: BLE001
        return dict(key=key, error=repr(e))
    cat = lambda k: np.concatenate([np.real(r[k]).ravel() for r in res])
    Ra, Rh, Da, Dh, ind, Rf = cat('ru_awe'), cat('ru'), cat('ee_awe'), cat('ee'), cat('indicator'), cat('ru_flat')
    ns = np.array([complex(r['n_w']) for r in res])
    lossless = bool(info['lossless'])
    epsr = (ns ** 2).real
    # truth on the edge columns (j = 0, -1) of every window: pointwise HOPS, same grid
    tinfo = dict(profile=info['profile'], mode=info['mode'], a=info['a'], b=info['b'], Nx=info['Nx'], Ny=info['Ny'],
                 Nz=info['Nz'])
    ea, eh = [], []
    for r in res:
        for j in (0, -1):
            s_ = 1 + r['delta'][j]
            Rt, _, _ = ref.hops_pointwise(tinfo, r['n_u'], r['n_w'], float(r['omega'][j]), r['alpha'] * s_,
                                          r['beta'] * s_, r['Eps'])
            ea.append(np.abs(np.real(r['ru_awe'][:, j]) - Rt))
            eh.append(np.abs(np.real(r['ru'][:, j]) - Rt))
    ea, eh = np.concatenate(ea), np.concatenate(eh)
    s = dict(key=key, category=mat.MATERIALS[key].get('category', ''), period=mat.MATERIALS[key]['period'],
             windows=len(res), n_re_mean=float(ns.real.mean()), n_im_mean=float(ns.imag.mean()),
             n_im_max=float(ns.imag.max()), epsr_min=float(epsr.min()), epsr_max=float(epsr.max()),
             spp_dist=float(np.min(np.abs(epsr + 1))), lossless=lossless, rayleigh=rayleigh_count(res),
             unreliable=float(np.mean(ind > 1e-10)), solved=float(np.mean(np.concatenate([r['solved'].ravel() for r in res]))),
             awe_nonphys=float(np.mean((Ra < -1e-6) | (Ra > 1 + 1e-6))),
             hyb_nonphys=float(np.mean((Rh < -1e-6) | (Rh > 1 + 1e-6))),
             dR_max=float(np.max(np.abs(Rh - Ra))), dR_median=float(np.median(np.abs(Rh - Ra))),
             frac_dR6=float(np.mean(np.abs(Rh - Ra) > 1e-6)), frac_dR3=float(np.mean(np.abs(Rh - Ra) > 1e-3)),
             err_awe_max=float(ea.max()), err_hyb_max=float(eh.max()),
             err_awe_med=float(np.median(ea)), err_hyb_med=float(np.median(eh)),
             hyb_better=float(np.mean(eh < 0.5 * ea)), hyb_worse=float(np.mean(eh > 2 * ea)),
             Rrel_min=float(np.percentile(Rh / Rf, 1)), Rrel_spread=float(np.percentile(Rh / Rf, 95) - np.percentile(Rh / Rf, 5)),
             A_min=float(np.min(Dh)) if not lossless else np.nan, A_max=float(np.max(Dh)) if not lossless else np.nan,
             D_awe_max=float(np.log10(np.max(np.abs(Da)) + 1e-300)) if lossless else np.nan,
             D_hyb_max=float(np.log10(np.max(np.abs(Dh)) + 1e-300)) if lossless else np.nan,
             t_awe=info['t_awe'], t_hybrid=info['t_hybrid'], time=time.time() - t0)
    os.makedirs(OUT, exist_ok=True)
    json.dump(s, open(fn, 'w'))
    np.savez_compressed(os.path.join(OUT, f'{key}.npz'), **{f'w{k}_{q}': np.asarray(r[q]) for k, r in enumerate(res)
                                                              for q in ('lam', 'Eps', 'ru', 'ru_awe', 'ru_flat', 'ee', 'ee_awe', 'indicator')},
             period=info['period'] or np.nan)
    print(f"{key:12s} n~{s['n_re_mean']:.2f}+{s['n_im_mean']:.2f}i  win {s['windows']:2d}  max|dR| {s['dR_max']:.1e}  "
          f"true err AWE {s['err_awe_max']:.1e} -> hybrid {s['err_hyb_max']:.1e}  nonphys AWE {100 * s['awe_nonphys']:.1f} % "
          f"-> {100 * s['hyb_nonphys']:.1f} %  ({s['time']:.0f} s)", flush=True)
    return s


def analyse(rows):
    rows = [r for r in rows if 'error' not in r]
    os.makedirs(FIG, exist_ok=True)
    keys = list(rows[0].keys())
    with open(os.path.join(HERE, 'results', 'materials.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    A = lambda k: np.array([r[k] for r in rows], float)
    y = np.log10(A('err_awe_max') + 1e-17)
    feats = {'Re n': A('n_re_mean'), 'Im n': A('n_im_mean'), 'Re eps': A('epsr_min'),
             'substrate Rayleigh anomalies in range': A('rayleigh'), 'lossless': A('lossless'),
             '|Re eps + 1| (SPP)': A('spp_dist')}
    L = [f'{len(rows)} materials, crossed grating f = cos x cos y (3D), vacuum above, Pade, N = M = 10, 16 x 16 x 16, '
         'joint windows, omega in [1, 3].', '',
         '| feature | Spearman rho with log10 (true max error of AWE) | p-value | rho with log10 max\\|R_hyb - R_AWE\\| | p-value |',
         '|---|---|---|---|---|']
    for k, v in feats.items():
        r1, p1 = sps.spearmanr(v, y)
        r2, p2 = sps.spearmanr(v, np.log10(A('dR_max') + 1e-17))
        L.append(f'| {k} | {r1:.2f} | {p1:.0e} | {r2:.2f} | {p2:.0e} |')
    ea, eh = A('err_awe_max'), A('err_hyb_max')
    L += ['', f"True summation error on the window edges (pointwise HOPS, same grid), pooled over materials: "
          f"max AWE {ea.max():.1e} -> hybrid {eh.max():.1e}; median of the per-material maxima AWE {np.median(ea):.1e} -> "
          f"hybrid {np.median(eh):.1e}; hybrid better (> 2x) on {100 * np.median(A('hyb_better')):.0f} % of the samples of the "
          f"median material, worse (> 2x) on {100 * np.median(A('hyb_worse')):.0f} %."]
    groups = {'lossless dielectrics': [r for r in rows if r['lossless']],
              'absorbing, Re eps > 0': [r for r in rows if not r['lossless'] and r['epsr_min'] > 0],
              'metallic somewhere (Re eps < 0)': [r for r in rows if r['epsr_min'] < 0]}
    L += ['', '| group | materials | median true max err AWE | median true max err hybrid | median max\\|R_hyb - R_AWE\\| | '
          'AWE non-physical R | hybrid non-physical R |', '|---|---|---|---|---|---|---|']
    for g, rs in groups.items():
        if rs:
            L.append(f"| {g} | {len(rs)} | {np.median([r['err_awe_max'] for r in rs]):.1e} | "
                     f"{np.median([r['err_hyb_max'] for r in rs]):.1e} | {np.median([r['dR_max'] for r in rs]):.1e} | "
                     f"{100 * np.mean([r['awe_nonphys'] for r in rs]):.2f} % | {100 * np.mean([r['hyb_nonphys'] for r in rs]):.2f} % |")
    try:
        from sklearn.tree import DecisionTreeClassifier, export_text
        X = np.column_stack(list(feats.values()))
        lab = (ea > 1e-2).astype(int)
        clf = DecisionTreeClassifier(max_depth=2, random_state=0).fit(X, lab)
        L += ['', f'Decision tree (depth 2) for "3D AWE off by more than 0.01 in R (true summation error)" '
                  f'(training accuracy {clf.score(X, lab):.0%}, {lab.sum()} of {lab.size} materials positive):',
              '```', export_text(clf, feature_names=list(feats)), '```']
    except ImportError:
        pass
    srt = sorted(rows, key=lambda r: -r['err_awe_max'])
    L += ['', '| material | n (window mean) | windows | true max err AWE -> hybrid | max \\|R_hyb - R_AWE\\| | AWE off by > 1e-6 | '
          'AWE / hybrid non-physical R | max log10\\|D\\| AWE -> hybrid | hybrid min R/R_flat | absorptance range |',
          '|---|---|---|---|---|---|---|---|---|---|']
    for r in srt:
        dd = f"{r['D_awe_max']:.1f} -> {r['D_hyb_max']:.1f}" if r['lossless'] else 'absorbing'
        ar = f"{r['A_min']:.2f}-{r['A_max']:.2f}" if not r['lossless'] else '-'
        L.append(f"| {r['key']} | {r['n_re_mean']:.2f}+{r['n_im_mean']:.2f}i | {r['windows']} | {r['err_awe_max']:.1e} -> "
                 f"{r['err_hyb_max']:.1e} | {r['dR_max']:.1e} | {100 * r['frac_dR6']:.0f} % | "
                 f"{100 * r['awe_nonphys']:.1f} % / {100 * r['hyb_nonphys']:.1f} % | {dd} | {r['Rrel_min']:.3f} | {ar} |")
    open(os.path.join(HERE, 'results', 'materials.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L[:40]))
    fig, axs = plt.subplots(1, 3, figsize=(19, 5.5))
    sc = axs[0].scatter(A('n_re_mean'), A('n_im_mean') + 1e-3, c=y, cmap='viridis', s=40)
    for r in rows:
        axs[0].annotate(r['key'], (r['n_re_mean'], r['n_im_mean'] + 1e-3), fontsize=6)
    axs[0].set_xscale('log'); axs[0].set_yscale('log')
    axs[0].set_xlabel('Re n (window mean)'); axs[0].set_ylabel('Im n (+1e-3)')
    fig.colorbar(sc, ax=axs[0], label='log10 true max error of 3D AWE')
    axs[0].set_title('where 3D HOPS/AWE needs the PINN correction')
    axs[1].loglog(ea + 1e-17, eh + 1e-17, 'o', ms=5)
    for r in rows:
        axs[1].annotate(r['key'], (r['err_awe_max'] + 1e-17, r['err_hyb_max'] + 1e-17), fontsize=6)
    lim = [1e-16, 1]
    axs[1].plot(lim, lim, 'k:', lw=1)
    axs[1].set_xlabel('true max error, 3D HOPS/AWE'); axs[1].set_ylabel('true max error, 3D HOPS/AWE + PINN')
    axs[1].set_title('per material (below the diagonal: hybrid better)')
    top = srt[:25]
    yy = np.arange(len(top))
    axs[2].barh(yy + 0.2, [max(r['err_awe_max'], 1e-17) for r in top], 0.4, log=True, color='C0', label='AWE')
    axs[2].barh(yy - 0.2, [max(r['err_hyb_max'], 1e-17) for r in top], 0.4, log=True, color='C3', label='AWE + PINN')
    axs[2].set_yticks(yy); axs[2].set_yticklabels([r['key'] for r in top], fontsize=7); axs[2].invert_yaxis()
    axs[2].set_xlabel('true max |R error| (window edges)'); axs[2].legend(fontsize=8)
    axs[2].set_title('25 materials where 3D AWE is least accurate')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'materials_overview.png'), dpi=105)
    plt.close(fig)


def gallery(keys=('Ag', 'Au', 'Al', 'TiN', 'Si', 'TiO2', 'SiC', 'water')):
    """hybrid R/R_flat and log10 |R_hyb - R_AWE| maps for a selection of materials (saved arrays)"""
    keys = [k for k in keys if os.path.exists(os.path.join(OUT, f'{k}.npz'))]
    fig, axs = plt.subplots(int(np.ceil(len(keys) / 2)), 4, figsize=(20, 3.7 * np.ceil(len(keys) / 2)), squeeze=False)
    for n_, key in enumerate(keys):
        d = np.load(os.path.join(OUT, f'{key}.npz'))
        per = float(d['period'])
        nwin = len({k.split('_')[0] for k in d.files if k.startswith('w')})
        axR, axE = axs[n_ // 2, 2 * (n_ % 2)], axs[n_ // 2, 2 * (n_ % 2) + 1]
        Rr_all = np.concatenate([(np.real(d[f'w{k}_ru']) / np.real(d[f'w{k}_ru_flat'])).ravel() for k in range(nwin)])
        lo = np.percentile(Rr_all, 1)
        for k in range(nwin):
            lam = d[f'w{k}_lam'] * per / (2 * np.pi)
            Rr = np.real(d[f'w{k}_ru']) / np.real(d[f'w{k}_ru_flat'])
            m1 = axR.pcolormesh(lam, d[f'w{k}_Eps'], Rr, cmap='hot', shading='auto', vmin=lo, vmax=1.0)
            m2 = axE.pcolormesh(lam, d[f'w{k}_Eps'], np.log10(np.abs(np.real(d[f'w{k}_ru']) - np.real(d[f'w{k}_ru_awe'])) + 1e-17),
                                cmap='viridis', vmin=-14, vmax=0, shading='auto')
        axR.set_title(f'{key} (3D, cos x cos y): HOPS/AWE + PINN, R/R_flat', fontsize=10)
        axE.set_title(f'{key}: log10 |R_hyb - R_AWE|', fontsize=10)
        for ax in (axR, axE):
            ax.set_xlabel(r'$\lambda$ ($\mu$m)'); ax.set_ylabel(r'$\varepsilon$')
        fig.colorbar(m1, ax=axR); fig.colorbar(m2, ax=axE)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'materials_gallery.png'), dpi=90)
    plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--keys', nargs='*')
    ap.add_argument('--procs', type=int, default=2)
    ap.add_argument('--analyse-only', action='store_true')
    ap.add_argument('--gallery', action='store_true')
    a = ap.parse_args()
    keys = a.keys or [k for k, v in mat.MATERIALS.items() if v.get('category') != 'superstrate']
    if a.analyse_only:
        rows = [json.load(open(os.path.join(OUT, f'{k}.json'))) for k in keys if os.path.exists(os.path.join(OUT, f'{k}.json'))]
    elif a.procs > 1:
        with ProcessPoolExecutor(max_workers=a.procs) as ex:
            rows = list(ex.map(one, keys))
    else:
        rows = [one(k) for k in keys]
    analyse(rows)
    if a.gallery:
        gallery()
