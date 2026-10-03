"""materials_hybrid.py -- what do the reflectivity maps and energy defects of the curated material library
(hops/materials.py, 64 refractiveindex.info materials) teach us, once HOPS/AWE + PINN is available as an
accurate reference and the PINN indicator as a reliability map?

Setup (as material_survey.py, but f = cos x, which Nx = 32 resolves -- with cos 4x on Nx = 32 the
discretisation error dominates for high-index and metallic substrates, see results/materials_cos4x_nx32.md):
vacuum over the material, TM, eps_max = 0.2, Pade, the
material's suggested grating period, the index evaluated at the centre of each band q = 1..6, BUT one AWE
window per band ('paper' windows -- no splitting at the substrate's Wood anomalies), 15 x 15 per band,
N = M = 12.  For every material:
  features   Re n, Im n, Re eps = Re n^2, lossless?, number of lower-layer Wood anomalies inside the bands,
             distance of Re eps to the surface-plasmon condition eps = -1
  AWE        fraction of the map where |R_AWE - R_hybrid| > 1e-6 (> 1e-3), non-physical R (< 0 or > 1),
             fraction where the indicator exceeds 1e-10,
             max log10|D| (lossless materials)
  hybrid     the same map from the physics-informed weights; max |R_hyb - R_AWE|; max log10|D|
  physics    hybrid: min R/R_flat (resonance depth), spread of R/R_flat, absorptance range
Outputs results/materials.csv, results/materials.md (with Spearman correlations and a small decision tree),
figures/materials/*.png.

    python materials_hybrid.py                # all 64 materials, 2 processes (~40 min)
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
from awepinn import core                    # noqa: E402
from awepinn.core import rm                 # noqa: E402
from hops import materials as mat           # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PROFILE = 'cosx'
OUT = os.path.join(HERE, 'results', 'materials')
FIG = os.path.join(HERE, 'figures', 'materials')


def wood_count(res):
    n = 0
    for r in res:
        nw = complex(r['n_w'])
        if abs(nw.imag) > 1e-3 * abs(nw) or nw.real <= 0:
            continue
        lo, hi = r['omega'].min() * nw.real, r['omega'].max() * nw.real
        n += sum(1 for p in range(-60, 61) if lo < abs(p) < hi)
    return n


def one(key, n=15, M=12):
    fn = os.path.join(OUT, f'{key}.json')
    if os.path.exists(fn):
        return json.load(open(fn))
    warnings.filterwarnings('ignore')
    t0 = time.time()
    try:
        res, info = core.run('silver', N_Eps=n, N_delta=n, workers=1, verbose=False, n_w=key, profile=PROFILE, M=M,
                             period=mat.MATERIALS[key]['period'], relative=True, windows='paper')
    except Exception as e:                    # noqa: BLE001
        return dict(key=key, error=repr(e))
    cat = lambda k: np.concatenate([np.real(r[k]).ravel() for r in res])
    Ra, Rh, Da, Dh, ind, Rf = cat('ru_awe'), cat('ru'), cat('ee_awe'), cat('ee'), cat('indicator'), cat('ru_flat')
    ns = np.array([complex(r['n_w']) for r in res])
    lossless = bool(info['lossless'])
    epsr = (ns ** 2).real
    s = dict(key=key, category=mat.MATERIALS[key].get('category', ''), period=mat.MATERIALS[key]['period'],
             n_re_mean=float(ns.real.mean()), n_im_mean=float(ns.imag.mean()), n_im_max=float(ns.imag.max()),
             epsr_min=float(epsr.min()), epsr_max=float(epsr.max()),
             spp_dist=float(np.min(np.abs(epsr + 1))), lossless=lossless, wood=wood_count(res),
             unreliable=float(np.mean(ind > 1e-10)), solved=float(np.mean(np.concatenate([r['solved'].ravel() for r in res]))),
             awe_nonphys=float(np.mean((Ra < -1e-6) | (Ra > 1 + 1e-6))),
             hyb_nonphys=float(np.mean((Rh < -1e-6) | (Rh > 1 + 1e-6))),
             dR_max=float(np.max(np.abs(Rh - Ra))), dR_median=float(np.median(np.abs(Rh - Ra))),
             frac_dR6=float(np.mean(np.abs(Rh - Ra) > 1e-6)), frac_dR3=float(np.mean(np.abs(Rh - Ra) > 1e-3)),
             hyb_D_bad=float(np.mean(np.abs(Dh) > 1e-6)) if lossless else np.nan,
             awe_D_bad=float(np.mean(np.abs(Da) > 1e-6)) if lossless else np.nan,
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
    print(f"{key:12s} n~{s['n_re_mean']:.2f}+{s['n_im_mean']:.2f}i wood {s['wood']:3d} unreliable {100 * s['unreliable']:3.0f} % "
          f"max|dR| {s['dR_max']:.1e} nonphys AWE {100 * s['awe_nonphys']:.1f} % -> {100 * s['hyb_nonphys']:.1f} % "
          f"({s['time']:.0f} s)", flush=True)
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
    y = np.log10(A('dR_max') + 1e-17)
    feats = {'Re n': A('n_re_mean'), 'Im n': A('n_im_mean'), 'Re eps': A('epsr_min'), 'Wood anomalies in bands': A('wood'),
             'lossless': A('lossless'), '|Re eps + 1| (SPP)': A('spp_dist')}
    L = ['| feature | Spearman rho with log10 max\\|R_hyb - R_AWE\\| | p-value | Spearman rho with the fraction of the map where AWE is off by > 1e-6 | p-value |',
         '|---|---|---|---|---|']
    for k, v in feats.items():
        r1, p1 = sps.spearmanr(v, y)
        r2, p2 = sps.spearmanr(v, A('frac_dR6'))
        L.append(f'| {k} | {r1:.2f} | {p1:.0e} | {r2:.2f} | {p2:.0e} |')
    groups = {'lossless dielectrics': [r for r in rows if r['lossless']],
              'absorbing, Re eps > 0': [r for r in rows if not r['lossless'] and r['epsr_min'] > 0],
              'metallic somewhere (Re eps < 0)': [r for r in rows if r['epsr_min'] < 0]}
    L += ['', '| group | materials | median fraction with AWE off by > 1e-6 | median max\\|R_hyb - R_AWE\\| | AWE non-physical R | hybrid non-physical R |',
          '|---|---|---|---|---|---|']
    for g, rs in groups.items():
        if rs:
            L.append(f"| {g} | {len(rs)} | {100 * np.median([r['frac_dR6'] for r in rs]):.0f} % | "
                     f"{np.median([r['dR_max'] for r in rs]):.1e} | {100 * np.mean([r['awe_nonphys'] for r in rs]):.1f} % | "
                     f"{100 * np.mean([r['hyb_nonphys'] for r in rs]):.1f} % |")
    try:
        from sklearn.tree import DecisionTreeClassifier, export_text
        X = np.column_stack(list(feats.values()))
        lab = (A('dR_max') > 1e-2).astype(int)            # AWE off by more than 0.01 in R somewhere on the map
        clf = DecisionTreeClassifier(max_depth=2, random_state=0).fit(X, lab)
        L += ['', f'Decision tree (depth 2) for "AWE off by more than 0.01 in R somewhere on the map" '
                  f'(training accuracy {clf.score(X, lab):.0%}, {lab.sum()} of {lab.size} materials positive):',
              '```', export_text(clf, feature_names=list(feats)), '```']
    except ImportError:
        pass
    srt = sorted(rows, key=lambda r: -r['dR_max'])
    L += ['', '| material | n (band mean) | Wood anomalies | AWE off by > 1e-6 | max \\|R_hyb - R_AWE\\| | AWE / hybrid non-physical R | '
          'max log10\\|D\\| AWE -> hybrid | hybrid: min R/R_flat (resonance) | absorptance range |', '|---|---|---|---|---|---|---|---|---|']
    for r in srt:
        dd = f"{r['D_awe_max']:.1f} -> {r['D_hyb_max']:.1f}" if r['lossless'] else 'absorbing'
        ar = f"{r['A_min']:.2f}-{r['A_max']:.2f}" if not r['lossless'] else '-'
        L.append(f"| {r['key']} | {r['n_re_mean']:.2f}+{r['n_im_mean']:.2f}i | {r['wood']} | {100 * r['frac_dR6']:.0f} % | "
                 f"{r['dR_max']:.1e} | {100 * r['awe_nonphys']:.1f} % / {100 * r['hyb_nonphys']:.1f} % | {dd} | "
                 f"{r['Rrel_min']:.3f} | {ar} |")
    open(os.path.join(HERE, 'results', 'materials.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L[:30]))
    # figures
    fig, axs = plt.subplots(1, 3, figsize=(19, 5.5))
    sc = axs[0].scatter(A('n_re_mean'), A('n_im_mean') + 1e-3, c=y, cmap='viridis', s=40)
    for r in rows:
        axs[0].annotate(r['key'], (r['n_re_mean'], r['n_im_mean'] + 1e-3), fontsize=6)
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_xlabel('Re n (band mean)')
    axs[0].set_ylabel('Im n (+1e-3)')
    fig.colorbar(sc, ax=axs[0], label='log10 max |R_hyb - R_AWE|')
    axs[0].set_title('where HOPS/AWE needs the PINN correction')
    axs[1].scatter(A('wood'), A('frac_dR6') * 100, c=A('lossless'), cmap='coolwarm', s=40)
    for r in rows:
        axs[1].annotate(r['key'], (r['wood'], 100 * r['frac_dR6']), fontsize=6)
    axs[1].set_xlabel('number of lower-layer Wood anomalies inside the bands')
    axs[1].set_ylabel('% of the map where AWE is off by > 1e-6')
    axs[1].set_title('blue: absorbing, red: lossless')
    srt = srt[:25]
    yy = np.arange(len(srt))
    axs[2].barh(yy, [max(r['dR_max'], 1e-17) for r in srt], log=True, color='C3')
    axs[2].set_yticks(yy)
    axs[2].set_yticklabels([r['key'] for r in srt], fontsize=7)
    axs[2].invert_yaxis()
    axs[2].set_xlabel('max |R_hyb - R_AWE| over the six-band map')
    axs[2].set_title('25 materials where AWE and AWE + PINN differ most')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'materials_overview.png'), dpi=105)
    plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--keys', nargs='*')
    ap.add_argument('--procs', type=int, default=2)
    ap.add_argument('--analyse-only', action='store_true')
    ap.add_argument('--gallery', action='store_true')
    a = ap.parse_args()
    keys = a.keys or list(mat.MATERIALS)
    if a.analyse_only:
        rows = [json.load(open(os.path.join(OUT, f'{k}.json'))) for k in keys if os.path.exists(os.path.join(OUT, f'{k}.json'))]
    elif a.procs > 1:
        with ProcessPoolExecutor(max_workers=a.procs) as ex:
            rows = list(ex.map(one, keys))
    else:
        rows = [one(k) for k in keys]
    analyse(rows)


def gallery(keys=('Al', 'Ag', 'Mo', 'W', 'GaN', 'TiO2', 'ZnGeP2', 'water')):
    """hybrid R/R_flat maps and |R_hyb - R_AWE| for a selection of materials (from the saved arrays)"""
    fig, axs = plt.subplots(4, 4, figsize=(20, 15))
    for n_, key in enumerate(keys):
        fn = os.path.join(OUT, f'{key}.npz')
        if not os.path.exists(fn):
            continue
        d = np.load(fn)
        per = float(d['period'])
        nwin = len({k.split('_')[0] for k in d.files if k.startswith('w')})
        axR, axE = axs[(n_ // 2) * 1 + 0 if False else n_ // 2, 2 * (n_ % 2)], axs[n_ // 2, 2 * (n_ % 2) + 1]
        for k in range(nwin):
            lam = d[f'w{k}_lam'] * per / (2 * np.pi)
            Rr = np.real(d[f'w{k}_ru']) / np.real(d[f'w{k}_ru_flat'])
            m1 = axR.pcolormesh(lam, d[f'w{k}_Eps'], Rr, cmap='hot', shading='auto',
                                vmin=0.85, vmax=1.0)
            m2 = axE.pcolormesh(lam, d[f'w{k}_Eps'], np.log10(np.abs(np.real(d[f'w{k}_ru']) - np.real(d[f'w{k}_ru_awe'])) + 1e-17),
                                cmap='viridis', vmin=-14, vmax=0, shading='auto')
        axR.set_title(f'{key}: HOPS/AWE + PINN, R/R_flat', fontsize=10)
        axE.set_title(f'{key}: log10 |R_hyb - R_AWE|', fontsize=10)
        for ax in (axR, axE):
            ax.set_xlabel(r'$\lambda$ ($\mu$m)')
            ax.set_ylabel(r'$\varepsilon$')
        fig.colorbar(m1, ax=axR)
        fig.colorbar(m2, ax=axE)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'materials_gallery.png'), dpi=90)
    plt.close(fig)


if __name__ == '__main__' and '--gallery' in sys.argv:
    gallery()
