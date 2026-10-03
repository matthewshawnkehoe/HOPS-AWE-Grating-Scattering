"""material_survey.py -- screen every material of the curated refractiveindex.info library
(hops/materials.py) with a fast, reduced-resolution reflectivity map and rank them.

For each material (lower layer, vacuum above, TM, f = cos(4x), eps_max = 0.2, Pade, the
material's suggested grating period, dispersion evaluated per band) it records
  Rmin      1st percentile of R/R_flat over the map (depth of grating resonances; 1 = featureless)
  Rspread   5-95 % spread of R/R_flat
  Rabs      1-99 % range of the absolute reflectivity R
  bad       fraction of (eps, delta) points with a non-physical R (< 0 or > 1.05) -> Pade trouble
  Dsmall    median log10|D| for eps <= 0.05, Dmed over the whole map (lossless materials only)
  absorb    median absorptance D = 1 - R (absorbing substrates, no transmitted propagating modes)
Windows: 'joint' (bands split at the lower layer's Wood anomalies, see refl_map.py --windows).
and writes figures/material_survey/material_survey.csv plus a gallery figures/material_survey/material_survey.png.

    python material_survey.py                 # all curated materials (~10-20 min)
    python material_survey.py --keys Ag Au Al Na TiN ITO --M 15 --n 60
"""
import argparse
import csv
import os
import time
import warnings

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import refl_map as rm
from hops import materials as mat
from hops.plotting import matlab_contourf
from matplotlib.colors import Normalize

OUT = os.path.join(os.path.dirname(rm.OUTDIR), 'material_survey')


def survey_one(key, M=12, n=40, profile='cos4x', workers=1, n_u=1.0, windows='joint'):
    warnings.simplefilter('ignore')
    t0 = time.time()
    res, info = rm.run('silver', qq=None, N_Eps=n, N_delta=n, verbose=False, workers=workers,
                       n_u=n_u, n_w=key, profile=profile, M=M, period=mat.MATERIALS[key]['period'],
                       relative=True, windows=windows)
    RR = np.concatenate([np.real(r['RR']).ravel() for r in res])
    Rabs = np.concatenate([np.real(r['ru']).ravel() for r in res])
    Dall = np.concatenate([np.abs(r['ee']).ravel() for r in res])
    D5 = np.concatenate([np.abs(r['ee'][r['Eps'] <= 0.05]).ravel() for r in res])
    Dsmall = np.log10(np.nanmedian(D5) + 1e-300)
    fin = np.isfinite(RR)
    bad = np.mean(~fin | (Rabs < -1e-3) | (Rabs > 1.05))
    row = dict(key=key, category=mat.MATERIALS[key]['category'], period=info['period'],
               n_q1=f"{complex(res[0]['n_w']):.3g}", n_q6=f"{complex(res[-1]['n_w']):.3g}",
               windows=len(res),
               Rmin=np.nanpercentile(RR[fin], 1), Rspread=np.nanpercentile(RR[fin], 95) - np.nanpercentile(RR[fin], 5),
               Rabs_min=np.nanpercentile(Rabs[fin], 1), Rabs_max=np.nanpercentile(Rabs[fin], 99), bad=bad,
               Dsmall=Dsmall if info['lossless'] else np.nan,
               Dmed=np.log10(np.nanmedian(Dall) + 1e-300) if info['lossless'] else np.nan,
               absorb=np.nan if info['lossless'] else np.nanmedian(np.real(Dall)),
               lossless=info['lossless'],
               seconds=time.time() - t0)
    return row, res, info


def gallery(items, fn, ncol=4):
    nrow = int(np.ceil(len(items) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow), squeeze=False)
    for ax, (row, res, info) in zip(axs.ravel(), items):
        Zs = [np.clip(np.real(r['RR']), 0, 1.0) for r in res]
        allz = np.concatenate([z.ravel() for z in Zs])
        norm = Normalize(np.nanmin(allz), 1.0)
        for r, Z in zip(res, Zs):
            matlab_contourf(ax, r['lam'] * info['period'] / (2 * np.pi), r['Eps'], Z, 'hot', norm)
        ax.set_title(f"{row['key']} ({row['category']})\nP = {info['period']:g} um, "
                     f"min R/R0 = {row['Rmin']:.2f}", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_xlabel(r'$\lambda$ ($\mu$m)', fontsize=7)
    for ax in axs.ravel()[len(items):]:
        ax.axis('off')
    fig.tight_layout()
    fig.savefig(fn, dpi=90)
    return fn


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--keys', nargs='*', help='materials (default: all non-superstrate curated keys)')
    ap.add_argument('--M', type=int, default=12)
    ap.add_argument('--n', type=int, default=40, help='grid N_eps = N_delta')
    ap.add_argument('--profile', default='cos4x')
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    a = ap.parse_args()
    keys = a.keys or [k for k, v in mat.MATERIALS.items() if v['category'] != 'superstrate']
    rows, items = [], []
    for k in keys:
        try:
            row, res, info = survey_one(k, a.M, a.n, a.profile, a.workers)
        except Exception as exc:          # keep going
            print(f'{k}: FAILED ({exc})')
            continue
        rows.append(row)
        items.append((row, res, info))
        print(f"{k:13s} {row['category']:22s} P={row['period']:<5g} n(q1)={row['n_q1']:>14s} "
              f"Rmin={row['Rmin']:.3f} spread={row['Rspread']:.3f} R=[{row['Rabs_min']:.3f},{row['Rabs_max']:.3f}] "
              f"bad={row['bad']:.3f} D(eps<.05)={row['Dsmall']:.1f} D={row['Dmed']:.1f} A={row['absorb']:.2f} "
              f"win={row['windows']} ({row['seconds']:.0f}s)", flush=True)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'material_survey.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print('saved', gallery(items, os.path.join(OUT, 'material_survey.png')))
