"""refl_map_hybrid.py -- the HOPS/AWE + PINN version of refl_map.py (one scenario, plotted).

Works like HOPS_Python/refl_map.py: pick a scenario (edit SCENARIO below, or pass --scenario) and run.
It computes the standard HOPS/AWE map AND the HOPS/AWE + PINN (physics-informed summation) map from the
same coefficient fields and saves, in figures/refl_map_hybrid/:
    refl_map_<tag>_awe_R.png,    refl_map_<tag>_awe_D.png       standard HOPS/AWE (as refl_map.py draws it)
    refl_map_<tag>_hybrid_R.png, refl_map_<tag>_hybrid_D.png    HOPS/AWE + PINN, same rendering
    refl_map_<tag>_compare.png                                  side by side, |R_hyb - R_AWE|, indicator
    refl_map_<tag>_hybrid.npz                                   the arrays

    python refl_map_hybrid.py                                   # SCENARIO below (gold), 21 x 21 per band
    python refl_map_hybrid.py --scenario dielectric --neps 31 --ndelta 31
    python refl_map_hybrid.py --scenario silver --Nx 64         # resolved grid for the metal (slower)
    python refl_map_hybrid.py --scenario thesis22 --q 1         # only frequency band q = 1
    python refl_map_hybrid.py --nw Au --nu water                # any materials (refl_map.py syntax)
    python refl_map_hybrid.py --list-scenarios

Cost: the least-squares step is ~0.05-0.5 s per corrected point, so keep the grid modest
(21 x 21 x 6 bands ~ 2-5 min with --workers 4; the paper's 100 x 100 would take hours).
"""
import argparse
import os
import sys
import time
import warnings

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')     # one BLAS thread per process: the bands run in parallel instead
import numpy as np
import matplotlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from awepinn import core                 # noqa: E402
from awepinn.core import rm              # noqa: E402

SCENARIO = 'gold'                        # <- edit to run another scenario from the IDE
N_EPS = 21                               # points in eps per band
N_DELTA = 21                             # points in frequency per band
WORKERS = max(1, min(4, (os.cpu_count() or 2) - 1))   # bands in parallel processes
SHOW = True                              # open the figures when done
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'refl_map_hybrid')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scenario', default=SCENARIO, choices=list(rm.SCENARIOS))
    ap.add_argument('--list-scenarios', action='store_true')
    ap.add_argument('--neps', type=int, default=N_EPS)
    ap.add_argument('--ndelta', type=int, default=N_DELTA)
    ap.add_argument('--q', type=int, nargs='*', help='frequency bands (default: the scenario\'s)')
    ap.add_argument('--tol', type=float, default=core.DEFAULT_TOL,
                    help='indicator threshold; inf = AWE + indicator only (default 1e-14)')
    ap.add_argument('--basis-order', type=int, default=core.DEFAULT_ORDER,
                    help='orders n, m <= this in the least-squares basis (use N for PEC-like |n_w| >= 20)')
    ap.add_argument('--workers', type=int, default=WORKERS)
    ap.add_argument('--tag', help='name for the output files (default: scenario)')
    ap.add_argument('--no-show', action='store_true')
    g = ap.add_argument_group('overrides of the scenario (same meaning as in refl_map.py)')
    g.add_argument('--nu', help='upper index: number (1.5, 0.05+2.275i) or material key')
    g.add_argument('--nw', help='lower index: number or material key')
    g.add_argument('--mode', choices=['TE', 'TM'])
    g.add_argument('--profile')
    g.add_argument('--eps-max', type=float)
    g.add_argument('--alpha', type=float)
    g.add_argument('--period', type=float)
    g.add_argument('--M', type=int)
    g.add_argument('--N', type=int)
    g.add_argument('--Nx', type=int)
    g.add_argument('--Nz', type=int)
    g.add_argument('--summation', choices=['taylor', 'pade'])
    a = ap.parse_args(argv)

    if a.list_scenarios:
        w = max(map(len, rm.SCENARIOS))
        for k, v in rm.SCENARIOS.items():
            print(f'{k:{w}s}  {v["desc"]}')
        return
    if a.no_show or not SHOW:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    backend = matplotlib.get_backend()
    import run_scenarios as rs                       # compare figure + statistics (it selects Agg ...)
    if SHOW and not a.no_show:
        plt.switch_backend(backend)                  # ... so restore the interactive backend

    over = {k: v for k, v in dict(n_u=a.nu, n_w=a.nw, mode=a.mode, profile=a.profile, eps_max=a.eps_max,
                                  alpha=a.alpha, period=a.period, M=a.M, N=a.N, Nx=a.Nx, Nz=a.Nz,
                                  taylor=None if a.summation is None else a.summation == 'taylor').items()
            if v is not None}
    tag = a.tag or (a.scenario + ''.join(f'_{k}{v}' for k, v in over.items()).replace('+', 'p').replace('.', 'd'))
    print(f'{a.scenario}: {a.neps} x {a.ndelta} per band, tol {a.tol:g}, basis order {a.basis_order}, '
          f'{a.workers} workers ...', flush=True)
    t0 = time.time()
    warnings.filterwarnings('ignore')
    res, info = core.run(a.scenario, qq=a.q, N_Eps=a.neps, N_delta=a.ndelta, tol=a.tol,
                         basis_order=a.basis_order, workers=a.workers, verbose=True, **over)
    s = rs.stats(res, info)
    s.update(name=tag, grid=a.neps, desc=info.get('desc', ''))
    info['stats'] = s

    os.makedirs(OUTDIR, exist_ok=True)
    awe = [dict(r, ru=r['ru_awe'], rl=r['rl_awe'], ee=r['ee_awe'], RR=r['RR_awe']) for r in res]
    saved = list(rm.plot(awe, info, outdir=OUTDIR, tag=f'{tag}_awe') or [])
    saved += list(rm.plot(res, info, outdir=OUTDIR, tag=f'{tag}_hybrid') or [])
    cmp_png = os.path.join(OUTDIR, f'refl_map_{tag}_compare.png')
    rs.compare_figure(tag, res, info, cmp_png)
    saved.append(cmp_png)
    np.savez_compressed(os.path.join(OUTDIR, f'refl_map_{tag}_hybrid.npz'),
                        **{f'w{k}_{q}': np.asarray(r[q]) for k, r in enumerate(res)
                           for q in ('lam', 'omega', 'Eps', 'ru', 'rl', 'ee', 'ru_awe', 'ee_awe', 'ru_flat',
                                     'indicator', 'solved')})
    for fn in saved:
        print('saved', fn)
    print(f"max |R_hyb - R_AWE| = {s['dR_max']:.1e},  {100 * s['solved_frac']:.0f} % of points corrected" +
          (f",  max log10|D|: AWE {s['D_awe_max_log10']:.1f} -> hybrid {s['D_hyb_max_log10']:.1f}"
           if s['lossless'] else ''))
    print(f'total {time.time() - t0:.1f} s')
    if SHOW and not a.no_show:
        plt.show()


if __name__ == '__main__':
    main()
