"""refl_map_hybrid_3D.py -- the 3D HOPS/AWE + PINN version of refl_map_3D.py (one scenario, plotted).

Works like HOPS_Python/refl_map_3D.py: pick a scenario (edit SCENARIO below, or pass --scenario) and run.
It computes the standard 3D HOPS/AWE map AND the HOPS/AWE + PINN map from the same coefficient fields and saves,
in figures/refl_map_hybrid/:
    refl_map_3D_<tag>_awe_R.png,    refl_map_3D_<tag>_awe_D.png       standard 3D HOPS/AWE
    refl_map_3D_<tag>_hybrid_R.png, refl_map_3D_<tag>_hybrid_D.png    3D HOPS/AWE + PINN, same rendering
    refl_map_3D_<tag>_compare.png                                     side by side, |R_hyb - R_AWE|, indicator
    refl_map_3D_<tag>_hybrid.npz                                      the arrays

    python refl_map_hybrid_3D.py                                   # SCENARIO below (dielectric), 15 x 15 per window
    python refl_map_hybrid_3D.py --scenario gold --q 1 2           # only the first two bands (quick)
    python refl_map_hybrid_3D.py --scenario Au_crossed --neps 11 --ndelta 11
    python refl_map_hybrid_3D.py --nw Au --period 0.8 --profile cosx+cosy --q 1
    python refl_map_hybrid_3D.py --list-scenarios

Cost: ~0.05 s (16 x 16 grids) to ~0.4 s (32 x 32) per corrected point; keep the grid modest.
"""
import argparse
import os
import sys
import time
import warnings

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
import matplotlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from awepinn3d import core                 # noqa: E402
from awepinn3d.core import rm3             # noqa: E402

SCENARIO = 'dielectric'                    # <- edit to run another scenario from the IDE
N_EPS = 15
N_DELTA = 15
WORKERS = max(1, min(2, (os.cpu_count() or 2) - 1))   # windows in parallel (each 32x32 window needs ~1.5 GB)
SHOW = True
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'refl_map_hybrid')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scenario', default=SCENARIO, choices=list(rm3.SCENARIOS))
    ap.add_argument('--list-scenarios', action='store_true')
    ap.add_argument('--neps', type=int, default=N_EPS)
    ap.add_argument('--ndelta', type=int, default=N_DELTA)
    ap.add_argument('--q', type=int, nargs='*', help="frequency bands (default: the scenario's)")
    ap.add_argument('--tol', type=float, default=core.DEFAULT_TOL, help='indicator threshold (default 1e-14)')
    ap.add_argument('--basis-order', type=int, default=core.DEFAULT_ORDER)
    ap.add_argument('--workers', type=int, default=WORKERS)
    ap.add_argument('--tag')
    ap.add_argument('--no-show', action='store_true')
    g = ap.add_argument_group('overrides of the scenario (same meaning as in refl_map_3D.py)')
    for k in ('nu', 'nw', 'profile'):
        g.add_argument(f'--{k}')
    g.add_argument('--mode', choices=['TE', 'TM'])
    for k in ('eps-max', 'alpha', 'beta', 'theta', 'phi', 'period', 'max-delta'):
        g.add_argument(f'--{k}', type=float)
    for k in ('M', 'N', 'Nx', 'Ny', 'Nz'):
        g.add_argument(f'--{k}', type=int)
    g.add_argument('--windows', choices=['paper', 'joint'])
    g.add_argument('--summation', choices=['taylor', 'pade'])
    a = ap.parse_args(argv)
    if a.list_scenarios:
        w = max(map(len, rm3.SCENARIOS))
        for k, v in rm3.SCENARIOS.items():
            print(f'{k:{w}s}  {v["desc"]}')
        return
    if a.no_show or not SHOW:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    backend = matplotlib.get_backend()
    import run_scenarios as rs                       # compare figure + statistics (it selects Agg ...)
    if SHOW and not a.no_show:
        plt.switch_backend(backend)
    over = dict(n_u=a.nu, n_w=a.nw, profile=a.profile, mode=a.mode, eps_max=a.eps_max, alpha=a.alpha, beta=a.beta,
                theta=a.theta, phi=a.phi, period=a.period, max_delta=a.max_delta, M=a.M, N=a.N, Nx=a.Nx, Ny=a.Ny,
                Nz=a.Nz, windows=a.windows, taylor=None if a.summation is None else a.summation == 'taylor')
    over = {k: v for k, v in over.items() if v is not None}
    tag = a.tag or (a.scenario + ''.join(f'_{k}{v}' for k, v in over.items()).replace('+', 'p').replace('.', 'd'))
    print(f'{a.scenario}: {a.neps} x {a.ndelta} per window, tol {a.tol:g}, basis order {a.basis_order}, '
          f'{a.workers} workers ...', flush=True)
    t0 = time.time()
    warnings.filterwarnings('ignore')
    res, info = core.run(a.scenario, qq=a.q, N_Eps=a.neps, N_delta=a.ndelta, tol=a.tol, basis_order=a.basis_order,
                         workers=a.workers, verbose=True, **over)
    s = rs.stats(res, info)
    s.update(name=tag, grid=a.neps, desc=info.get('desc', ''))
    info['stats'] = s
    os.makedirs(OUTDIR, exist_ok=True)
    saved = list(rm3.plot(core.awe_view(res), info, outdir=OUTDIR, tag=f'{tag}_awe') or [])
    saved += list(rm3.plot(res, info, outdir=OUTDIR, tag=f'{tag}_hybrid') or [])
    cmp_png = os.path.join(OUTDIR, f'refl_map_3D_{tag}_compare.png')
    rs.compare_figure(tag, res, info, cmp_png)
    saved.append(cmp_png)
    np.savez_compressed(os.path.join(OUTDIR, f'refl_map_3D_{tag}_hybrid.npz'),
                        **{f'w{k}_{q}': np.asarray(r[q]) for k, r in enumerate(res) for q in rs.KEYS})
    for fn in saved:
        print('saved', fn)
    print(f"{len(res)} windows; max |R_hyb - R_AWE| = {s['dR_max']:.1e},  {100 * s['solved_frac']:.0f} % of points corrected" +
          (f",  max log10|D|: AWE {s['D_awe_max_log10']:.1f} -> hybrid {s['D_hyb_max_log10']:.1f}" if s['lossless'] else ''))
    print(f'total {time.time() - t0:.1f} s')
    if SHOW and not a.no_show:
        plt.show()


if __name__ == '__main__':
    main()
