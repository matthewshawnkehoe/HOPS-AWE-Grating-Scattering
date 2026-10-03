"""refl_map_3D.py -- 3D analogue of refl_map.py / refl_map.m: Reflectivity Map R(eps, delta) and energy
defect D of a DOUBLY PERIODIC (crossed) grating between two layers, computed with the joint
boundary/frequency (HOPS/AWE) expansion

        R(eps, delta) = sum_{n,m} R_{n,m} eps^n delta^m,   omega = omega_bar (1 + delta),

i.e. one field solve per frequency window, then Taylor/Pade summation for every (eps, omega).

Quick start (PyCharm: put the options in the run configuration)
  python refl_map_3D.py                                  # 'silver': crossed cos(4x)cos(4y) grating
  python refl_map_3D.py --scenario dielectric            # crossed analogue of paper Fig. 9
  python refl_map_3D.py --scenario silver_1d             # y-invariant: reproduces the 2D Fig. 10a
  python refl_map_3D.py --list-scenarios
  python refl_map_3D.py --nw Au --period 0.8 --profile cosx+cosy        # any material (refractiveindex.info)
  python refl_map_3D.py --nw 1.5 --theta 10 --phi 45     # oblique incidence at fixed angles (deg)
  python refl_map_3D.py --nw 2.5 --windows joint --q 1 2  # dielectric, split at lower-layer anomalies

Geometry and units: periods d_x = d_y = 2 pi, c0 = 1, interface z = eps f(x, y), artificial boundaries
z = a and z = -b.  Incident wave exp(i(alpha x + beta y - gamma^u z)).  Windows: see hops3d/windows.py
(the 3D bands lie between consecutive Rayleigh frequencies |(alpha+p, beta+q)| / n^u: 1, sqrt 2, 2,
sqrt 5, ... for normal incidence -- the 2D bands [q, q+1] refined by the crossed diffraction orders).
Named materials are evaluated at lambda = P / omega_bar of every window (P = --period in um).
Scalar model: tau^2 = (n^u/n^w)^2 ('TM', default) or 1 ('TE') in the transmission condition, as in 2D.
"""
import argparse
import os
import time
import warnings

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import hops3d as h3
from hops3d.windows import make_windows
from hops3d.setup_data import PROFILES_3D
from hops import materials as mat
from hops.plotting import safe_log10, matlab_contourf, shared_colorbar

SCENARIO = 'silver'
PlotLambda = 1
PlotRelative = 1
SHOW = True
SOLVER = 'coupled'           # 'coupled' (fast, default) | 'lean' (= two_layer_solve_fast.m ordering) | 'operator'
WORKERS = os.cpu_count()
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures_3d', 'refl_map')

_P = dict(n_u=1.0, profile='cosxcosy', eps_max=0.2, alpha=0.0, beta=0.0, theta=None, phi=0.0, M=12, Nx=16,
          Ny=16, Nz=16, a=1.0, b=1.0, mode='TM', taylor=False, q=(1, 2, 3), period=None, windows='paper',
          sigma=0.99, sum_domain='fourier')


def _s(desc, **kw):
    d = dict(_P)
    d.update(kw)
    d['desc'] = desc
    return d


SCENARIOS = {
    # ---- 3D analogues of the paper / refl_map.m ------------------------------------------------
    'silver': _s('crossed analogue of paper Fig. 10a / refl_map.m: vacuum / silver 0.05+2.275i, '
                 'f = cos(4x)cos(4y), Pade, omega in [1, 7] (q = 1..6)', n_w=0.05 + 2.275j, profile='cos4xcos4y',
                 Nx=32, Ny=32, q=(1, 2, 3, 4, 5, 6)),
    'gold': _s('crossed analogue of paper Fig. 10b: vacuum / gold 1.48+1.883i, f = cos(4x)cos(4y), Pade',
               n_w=1.48 + 1.883j, profile='cos4xcos4y', Nx=32, Ny=32, q=(1, 2, 3, 4, 5, 6)),
    'silver_lattice': _s('as silver, windows cut only at the orders cos(4x)cos(4y) excites (|(4,4)| = 5.66, 8, ...) '
                         'and |delta| <= 0.15', n_w=0.05 + 2.275j, profile='cos4xcos4y', Nx=32, Ny=32,
                         q=(1, 2, 3, 4, 5, 6), lattice='profile', max_delta=0.15),
    'dielectric': _s('crossed analogue of paper Fig. 9: vacuum / n = 1.1, f = cos(x)cos(y), Taylor',
                     n_w=1.1, taylor=True, M=14),
    'dielectric_joint': _s('as dielectric, windows also split at the lower-layer Rayleigh anomalies',
                           n_w=1.1, taylor=True, M=14, windows='joint'),
    'dielectric_alpha': _s('crossed analogue of paper Fig. 14: dielectric, oblique alpha = beta = 0.01',
                           n_w=1.1, taylor=True, M=14, alpha=0.01, beta=0.01),
    'egg_silver': _s('silver "egg-crate" (cos x + cos y + cos x cos y)/3, Pade', n_w=0.05 + 2.275j, profile='egg'),
    'sum_gold': _s('gold, f = (cos 4x + cos 4y)/2 (two crossed 1D gratings), Pade', n_w=1.48 + 1.883j,
                   profile='cos4x+cos4y', Nx=32, Ny=32),
    # ---- 2D reductions (y-invariant gratings, Ny = 1): must reproduce refl_map.py ----------------
    'silver_1d': _s('y-invariant cos(4x), Ny = 1: identical to the 2D refl_map.py silver (paper Fig. 10a)',
                    n_w=0.05 + 2.275j, profile='cos4x', Ny=1, Nx=32, Nz=32, M=15, q=(1, 2, 3, 4, 5, 6)),
    'dielectric_1d': _s('y-invariant cos(x), Ny = 1: identical to the 2D dielectric (paper Fig. 9)',
                        n_w=1.1, profile='cosx', Ny=1, Nx=32, Nz=32, M=16, taylor=True, q=(1, 2, 3, 4, 5, 6)),
    'conical_1d': _s('y-invariant cos(x) grating under CONICAL incidence beta = 0.3 (out of the plane of '
                     'periodicity) -- beyond the 2D code', n_w=1.1, profile='cosx', Ny=1, Nx=32, Nz=32, M=14,
                     taylor=True, beta=0.3, windows='joint'),
    # ---- materials (refractiveindex.info library, hops/materials.py) -----------------------------
    'Au_crossed': _s('gold crossed grating cos(x)cos(y), period 0.8 um: grating-coupled SPPs (dispersive)',
                     n_w='Au', period=0.8, q=(0, 1, 2), max_delta=0.05, omega_range=(0.5, 3.0)),
    'Ag_crossed': _s('silver crossed grating (cos x + cos y)/2, period 0.5 um', n_w='Ag', period=0.5,
                     profile='cosx+cosy', q=(0, 1, 2), max_delta=0.05, omega_range=(0.5, 3.0)),
    'Al_uv_crossed': _s('aluminium, period 0.25 um: UV plasmonics on a crossed grating', n_w='Al', period=0.25,
                        q=(0, 1, 2), max_delta=0.05, omega_range=(0.5, 3.0)),
    'water_over_gold_crossed': _s('gold crossed grating under water, period 0.6 um (SPR sensing)', n_u='water',
                                  n_w='Au', period=0.6, q=(0, 1, 2), max_delta=0.05, omega_range=(0.5, 3.0)),
    'TiO2_crossed': _s('rutile TiO2 crossed grating, period 1 um (lossless high index, joint windows)',
                       n_w='TiO2', period=1.0, q=(1, 2), windows='joint', max_delta=0.05),
    'Si_crossed': _s('silicon crossed grating, period 1 um', n_w='Si', period=1.0, q=(1, 2), windows='joint',
                     max_delta=0.05),
    'SiC_reststrahlen_crossed': _s('SiC crossed grating, period 10.5 um, q = 0: surface phonon polaritons',
                                   n_w='SiC', period=10.5, q=(0,), max_delta=0.02),
    'oblique_gold': _s('gold, cos(x)cos(y), fixed incidence theta = 20 deg, phi = 30 deg', n_w=1.48 + 1.883j,
                       theta=20.0, phi=30.0),
}


def _index(spec, lam, db=None):
    lit = mat.parse_index(spec)
    return lit if lit is not None else mat.refractive_index(spec, lam, db=db)


def _window(win, c):
    """Solve one frequency window (module level: runs in a worker process)."""
    key, omega_bar, dmax, edges = win
    n_u, n_w = c['n_u_w'][key], c['n_w_w'][key]
    if c['theta'] is not None:
        th, ph = np.deg2rad(c['theta']), np.deg2rad(c['phi'])
        k_u = np.real(n_u) * omega_bar
        alpha, beta = k_u * np.sin(th) * np.cos(ph), k_u * np.sin(th) * np.sin(ph)
    else:
        alpha, beta = c['alpha'], c['beta']
    N, M, Eps = c['N'], c['M'], c['Eps']
    delta = np.array([0.0]) if c['N_delta'] == 1 else np.linspace(-dmax, dmax, c['N_delta'])
    omega = omega_bar * (1 + delta)
    t0 = time.time()
    P = h3.make_problem(c['Nx'], c['Ny'], c['Nz'], N, M, n_u, n_w, omega_bar, alpha, beta, f=c['f'],
                        f_x=c['f_x'], f_y=c['f_y'], a=c['a'], b=c['b'], Mode=c['Mode'])
    zeta, psi = h3.setup_zeta_psi_n_m_3d(alpha, beta, P.gamma_u_bar, P.f, P.f_x, P.f_y, N, M)
    extra = {}
    if c.get('keep_fields'):
        # volume fields for pictures / movies (refl_movie_3D.py): keep the interface data U, W and one
        # vertical (x, z') slice of u_{n,m}, w_{n,m} at y = yy[iy]
        U, W, ubar, wbar, vu, vw = h3.two_layer_solve_3d_coupled(P, zeta, psi, keep_volume=True)
        iy = c.get('y_slice', 0) % c['Ny']
        extra = dict(U_n_m=U, W_n_m=W, u_xz=vu[:, iy].copy(), w_xz=vw[:, iy].copy(), iy=iy,
                     gamma_u_bar=P.gamma_u_bar, gamma_w_bar=P.gamma_w_bar, tau2=P.tau2, kx=P.kx, ky=P.ky,
                     dmax=dmax)
        del vu, vw
    else:
        U, W, ubar, wbar = h3.SOLVERS_3D[c['solver']](P, zeta, psi)
    t1 = time.time()
    kw = dict(sum_domain=c['sum_domain'])
    args = (P.tau2, ubar, wbar, P.kx, P.ky, alpha, beta, P.gamma_u_bar, P.gamma_w_bar, Eps, delta)
    ee_f, ru_f, rl_f = h3.energy_defect_3d(*args, 0, 0, 1, **kw)
    ee, ru, rl = h3.energy_defect_3d(*args, N, M, 1 if c['Taylor'] else 2,
                                     taylor_full_order=c.get('taylor_full_order', False), **kw)
    if c['verbose']:
        print(f'{key}: omega in [{omega.min():.3f}, {omega.max():.3f}] (window {edges[0]:.3f}-{edges[1]:.3f}), '
              f'n_u = {complex(n_u):.4g}, n_w = {complex(n_w):.4g}, alpha = {alpha:.3g}, beta = {beta:.3g}: '
              f'solve ({c["solver"]}) {t1 - t0:.1f} s, energy {time.time() - t1:.1f} s', flush=True)
    RR = ru / ru_f if c['relative'] else ru
    return dict(key=key, omega_bar=omega_bar, delta=delta, omega=omega, lam=2 * np.pi / omega, Eps=Eps,
                ee=ee, ru=ru, rl=rl, ru_flat=ru_f, RR=RR, n_u=n_u, n_w=n_w, alpha=alpha, beta=beta,
                ubar_n_m=ubar, wbar_n_m=wbar, edges=edges, **extra)


def run(scenario=SCENARIO, qq=None, N_Eps=None, N_delta=None, verbose=True, solver=None, workers=None, **over):
    """Compute the 3D reflectivity map for a scenario (+ keyword overrides of any scenario key:
    n_u, n_w, profile, eps_max, alpha, beta, theta, phi, M, N, Nx, Ny, Nz, a, b, mode, taylor, period,
    windows, max_delta, sigma, omega_range, relative, sum_domain, rii_db)."""
    solver = solver or SOLVER
    workers = WORKERS if workers is None else workers
    sc = dict(SCENARIOS[scenario])
    sc.update({k: v for k, v in over.items() if v is not None})
    qq = tuple(qq) if qq is not None else tuple(sc['q'])
    N_Eps = N_Eps or sc.get('n_eps', 100)
    N_delta = N_delta or sc.get('n_delta', 100)
    M = sc['M']
    N = sc.get('N') or M
    Nx, Ny, Nz = sc['Nx'], sc['Ny'], sc['Nz']
    period, db = sc.get('period'), sc.get('rii_db')
    named = [s for s in (sc['n_u'], sc['n_w']) if mat.parse_index(s) is None]
    if named and period is None:
        period = mat.MATERIALS.get(named[-1], {}).get('period', 1.0)
        print(f"note: named material(s) {named}: using grating period {period} um (set --period to change)")
    # frequency range from the bands q (as in 2D: band q = [q, q+1] / Re n^u for normal incidence)
    lam0 = period / (qq[0] + 0.5) if period else None
    n_u_ref = np.real(_index(sc['n_u'], lam0, db)) if named else np.real(mat.parse_index(sc['n_u']))
    if sc.get('omega_range'):
        o_lo, o_hi = sc['omega_range']
    else:
        o_lo, o_hi = min(qq) / n_u_ref, (max(qq) + 1) / n_u_ref
    n_w_ref = _index(sc['n_w'], period / (0.5 * (o_lo + o_hi)) if period else None, db)
    angle = None
    if sc.get('theta') is not None:
        th, ph = np.deg2rad(sc['theta']), np.deg2rad(sc.get('phi', 0.0))
        angle = (n_u_ref * np.sin(th) * np.cos(ph), n_u_ref * np.sin(th) * np.sin(ph))
    xx = (2 * np.pi / Nx) * np.arange(Nx)
    yy = (2 * np.pi / Ny) * np.arange(Ny)
    f, f_x, f_y = h3.profile_fn_3d(sc['profile'], xx, yy)
    lattice = None
    if sc.get('lattice', 'full') == 'profile':
        from hops3d.windows import profile_lattice
        nmax = max(abs(n_u_ref), abs(n_w_ref))
        lattice = profile_lattice(f, o_hi * abs(nmax) + abs(sc['alpha']) + abs(sc['beta']) + 2)
    wins = make_windows(o_lo, o_hi, n_u_ref, n_w_ref, sc['alpha'], sc['beta'], sc['sigma'], sc['windows'],
                        sc.get('min_window', 0.02), sc.get('max_delta'), Nx=Nx, Ny=Ny, angle=angle, lattice=lattice)
    if verbose:
        print(f'{len(wins)} frequency windows in omega = [{o_lo:.3f}, {o_hi:.3f}] '
              f'(cuts at the Rayleigh frequencies{" of both layers" if sc["windows"] == "joint" else ""})')
    per_win = lambda spec: {w[0]: _index(spec, period / w[1] if period else None, db) for w in wins}
    common = dict(n_u_w=per_win(sc['n_u']), n_w_w=per_win(sc['n_w']), theta=sc.get('theta'),
                  phi=sc.get('phi', 0.0), alpha=sc['alpha'], beta=sc['beta'], N=N, M=M, Nx=Nx, Ny=Ny, Nz=Nz,
                  f=f, f_x=f_x, f_y=f_y, a=sc['a'], b=sc['b'], Mode=1 if sc['mode'] == 'TE' else 2,
                  Eps=np.linspace(0, sc['eps_max'], N_Eps), N_delta=N_delta, Taylor=sc['taylor'], solver=solver,
                  verbose=verbose, relative=sc.get('relative', PlotRelative), sum_domain=sc['sum_domain'],
                  keep_fields=sc.get('keep_fields', False), y_slice=sc.get('y_slice', 0),
                  taylor_full_order=sc.get('taylor_full_order', False))
    for k, v in common['n_u_w'].items():
        if abs(np.imag(v)) > 1e-3:
            warnings.warn(f'window {k}: absorbing superstrate n_u = {v:.4g}')
    workers = min(workers or 1, len(wins))
    if workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(_window, wins, [common] * len(wins)))
    else:
        results = [_window(w, common) for w in wins]
    lossless = all(abs(np.imag(r['n_w'])) < 1e-6 and abs(np.imag(r['n_u'])) < 1e-6 for r in results)
    info = dict(scenario=scenario, N=N, M=M, Nx=Nx, Ny=Ny, Nz=Nz, n_u_spec=sc['n_u'], n_w_spec=sc['n_w'],
                Taylor=sc['taylor'], mode=sc['mode'], period=period, profile=sc['profile'], eps_max=sc['eps_max'],
                alpha=sc['alpha'], beta=sc['beta'], theta=sc.get('theta'), phi=sc.get('phi'), lossless=lossless,
                relative=common['relative'], desc=sc.get('desc', ''), windows=sc['windows'],
                a=sc['a'], b=sc['b'], f=f, f_x=f_x, f_y=f_y)
    return results, info


def plot(results, info, outdir=OUTDIR, tag=None):
    """Same rendering as refl_map.py: one contourf per window, shared colour scale ('hot')."""
    from matplotlib.colors import Normalize
    os.makedirs(outdir, exist_ok=True)
    tag = tag or info['scenario']
    period = info.get('period')
    if PlotLambda:
        xfun = (lambda r: r['lam'] * period / (2 * np.pi)) if period else (lambda r: r['lam'])
        xlabel = (r'$\lambda$ ($\mu$m)' + f'   (period {period:g} $\\mu$m)') if period else r'$\lambda$'
    else:
        xfun, xlabel = (lambda r: r['omega']), r'$\omega$'
    fmt = lambda v: (f'{mat.parse_index(v):.4g}' if mat.parse_index(v) is not None else str(v))
    inc = (f"$\\theta$={info['theta']:g}°, $\\phi$={info['phi']:g}°" if info.get('theta') is not None
           else f"$\\alpha$={info['alpha']:g}, $\\beta$={info['beta']:g}")
    sub = (f"  3D {info['profile']}, $n^u$={fmt(info['n_u_spec'])}, $n^w$={fmt(info['n_w_spec'])}, {inc}, "
           f"{'Taylor' if info.get('Taylor') else 'Pade'}")
    relative = info.get('relative', PlotRelative)
    figs = []
    for num, (title, getZ) in enumerate([('$D$', lambda r: safe_log10(r['ee'])),
                                         ('$R/R_{flat}$' if relative else '$R$', lambda r: np.real(r['RR']))],
                                        start=1):
        Zs = [getZ(r) for r in results]
        allz = np.concatenate([z[np.isfinite(z)].ravel() for z in Zs])
        vmin, vmax = allz.min(), allz.max()
        if num == 2 and np.mean(allz > 1.0 + 1e-9) < 1e-2 and vmax > 1.0:
            vmax = 1.0
            Zs = [np.minimum(z, 1.0) for z in Zs]
        norm = Normalize(vmin=vmin, vmax=vmax)
        fig, ax = plt.subplots(num=num, figsize=(7.5, 5), clear=True)
        for r, Z in zip(results, Zs):
            # D of a lossless map spans ~10 decades: one level per decade (as paper Fig. 9b).  For absorbing
            # layers D is the absorptance (< 1 decade of range): MATLAB-style automatic levels.
            matlab_contourf(ax, xfun(r), r['Eps'], Z, 'hot', norm, step=1.0 if (num == 1 and vmax - vmin >= 4) else None)
        shared_colorbar(fig, ax, 'hot', norm)
        ax.set_xlabel(xlabel, fontsize=13)
        ax.set_ylabel(r'$\varepsilon$', fontsize=16)
        note = '' if (num == 2 or info.get('lossless', True)) else '  (absorbing: D = absorption)'
        ax.set_title(title + sub + note, fontsize=9)
        fig.tight_layout()
        fname = os.path.join(outdir, f'refl_map_3D_{tag}_{"D" if num == 1 else "R"}.png')
        fig.savefig(fname, dpi=130)
        figs.append(fname)
    return figs


def save(results, info, tag, outdir=OUTDIR):
    os.makedirs(outdir, exist_ok=True)
    fn = os.path.join(outdir, f'refl_map_3D_{tag}.npz')
    np.savez_compressed(fn, period=np.nan if info['period'] is None else info['period'],
                        **{f'{r["key"]}_{k}': r[k] for r in results
                           for k in ('lam', 'omega', 'Eps', 'ee', 'ru', 'rl', 'RR')})
    return fn


def load_results(tag, outdir=OUTDIR):
    d = np.load(os.path.join(outdir, f'refl_map_3D_{tag}.npz'))
    keys = sorted({k.split('_')[0] for k in d.files if k.startswith('w')})
    res = [dict(key=kk, **{k: d[f'{kk}_{k}'] for k in ('lam', 'omega', 'Eps', 'ee', 'ru', 'rl', 'RR')}) for kk in keys]
    info = dict(scenario=tag, relative=True)
    if np.isfinite(d['period']):
        info['period'] = float(d['period'])
    return res, info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scenario', default=SCENARIO, choices=list(SCENARIOS))
    ap.add_argument('--list-scenarios', action='store_true')
    ap.add_argument('--list-materials', action='store_true')
    g = ap.add_argument_group('physics (override the scenario)')
    g.add_argument('--nu', help='upper index n^u: number or material key/path')
    g.add_argument('--nw', help='lower index n^w: number (1.5, 0.05+2.275i) or material key/path')
    g.add_argument('--mode', choices=['TE', 'TM'], help="tau^2 = 1 ('TE') or (n^u/n^w)^2 ('TM')")
    g.add_argument('--alpha', type=float, help='lateral Bloch wavenumber alpha (x)')
    g.add_argument('--beta', type=float, help='lateral Bloch wavenumber beta (y)')
    g.add_argument('--theta', type=float, help='polar incidence angle (deg): alpha, beta follow omega')
    g.add_argument('--phi', type=float, help='azimuthal incidence angle (deg), with --theta')
    g.add_argument('--period', type=float, help='grating period in micrometres')
    g.add_argument('--rii-db', help='path to refractiveindex.info database/data')
    g = ap.add_argument_group('geometry / numerics')
    g.add_argument('--profile', help=f'f(x, y): {", ".join(PROFILES_3D)}')
    g.add_argument('--eps-max', type=float)
    g.add_argument('--a', type=float)
    g.add_argument('--b', type=float)
    g.add_argument('--M', type=int, help='frequency order M (N = M unless --N)')
    g.add_argument('--N', type=int)
    g.add_argument('--Nx', type=int)
    g.add_argument('--Ny', type=int, help='Ny = 1: y-invariant grating (2D reduction, conical incidence)')
    g.add_argument('--Nz', type=int)
    g.add_argument('--summation', choices=['taylor', 'pade'])
    g.add_argument('--sum-domain', choices=['fourier', 'physical'])
    g.add_argument('--taylor-full-order', action='store_true',
                   help='Taylor: sum all min(N,M)+1 polar orders (MATLAB taylorsum_2_coeff sums only half)')
    g.add_argument('--windows', choices=['paper', 'joint'])
    g.add_argument('--max-delta', type=float, help='split windows so that |delta| <= this')
    g.add_argument('--lattice', choices=['full', 'profile'],
                   help="Rayleigh cuts at all orders (p,q) ('full', 2D-like) or only at orders the profile excites")
    g.add_argument('--sigma', type=float)
    g.add_argument('--q', type=int, nargs='*', help='bands: omega in [min q, max q + 1] / n^u')
    g.add_argument('--omega-range', type=float, nargs=2)
    g.add_argument('--neps', type=int)
    g.add_argument('--ndelta', type=int)
    g.add_argument('--absolute', action='store_true', help='plot R instead of R/R_flat')
    g = ap.add_argument_group('run')
    g.add_argument('--solver', default=SOLVER, choices=list(h3.SOLVERS_3D))
    g.add_argument('--workers', type=int, default=WORKERS)
    g.add_argument('--tag')
    g.add_argument('--no-show', action='store_true')
    args = ap.parse_args(argv)
    if args.list_scenarios:
        w = max(map(len, SCENARIOS))
        for k, v in SCENARIOS.items():
            print(f'{k:{w}s}  {v["desc"]}')
        return
    if args.list_materials:
        import refl_map
        refl_map.print_materials()
        return
    if args.no_show:
        matplotlib.use('Agg')
    over = dict(n_u=args.nu, n_w=args.nw, mode=args.mode, alpha=args.alpha, beta=args.beta, theta=args.theta,
                phi=args.phi, period=args.period, rii_db=args.rii_db, profile=args.profile, eps_max=args.eps_max,
                a=args.a, b=args.b, M=args.M, N=args.N, Nx=args.Nx, Ny=args.Ny, Nz=args.Nz,
                windows=args.windows, max_delta=args.max_delta, lattice=args.lattice, sigma=args.sigma, omega_range=args.omega_range,
                sum_domain=args.sum_domain, taylor=None if args.summation is None else args.summation == 'taylor',
                relative=False if args.absolute else None,
                taylor_full_order=True if args.taylor_full_order else None)
    t0 = time.time()
    res, info = run(args.scenario, args.q, args.neps, args.ndelta, solver=args.solver, workers=args.workers, **over)
    custom = any(v is not None for k, v in over.items() if k != 'relative')
    clean = lambda v: os.path.splitext(os.path.basename(str(v)))[0].replace('+', 'p').replace('.', 'd')
    tag = args.tag or (args.scenario if not custom else
                       f"{clean(info['n_u_spec'])}_over_{clean(info['n_w_spec'])}_{info['mode']}_{info['profile']}")
    for fn in plot(res, info, tag=tag):
        print('saved', fn)
    print('saved', save(res, info, tag))
    print(f'total {time.time() - t0:.1f} s')
    if SHOW and not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
