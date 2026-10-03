"""hybrid_movie_3D.py -- movies with the 3D HOPS/AWE + PINN method.

1. The seven refl_movie_3D.py presets (au_crossed_spp, au_crossed_spp_eps, silver_crossed, egg_silver,
   dielectric_crossed, conical, oblique_gold), re-made with the hybrid: maps, spectra, diffraction orders and
   the field pictures (x-z cut, surface top view) all come from the physics-informed weights.
   refl_movie_3D.make_movie is reused unchanged; its refl_map_3D module and field renderer are swapped.
2. New comparison movies (AWE vs hybrid vs pointwise HOPS truth, frame by frame):
     compare_dielectric_eps     crossed dielectric cos x cos y, eps swept 0 -> 0.4 (twice the map range)
     compare_dielectric_lambda  joint windows, lambda swept through the Wood anomalies of both layers
     compare_TiO2_lambda        lossless high-index TiO2 crossed grating (period 1 um)
   Panels: R, energy defect |D| and the indicator along the sweep (AWE / hybrid / truth), the hybrid field on
   the crossed surface, log10 |u_hybrid - u_AWE| on the surface, and the hybrid x-z field cut.

    python hybrid_movie_3D.py                    # everything (~1 h)
    python hybrid_movie_3D.py --only conical compare_dielectric_eps
"""
import argparse
import os
import sys
import time
import types
import warnings

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from awepinn3d import core, reference as ref     # noqa: E402
from awepinn3d.core import rm3                   # noqa: E402
import refl_movie_3D as rmv                      # noqa: E402  (HOPS_Python, on sys.path via core)
from hops.plotting import safe_log10             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'figures', 'movies')
TOKENS = ('u_xz', 'w_xz', 'ubar_n_m', 'wbar_n_m', 'W_n_m')


# ---------------------------------------------------------------------------------------------------
class HybridFieldRenderer3D(rmv.FieldRenderer3D):
    """FieldRenderer3D whose fields come from the PI-sum weights of the window object r['W']."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self._cache, self._cur = {}, None

    def weights(self, r, eps, delta):
        key = (r['key'], round(float(eps), 12), round(float(delta), 12))
        if key not in self._cache:
            if len(self._cache) > 64:
                self._cache.clear()
            self._cache[key] = r['W'].solve(float(eps), float(r['omega_bar'] * (1 + delta)))['c']
        return self._cache[key]

    def _sum(self, c, eps, delta):
        if isinstance(c, str):
            r = self._cur
            W = r['W']
            cc = self.weights(r, eps, delta)
            cu, cw = cc[:W.nfeat['u']], cc[W.nfeat['u']:]
            iy = r.get('iy', 0)
            return {'u_xz': lambda: W.feat['u'][:, iy] @ cu, 'w_xz': lambda: W.feat['w'][:, iy] @ cw,
                    'ubar_n_m': lambda: W.feat['u'][:, :, 0] @ cu, 'wbar_n_m': lambda: W.feat['w'][:, :, -1] @ cw,
                    'W_n_m': lambda: W.feat['w'][:, :, 0] @ cw}[c]()
        return super()._sum(c, eps, delta)

    def xz_cut(self, r, eps, delta):
        self._cur = r
        return super().xz_cut(_tokenised(r), eps, delta)

    def surface(self, r, eps, delta):
        self._cur = r
        return super().surface(_tokenised(r), eps, delta)

    def orders(self, r, eps, delta):
        W = r['W']
        s = 1 + delta
        cc = self.weights(r, eps, delta)
        cu, cw = cc[:W.nfeat['u']], cc[W.nfeat['u']:]
        gu, pu = W._gamma_pq('u', s)
        gw, pw = W._gamma_pq('w', s)
        Ru = np.where(pu, np.real(gu / gu[0, 0]) * np.abs(W.uhat_top @ cu) ** 2, 0)
        Tw = np.where(pw, np.real(W.tau2 * gw / gu[0, 0]) * np.abs(W.what_bot @ cw) ** 2, 0)
        AP = r['alpha'] * s + r['kx'][:, None] + 0 * r['ky'][None, :]
        BQ = r['beta'] * s + r['ky'][None, :] + 0 * r['kx'][:, None]
        omega = r['omega_bar'] * s
        return AP, BQ, Ru, Tw, np.real(r['n_u']) * omega, np.real(r['n_w']) * omega


def _tokenised(r):
    d = dict(r)
    for k in TOKENS:
        d[k] = k
    return d


def _shim():
    """stand-in for the refl_map_3D module inside refl_movie_3D: run() returns the hybrid results"""
    m = types.SimpleNamespace(SCENARIOS=rm3.SCENARIOS, OUTDIR=rm3.OUTDIR, PlotRelative=rm3.PlotRelative)

    def run(scenario, qq=None, N_Eps=None, N_delta=None, verbose=True, solver=None, workers=None, **over):
        over = {k: v for k, v in over.items() if k not in ('keep_fields',)}
        res, info = core.run(scenario, qq=qq, N_Eps=N_Eps, N_delta=N_delta, keep_window='light', workers=1, **over)
        return res, info
    m.run = run
    return m


def preset_movie(name, n_grid=15, frames=120):
    p = dict(rmv.PRESETS[name])
    orig = (rmv.rm, rmv.FieldRenderer3D, rmv.OUTDIR)
    rmv.rm, rmv.FieldRenderer3D, rmv.OUTDIR = _shim(), HybridFieldRenderer3D, OUTDIR
    try:
        rmv.make_movie(p.pop('scenario'), p.pop('q'), frames=frames, n_grid=n_grid, out=f'hybrid_{name}', workers=1, **p)
    finally:
        rmv.rm, rmv.FieldRenderer3D, rmv.OUTDIR = orig


# ---------------------------------------------------------------------------------------------------
COMPARE = {
    'compare_dielectric_eps': dict(scenario='dielectric', q=(1,), sweep='eps', lam=None, eps_max=0.4, extra={}),
    'compare_dielectric_lambda': dict(scenario='dielectric_joint', q=(1, 2), sweep='lambda', eps=0.2, extra={}),
    'compare_TiO2_lambda': dict(scenario='TiO2_crossed', q=(1,), sweep='lambda', eps=0.15, extra={}),
}


def compare_movie(name, frames=100, fps=12, n_grid=11, dpi=80):
    t0 = time.time()
    p = COMPARE[name]
    sc0 = rm3.SCENARIOS[p['scenario']]
    over = dict(p['extra'])
    if p.get('eps_max'):
        over['eps_max'] = p['eps_max']
    res, info = core.run(p['scenario'], qq=p['q'], N_Eps=n_grid, N_delta=n_grid, keep_window=True, workers=1, **over)
    res = sorted(res, key=lambda r: r['omega_bar'])
    period = info.get('period')
    scale = period / (2 * np.pi) if period else 1.0
    unit = r' ($\mu$m)' if period else ''
    st = 1 if info['Taylor'] else 2
    FRa = rmv.FieldRenderer3D(info, st)               # AWE (refl_map_3D summation)
    FRh = HybridFieldRenderer3D(info, st)
    eps_max = info['eps_max']
    lam_all = np.concatenate([r['lam'] for r in res]) * scale
    if p['sweep'] == 'lambda':
        lams = np.linspace(lam_all.max() * 0.999, lam_all.min() * 1.001, frames)
        epss = np.full(frames, p['eps'])
    else:
        r0 = res[len(res) // 2]
        lams = np.full(frames, r0['lam'][len(r0['lam']) // 3] * scale)
        epss = np.linspace(0, eps_max, frames)
    # sweep values: AWE, hybrid, indicator, truth (pointwise HOPS on the same grid, N = 24)
    cur = []
    for k in range(frames):
        r, dl = rmv._locate(res, 2 * np.pi / (lams[k] / scale))
        cur.append((r, float(dl)))
    vals = {q: np.zeros(frames) for q in ('Ra', 'Rh', 'Rt', 'Da', 'Dh', 'Dt', 'ind')}
    awe_args = {}
    tinfo = dict(profile=info['profile'], mode=info['mode'], a=info['a'], b=info['b'], Nx=info['Nx'], Ny=info['Ny'],
                 Nz=info['Nz'])
    truth_cache = {}
    for k, (r, dl) in enumerate(cur):
        e, W = float(epss[k]), r['W']
        om = r['omega_bar'] * (1 + dl)
        ee, ru, rl = rm3.h3.energy_defect_3d(r['tau2'], r['ubar_n_m'], r['wbar_n_m'], r['kx'], r['ky'], r['alpha'],
                                             r['beta'], r['gamma_u_bar'], r['gamma_w_bar'], np.array([e]),
                                             np.array([dl]), info['N'], info['M'], st)
        s = W.point(e, om)
        vals['Ra'][k], vals['Da'][k] = np.real(ru[0, 0]), np.real(ee[0, 0])
        vals['Rh'][k], vals['Dh'][k], vals['ind'][k] = s['R'], s['D'], s['loss_awe_taylor']
        tk = (r['key'], round(dl, 12))
        if tk not in truth_cache:
            sv = 1 + dl
            grid_e = epss if p['sweep'] == 'eps' else np.array([e])
            truth_cache[tk] = ref.hops_pointwise(tinfo, r['n_u'], r['n_w'], om, r['alpha'] * sv, r['beta'] * sv, grid_e)
        Rt, Tt, Dt = truth_cache[tk]
        idx = k if p['sweep'] == 'eps' else 0
        vals['Rt'][k], vals['Dt'][k] = Rt[idx], Dt[idx]
    xs = epss if p['sweep'] == 'eps' else lams
    xlabel = r'$\varepsilon$' if p['sweep'] == 'eps' else r'$\lambda$' + unit
    o = np.argsort(xs)
    # ---- figure
    fig = plt.figure(figsize=(16, 9.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.2], hspace=0.32, wspace=0.28)
    aR, aD, aI = (fig.add_subplot(gs[0, i]) for i in range(3))
    aS, aE, aF = (fig.add_subplot(gs[1, i]) for i in range(3))
    aR.plot(xs[o], vals['Rt'][o], '-', color='0.6', lw=4, label='pointwise HOPS (truth, same grid)')
    aR.plot(xs[o], vals['Ra'][o], 'C0-', lw=1.4, label='3D HOPS/AWE')
    aR.plot(xs[o], vals['Rh'][o], 'C3--', lw=1.4, label='3D HOPS/AWE + PINN')
    aR.set_title('Reflectivity R along the sweep', fontsize=10)
    aR.legend(fontsize=8)
    lossless = info['lossless']
    Lg = lambda v: np.log10(np.abs(v) + 1e-17)
    if lossless:
        aD.plot(xs[o], Lg(vals['Da'][o]), 'C0-', lw=1.4, label='AWE')
        aD.plot(xs[o], Lg(vals['Dh'][o]), 'C3-', lw=1.4, label='AWE + PINN')
        aD.set_title(r'Energy defect $\log_{10}|D|$ (lossless: exact D = 0)', fontsize=10)
    else:
        aD.plot(xs[o], Lg(vals['Ra'][o] - vals['Rt'][o]), 'C0-', label='AWE')
        aD.plot(xs[o], Lg(vals['Rh'][o] - vals['Rt'][o]), 'C3-', label='AWE + PINN')
        aD.set_title(r'$\log_{10}|R - R_{true}|$', fontsize=10)
    aD.legend(fontsize=8)
    aI.semilogy(xs[o], np.abs(vals['Ra'][o] - vals['Rt'][o]) + 1e-17, 'C0-', lw=1.2, label='|R_AWE - R_true|')
    aI.semilogy(xs[o], np.abs(vals['Rh'][o] - vals['Rt'][o]) + 1e-17, 'C3-', lw=1.2, label='|R_hybrid - R_true|')
    aI.semilogy(xs[o], vals['ind'][o] + 1e-30, 'k:', lw=1.2, label='indicator (PINN loss of the AWE sum)')
    aI.set_title('error vs the reference, and the a-posteriori indicator', fontsize=10)
    aI.legend(fontsize=7)
    for ax in (aR, aD, aI):
        ax.set_xlabel(xlabel)
        ax.grid(alpha=0.3)
    curs = [ax.axvline(xs[0], color='c', lw=1.2) for ax in (aR, aD, aI)]
    smS = plt.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap='inferno')
    smE = plt.cm.ScalarMappable(norm=plt.Normalize(-16, 0), cmap='viridis')
    smF = plt.cm.ScalarMappable(norm=plt.Normalize(-1, 1), cmap='RdBu_r')
    fig.colorbar(smS, ax=aS, pad=0.01); fig.colorbar(smE, ax=aE, pad=0.01); fig.colorbar(smF, ax=aF, pad=0.01)
    aS.set_title(r'hybrid: $|u_{tot}|$ on the surface $z = \varepsilon f(x, y)$', fontsize=10)
    aE.set_title(r'$\log_{10}|u_{hybrid} - u_{AWE}|$ on the surface', fontsize=10)
    aF.set_title(r'hybrid: Re $u_{tot}$ in the plane y = 0', fontsize=10)
    title = fig.suptitle('', fontsize=11)
    art = []

    def update(k):
        for a_ in art:
            a_.remove()
        art.clear()
        r, dl = cur[k]
        e = float(epss[k])
        Xg, Yg, Sh, F = FRh.surface(r, e, dl)
        _, _, Sa, _ = FRa.surface(r, e, dl)
        A = np.abs(Sh)
        smS.set_clim(np.percentile(A, 0.5), np.percentile(A, 99.5) + 1e-12)
        art.append(aS.pcolormesh(Xg * scale, Yg * scale, A.T, cmap='inferno', vmin=smS.norm.vmin, vmax=smS.norm.vmax,
                                 shading='gouraud'))
        Ed = np.log10(np.abs(Sh - Sa) + 1e-17)
        art.append(aE.pcolormesh(Xg * scale, Yg * scale, Ed.T, cmap='viridis', vmin=-16, vmax=0, shading='gouraud'))
        X, Zu, Hu, Zw, Hw = FRh.xz_cut(r, e, dl)
        Xs = np.broadcast_to(X * scale, Zu.shape)
        vmax = np.percentile(np.abs(np.concatenate([np.real(Hu).ravel(), np.real(Hw).ravel()])), 99.5)
        smF.set_clim(-vmax, vmax)
        art.append(aF.pcolormesh(Xs, Zu * scale, np.real(Hu), cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='gouraud'))
        art.append(aF.pcolormesh(Xs, Zw * scale, np.real(Hw), cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='gouraud'))
        art.append(aF.plot(Xs[:, -1], Zu[:, -1] * scale, 'k-', lw=1.2)[0])
        for ax in (aS, aE):
            ax.set_aspect('equal'); ax.set_xlabel('x' + unit); ax.set_ylabel('y' + unit)
        aF.set_aspect('equal'); aF.set_xlabel('x' + unit); aF.set_ylabel('z' + unit)
        aF.set_ylim(-FRh.H * scale, FRh.H * scale)
        for c_ in curs:
            c_.set_xdata([xs[k], xs[k]])
        title.set_text(f"3D HOPS/AWE vs HOPS/AWE + PINN: {info.get('desc', '')[:100]}\n"
                       rf"$\lambda$ = {lams[k]:.4g}{unit},  $\varepsilon$ = {e:.3f}:   R AWE {vals['Ra'][k]:.6f}, "
                       f"hybrid {vals['Rh'][k]:.6f}, truth {vals['Rt'][k]:.6f};  indicator {vals['ind'][k]:.1e}")
        return art

    anim = animation.FuncAnimation(fig, update, frames=frames, interval=1000 / fps, blit=False)
    os.makedirs(OUTDIR, exist_ok=True)
    if not animation.writers.is_available('ffmpeg'):
        import imageio_ffmpeg
        matplotlib.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
    path = os.path.join(OUTDIR, name + '.mp4')
    anim.save(path, writer=animation.FFMpegWriter(fps=fps, bitrate=2800), dpi=dpi)
    plt.close(fig)
    print(f'saved {path} ({time.time() - t0:.0f} s)', flush=True)
    return path


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--frames', type=int, default=120)
    ap.add_argument('--n-grid', type=int, default=15)
    a = ap.parse_args()
    warnings.filterwarnings('ignore')
    names = a.only or (list(rmv.PRESETS) + list(COMPARE))
    for n in names:
        if n in COMPARE:
            compare_movie(n)
        else:
            preset_movie(n, n_grid=a.n_grid, frames=a.frames)
