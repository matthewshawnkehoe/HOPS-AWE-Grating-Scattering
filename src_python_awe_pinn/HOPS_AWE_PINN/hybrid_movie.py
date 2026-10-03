"""hybrid_movie.py -- movies made with HOPS/AWE + PINN.

1. The refl_movie.py presets (silver_spp, silver_spp_eps, gold_spr, dielectric, sic_sphp, tir), rendered
   from the HYBRID: the reflectivity map, the energy defect / absorptance, the spectrum and the field are
   all computed from the physics-informed weights (refl_movie.py is reused unchanged; its refl_map.run
   and its field reconstruction are swapped for the hybrid ones).
2. New comparison movies (--compare): HOPS/AWE and HOPS/AWE + PINN side by side while the wavelength
   (or eps) sweeps through a band where AWE breaks down:
     wood_anomalies   n_w = 1.8 dielectric, two lower-layer Wood anomalies inside band q = 1
     large_eps        thesis Fig. 23 (n_u = 5, n_w = 8.1), eps 0 -> 0.4 at fixed lambda
     paper_fig9_edges paper Fig. 9 band q = 1, lambda sweep at eps = 0.2 into both band edges
   panels: spectrum R (AWE / hybrid), log10|D| or the indicator, the hybrid field, the AWE field.

    python hybrid_movie.py --preset dielectric
    python hybrid_movie.py --compare wood_anomalies
    python hybrid_movie.py --all
Output: figures/movies/*.mp4
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
from awepinn import core                                       # noqa: E402
from awepinn.core import rm, WindowBand, HybridSum              # noqa: E402
import refl_movie as rmv                                        # noqa: E402  (HOPS_Python, on sys.path via core)

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'figures', 'movies')
_CACHE = {}


def _hyb(r, info):
    key = (id(r), float(r['omega_bar']))
    if key not in _CACHE:
        _CACHE[key] = HybridSum(WindowBand(r, info), compress=core.COMPRESS)
    return _CACHE[key]


class HybridFieldRenderer(rmv.FieldRenderer):
    """refl_movie.FieldRenderer with the nodal fields from the physics-informed weights"""

    def nodal(self, r, eps, delta):
        S = _hyb(r, self.info)
        return S.fields(eps, float(r['omega_bar'] * (1 + delta)))

    def field(self, r, eps, delta):
        u, w = self.nodal(r, eps, delta)
        return self._assemble(r, eps, delta, u, w)

    def _assemble(self, r, eps, delta, u, w, incident=True):
        u = rmv._fourier_upsample(u, self.nx) @ self.C.T
        w = rmv._fourier_upsample(w, self.nx) @ self.C.T
        g = eps * self.f[:, None]
        Zu = g + self.zu[None, :] * (self.a - g) / self.a
        Zw = g + self.zw[None, :] * (self.b + g) / self.b
        alpha = self.info['alpha'] * (1 + delta)
        gam_u = (1 + delta) * r['gamma_u_bar']
        omega = r['omega_bar'] * (1 + delta)
        if self.ze_u.size:
            ue = self._extend(u[:, 0], r['n_u'] * omega, alpha, self.a, self.ze_u, up=True)
            u = np.concatenate([ue, u], axis=1)
            Zu = np.concatenate([np.broadcast_to(self.ze_u, (self.nx, self.ze_u.size)), Zu], axis=1)
        if self.ze_w.size:
            we = self._extend(w[:, -1], r['n_w'] * omega, alpha, -self.b, self.ze_w, up=False)
            w = np.concatenate([w, we], axis=1)
            Zw = np.concatenate([Zw, np.broadcast_to(self.ze_w, (self.nx, self.ze_w.size))], axis=1)
        phase = np.exp(1j * alpha * self.x)[:, None]
        Hu = phase * (u + (np.exp(-1j * gam_u * Zu) if incident else 0))
        Hw = phase * w
        d = 2 * np.pi
        X = np.concatenate([self.x + k * d for k in range(self.periods)] + [[self.periods * d]])[:, None]
        tile = lambda A: np.concatenate([A * np.exp(1j * alpha * d * k) for k in range(self.periods)]
                                        + [A[:1] * np.exp(1j * alpha * d * self.periods)], axis=0)
        tileZ = lambda A: np.concatenate([A] * self.periods + [A[:1]], axis=0)
        return X, tileZ(Zu), tile(Hu), tileZ(Zw), tile(Hw)


def _shim():
    def run(scenario, qq=None, N_Eps=None, N_delta=None, verbose=False, workers=None, **over):
        over.pop('keep_fields', None)
        res, info = core.run(scenario, qq=qq, N_Eps=N_Eps, N_delta=N_delta, workers=1, keep_fields=True,
                             verbose=True, **over)
        return res, info
    return types.SimpleNamespace(run=run, SCENARIOS=rm.SCENARIOS, OUTDIR=rm.OUTDIR)


def preset_movie(name, frames=120, n_grid=40):
    rmv.rm = _shim()
    rmv.FieldRenderer = HybridFieldRenderer
    rmv.OUTDIR = OUTDIR
    p = dict(rmv.PRESETS[name])
    sc = p.pop('scenario')
    q = p.pop('q')
    paths, anim, fig = rmv.make_movie(sc, q, frames=frames, n_grid=n_grid, out=f'hybrid_{name}', **p)
    plt.close(fig)
    return paths


COMPARE = {
    'wood_anomalies': dict(scenario='dielectric', q=(1,), over=dict(n_w=1.8), sweep='lambda', eps=0.15,
                           title='n_w = 1.8: two lower-layer Wood anomalies inside band q = 1'),
    'paper_fig9_edges': dict(scenario='dielectric', q=(1,), over={}, sweep='lambda', eps=0.2,
                             title='paper Fig. 9, band q = 1, eps = 0.2: the frequency-series truncation at the band edges'),
    'large_eps': dict(scenario='thesis23', q=(1,), over={}, sweep='eps', lam=None,
                      title='thesis Fig. 23 (n_u = 5, n_w = 8.1): eps 0 -> 0.4, Pade breaks down'),
}


def compare_movie(name, frames=120, n_grid=41, fps=15):
    t0 = time.time()
    c = COMPARE[name]
    res, info = core.run(c['scenario'], qq=c['q'], N_Eps=n_grid, N_delta=n_grid, workers=1, keep_fields=True,
                         verbose=True, **c['over'])
    r = res[0]
    S = _hyb(r, info)
    FR = HybridFieldRenderer(info, 1 if info['Taylor'] else 2)
    lossless = info['lossless']
    lam = r['lam']
    if c['sweep'] == 'lambda':
        e_idx = int(np.argmin(abs(r['Eps'] - c['eps'])))
        ks = np.linspace(0, len(lam) - 1, frames)
        dl = np.interp(ks, np.arange(len(lam)), r['delta'])
        es = np.full(frames, r['Eps'][e_idx])
    else:
        j0 = len(lam) // 2 if c.get('lam') is None else int(np.argmin(abs(lam - c['lam'])))
        dl = np.full(frames, r['delta'][j0])
        es = np.linspace(0, info['eps_max'], frames)
    fig = plt.figure(figsize=(14, 8.4))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.28)
    axS, axD, axI = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
    axFh, axFa = fig.add_subplot(gs[1, 0:2]), fig.add_subplot(gs[1, 2])
    # spectra (lambda sweep: R(lambda) at eps; eps sweep: R(eps) at lambda)
    if c['sweep'] == 'lambda':
        xs, xl = lam, r'$\lambda$'
        Ra, Rh = np.real(r['ru_awe'][e_idx]), np.real(r['ru'][e_idx])
        Da, Dh = np.real(r['ee_awe'][e_idx]), np.real(r['ee'][e_idx])
        I_ = r['indicator'][e_idx]
    else:
        xs, xl = r['Eps'], r'$\varepsilon$'
        Ra, Rh = np.real(r['ru_awe'][:, j0]), np.real(r['ru'][:, j0])
        Da, Dh = np.real(r['ee_awe'][:, j0]), np.real(r['ee'][:, j0])
        I_ = r['indicator'][:, j0]
    axS.plot(xs, Ra, 'C0-', lw=1.6, label='HOPS/AWE (refl_map.py)')
    axS.plot(xs, Rh, 'C3--', lw=1.6, label='HOPS/AWE + PINN')
    axS.set_xlabel(xl)
    axS.set_title('reflectivity R')
    axS.legend(fontsize=8)
    if lossless:
        axD.semilogy(xs, np.abs(Da) + 1e-17, 'C0-', label='AWE')
        axD.semilogy(xs, np.abs(Dh) + 1e-17, 'C3--', label='AWE + PINN')
        axD.set_title('energy defect |D| = |1 - R - T| (must be 0)')
    else:
        axD.plot(xs, Da, 'C0-', label='AWE')
        axD.plot(xs, Dh, 'C3--', label='AWE + PINN')
        axD.set_title('absorptance D')
    axD.set_xlabel(xl)
    axD.legend(fontsize=8)
    axI.semilogy(xs, I_ + 1e-30, 'm-')
    axI.axhline(core.DEFAULT_TOL, color='k', ls=':', lw=0.8)
    axI.set_xlabel(xl)
    axI.set_title('indicator: PINN loss of the AWE sum')
    marks = [ax.axvline(xs[0], color='c', lw=1) for ax in (axS, axD, axI)]
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-14, 0), cmap='viridis'), ax=axFa, fraction=0.04, pad=0.02)
    title = fig.suptitle('', fontsize=11)
    art = []
    sumtype = 1 if info['Taylor'] else 2

    def awe_nodal(eps, delta):
        from hops.summation import sum_series
        cu = np.transpose(r['u_n_m'], (0, 1, 3, 2))
        cw = np.transpose(r['w_n_m'], (0, 1, 3, 2))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return (sum_series(sumtype, cu, eps, delta, info['N'], info['M']),
                    sum_series(sumtype, cw, eps, delta, info['N'], info['M']))

    order = np.argsort(xs)

    def update(k):
        for a in art:
            a.remove()
        art.clear()
        e, d = float(es[k]), float(dl[k])
        u, w = FR.nodal(r, e, d)
        X, Zu, Hu, Zw, Hw = FR._assemble(r, e, d, u, w, incident=False)
        ua, wa = awe_nodal(e, d)
        _, _, Hua, _, Hwa = FR._assemble(r, e, d, ua, wa, incident=False)
        v = np.percentile(np.abs(np.real(np.concatenate([Hu.ravel(), Hw.ravel()]))), 99.5)
        Xs = np.broadcast_to(X, Zu.shape)
        art.append(axFh.pcolormesh(Xs, Zu, np.real(Hu), cmap='RdBu_r', vmin=-v, vmax=v, shading='gouraud'))
        art.append(axFh.pcolormesh(Xs, Zw, np.real(Hw), cmap='RdBu_r', vmin=-v, vmax=v, shading='gouraud'))
        Eu = np.log10(np.abs(Hu - Hua) + 1e-16)
        Ew = np.log10(np.abs(Hw - Hwa) + 1e-16)
        art.append(axFa.pcolormesh(Xs, Zu, Eu, cmap='viridis', vmin=-14, vmax=0, shading='gouraud'))
        art.append(axFa.pcolormesh(Xs, Zw, Ew, cmap='viridis', vmin=-14, vmax=0, shading='gouraud'))
        for ax in (axFh, axFa):
            art.append(ax.plot(Xs[:, -1], Zu[:, -1], 'k-', lw=1.2)[0])
            ax.set_xlim(Xs.min(), Xs.max())
            ax.set_ylim(-FR.H, FR.H)
            ax.set_aspect('equal')
        axFh.set_title('HOPS/AWE + PINN: Re of the scattered field (u above, w below)', fontsize=10)
        axFa.set_title(r'$\log_{10}$|field(AWE + PINN) - field(AWE)|', fontsize=10)
        x0 = 2 * np.pi / (r['omega_bar'] * (1 + d)) if c['sweep'] == 'lambda' else e
        for m in marks:
            m.set_xdata([x0, x0])
        s = S.point(e, float(r['omega_bar'] * (1 + d)), tol=core.DEFAULT_TOL)
        ra = np.interp(x0, xs[order], Ra[order])
        title.set_text(f"{c['title']}\n" + rf"$\lambda$ = {2 * np.pi / (r['omega_bar'] * (1 + d)):.4f},  "
                       rf"$\varepsilon$ = {e:.3f}:  R(AWE) = {ra:.7g},  R(AWE + PINN) = {s['R']:.7g},  "
                       f"indicator = {s['loss_awe_taylor']:.1e}")
        return art

    anim = animation.FuncAnimation(fig, update, frames=frames, interval=1000 / fps)
    os.makedirs(OUTDIR, exist_ok=True)
    p = os.path.join(OUTDIR, f'compare_{name}.mp4')
    anim.save(p, writer=animation.FFMpegWriter(fps=fps, bitrate=2400), dpi=85)
    plt.close(fig)
    print(f'saved {p} ({time.time() - t0:.0f} s)', flush=True)
    return p


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--preset', nargs='*', choices=list(rmv.PRESETS))
    ap.add_argument('--compare', nargs='*', choices=list(COMPARE))
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--frames', type=int, default=120)
    a = ap.parse_args()
    warnings.filterwarnings('ignore')
    presets = list(rmv.PRESETS) if a.all else (a.preset or [])
    compares = list(COMPARE) if a.all else (a.compare or [])
    for c_ in compares:
        compare_movie(c_, frames=a.frames)
    for p_ in presets:
        preset_movie(p_, frames=a.frames)
