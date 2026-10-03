"""refl_movie.py -- movies of the HOPS/AWE reflectivity map, energy defect and the 2D field.

For a scenario of refl_map.py the HOPS/AWE coefficients are computed ONCE (per frequency
window); every frame is then only a Pade/Taylor summation, so a 150-frame movie costs about as
much as a single reflectivity map.  Each frame shows

  top-left    R(lambda, eps) (or R/R_flat) with a cursor at the current (lambda, eps)
  top-right   energy defect log10|D| (lossless) or absorptance D = 1 - R - T (absorbing substrate)
  bottom-left the spectrum R, T, D at the current eps, with a moving marker
  bottom-right the total field Re H_y (TM) / Re E_y (TE) in physical space over two grating periods:
              incident + reflected field above the grating, transmitted field below it,
              reconstructed from the HOPS/AWE volume expansions u_{n,m}, w_{n,m}

  --sweep lambda   (default) lambda runs through --lam-range at fixed eps = --eps
  --sweep eps      eps runs from 0 to eps_max at fixed lambda = --lam

Movies (presets, --list):
  silver_spp     Ag grating (P = 0.5 um): Rayleigh anomaly at 0.5 um and the surface-plasmon dip at
                 0.525 um -- the field locks onto the metal surface at resonance
  silver_spp_eps same, fixed lambda = 0.525 um, eps grows 0 -> 0.2 (dip deepens, field localizes)
  gold_spr       gold grating under water (P = 0.6 um): SPR biosensor dip at 0.825 um
  dielectric     paper Fig. 9 / thesis Fig. 17 (n_w = 1.1): Wood anomalies, grazing diffraction orders
  sic_sphp       SiC in the mid-IR (P = 10.5 um): surface phonon polariton in the Reststrahlen band
  tir            glass over air at 53 deg (alpha = 6): total internal reflection, evanescent field

Usage:  python refl_movie.py --preset silver_spp            (writes figures/movies/movie_silver_spp.mp4/.gif)
        python refl_movie.py --scenario Na --lam-range 0.4 0.7 --eps 0.1 --q 0 1 --max-delta 0.03
Output: MP4 if ffmpeg is available (conda/pip 'imageio-ffmpeg' or system ffmpeg), else GIF.
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
from matplotlib import animation

import refl_map as rm
from hops.summation import sum_series
from hops.plotting import safe_log10

OUTDIR = os.path.join(os.path.dirname(rm.OUTDIR), 'movies')

PRESETS = {
    'silver_spp': dict(scenario='Ag_disp', q=(0, 1), lam_range=(0.40, 0.70), eps=0.12, max_delta=0.03,
                       M=12, sweep='lambda'),
    'silver_spp_eps': dict(scenario='Ag_disp', q=(0, 1), lam=0.5255, lam_range=(0.40, 0.70), max_delta=0.03,
                           M=12, sweep='eps'),
    'gold_spr': dict(scenario='water_over_gold', q=(0, 1), lam_range=(0.62, 1.10), eps=0.12, max_delta=0.03,
                     M=12, sweep='lambda'),
    'dielectric': dict(scenario='dielectric', q=(1, 2, 3, 4, 5, 6), lam_range=None, eps=0.15, M=16,
                       sweep='lambda', map_relative=True),
    'sic_sphp': dict(scenario='SiC_reststrahlen', q=(0,), lam_range=(10.6, 12.6), eps=0.12, max_delta=0.005,
                     M=12, sweep='lambda'),
    'tir': dict(scenario='glass_over_air', q=(1,), lam_range=None, eps=0.12, M=12, sweep='lambda'),
}


# ----------------------------------------------------------------------------
# field reconstruction
# ----------------------------------------------------------------------------
def _cheb_interp_matrix(Nz, nfine):
    """Matrix mapping values at t_j = cos(pi j/Nz) to values at nfine equispaced t in [-1, 1]."""
    t = np.cos(np.pi * np.arange(Nz + 1) / Nz)
    tf = np.linspace(1, -1, nfine)
    V = np.polynomial.chebyshev.chebvander(t, Nz)
    Vf = np.polynomial.chebyshev.chebvander(tf, Nz)
    return Vf @ np.linalg.inv(V), tf


def _fourier_upsample(u, nx_fine):
    """Periodic spectral interpolation along axis 0."""
    Nx = u.shape[0]
    U = np.fft.fft(u, axis=0)
    Uf = np.zeros((nx_fine,) + u.shape[1:], dtype=complex)
    h = Nx // 2
    Uf[:h] = U[:h]
    Uf[-h:] = U[-h:]
    return np.fft.ifft(Uf, axis=0) * (nx_fine / Nx)


class FieldRenderer:
    def __init__(self, info, sumtype, nx_fine=128, nz_fine=60, periods=2, height=None, n_ext=50):
        self.info, self.sumtype = info, sumtype
        self.N, self.M, self.Nx, self.Nz = info['N'], info['M'], info['Nx'], info['Nz']
        self.a, self.b = info['a'], info['b']
        self.nx, self.periods = nx_fine, periods
        self.C, tf = _cheb_interp_matrix(self.Nz, nz_fine)
        self.zu = (self.a / 2) * (tf - 1) + self.a        # z' in [0, a]  (top -> interface)
        self.zw = (self.b / 2) * (tf - 1)                 # z' in [-b, 0]
        d = 2 * np.pi
        self.x = d * np.arange(nx_fine) / nx_fine
        self.f = np.real(_fourier_upsample(info['f'][:, None].astype(complex), nx_fine))[:, 0]
        # the picture extends beyond the artificial boundaries z = a, -b by the exact outgoing
        # Rayleigh expansions (the same radiation conditions the DNO boundary conditions encode)
        self.H = max(height if height is not None else d / 2, self.a, self.b)
        self.ze_u = np.linspace(self.H, self.a, n_ext) if self.H > self.a else np.zeros(0)
        self.ze_w = np.linspace(-self.b, -self.H, n_ext) if self.H > self.b else np.zeros(0)
        self.p = np.fft.fftfreq(nx_fine) * nx_fine

    @staticmethod
    def _gamma(k, alpha_p):
        g = np.lib.scimath.sqrt(k ** 2 - alpha_p ** 2 + 0j)
        return np.where(np.imag(g) < 0, -g, g)                        # Im gamma >= 0 (outgoing)

    def _extend(self, trace, k, alpha, z0, zs, up):
        c = np.fft.fft(trace) / trace.size
        gp = self._gamma(k, alpha + self.p)
        E = np.exp((1j if up else -1j) * np.outer(gp, zs - z0))       # (nx, n_ext)
        return np.fft.ifft(c[:, None] * E * trace.size, axis=0)

    def field(self, r, eps, delta):
        """Total field on the physical (x, z) grid for window result r at (eps, delta)."""
        N, M = self.N, self.M
        cu = np.transpose(r['u_n_m'], (0, 1, 3, 2))       # (Nx, Nz+1, N+1, M+1)
        cw = np.transpose(r['w_n_m'], (0, 1, 3, 2))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            u = sum_series(self.sumtype, cu, eps, delta, N, M)      # (Nx, Nz+1) on Chebyshev z'
            w = sum_series(self.sumtype, cw, eps, delta, N, M)
        u = _fourier_upsample(u, self.nx) @ self.C.T                   # (nx, nz_fine)
        w = _fourier_upsample(w, self.nx) @ self.C.T
        g = eps * self.f[:, None]
        Zu = g + self.zu[None, :] * (self.a - g) / self.a              # physical z (upper)
        Zw = g + self.zw[None, :] * (self.b + g) / self.b              # physical z (lower)
        alpha = self.info['alpha'] * (1 + delta)
        gam_u = (1 + delta) * r['gamma_u_bar']
        omega = r['omega_bar'] * (1 + delta)
        if self.ze_u.size:                                             # above z = a
            ue = self._extend(u[:, 0], r['n_u'] * omega, alpha, self.a, self.ze_u, up=True)
            u = np.concatenate([ue, u], axis=1)
            Zu = np.concatenate([np.broadcast_to(self.ze_u, (self.nx, self.ze_u.size)), Zu], axis=1)
        if self.ze_w.size:                                             # below z = -b
            we = self._extend(w[:, -1], r['n_w'] * omega, alpha, -self.b, self.ze_w, up=False)
            w = np.concatenate([w, we], axis=1)
            Zw = np.concatenate([Zw, np.broadcast_to(self.ze_w, (self.nx, self.ze_w.size))], axis=1)
        phase = np.exp(1j * alpha * self.x)[:, None]
        Hu = phase * (u + np.exp(-1j * gam_u * Zu))                    # incident + reflected
        Hw = phase * w                                                  # transmitted
        # tile over several periods (Bloch phase e^{i alpha d} per period)
        d = 2 * np.pi
        X = np.concatenate([self.x + k * d for k in range(self.periods)] + [[self.periods * d]])[:, None]
        tile = lambda A: np.concatenate([A * np.exp(1j * alpha * d * k) for k in range(self.periods)]
                                        + [A[:1] * np.exp(1j * alpha * d * self.periods)], axis=0)
        tileZ = lambda A: np.concatenate([A] * self.periods + [A[:1]], axis=0)   # no Bloch phase on z!
        return X, tileZ(Zu), tile(Hu), tileZ(Zw), tile(Hw)


# ----------------------------------------------------------------------------
def _locate(results, omega):
    """Window result and delta for a frequency omega (nearest window if in a gap)."""
    best, bd = None, np.inf
    for r in results:
        lo, hi = r['omega'].min(), r['omega'].max()
        dist = 0 if lo <= omega <= hi else min(abs(omega - lo), abs(omega - hi))
        if dist < bd:
            best, bd = r, dist
    delta = np.clip(omega / best['omega_bar'] - 1, best['delta'].min(), best['delta'].max())
    return best, delta


def _spectrum(results, eps, key):
    lam, vals = [], []
    for r in results:
        j = np.argmin(abs(r['Eps'] - eps))
        lam.append(r['lam'])
        vals.append(np.real(r[key][j]) if key != 'ee' else r[key][j])
    lam, vals = np.concatenate(lam), np.concatenate(vals)
    o = np.argsort(lam)
    return lam[o], vals[o]


def make_movie(scenario, q, lam_range=None, eps=0.1, lam=None, sweep='lambda', frames=150, fps=20,
               M=None, max_delta=None, n_grid=60, quantity='re', out=None, workers=None, dpi=90,
               fmt='auto', verbose=True, extra=None, field_scale='frame', map_relative=False):
    t0 = time.time()
    info0 = rm.SCENARIOS[scenario]
    period = info0.get('period')
    scale = (period / (2 * np.pi)) if period else 1.0                  # nondim length -> um
    over = dict(keep_fields=True, relative=False)
    if M:
        over['M'] = M
    if max_delta:
        over['max_delta'] = max_delta
    if extra:
        over.update(extra)
    if lam_range is not None:
        lam_nd = np.array(lam_range) / scale                           # nondimensional lambda
        over['omega_range'] = (2 * np.pi / lam_nd.max(), 2 * np.pi / lam_nd.min())
    res, info = rm.run(scenario, qq=q, N_Eps=n_grid, N_delta=n_grid, verbose=False, workers=workers, **over)
    res = sorted(res, key=lambda r: r['omega_bar'])
    if verbose:
        print(f'{scenario}: {len(res)} windows computed in {time.time() - t0:.0f} s', flush=True)
    sumtype = 1 if info['Taylor'] else 2
    FR = FieldRenderer(info, sumtype)
    lossless = info['lossless']

    lam_all = np.concatenate([r['lam'] for r in res]) * scale
    if lam_range is None:
        lam_range = (lam_all.min(), lam_all.max())
    lam_lo, lam_hi = lam_range
    unit = r' ($\mu$m)' if period else ''
    eps_max = info['eps_max']

    # frame parameters
    if sweep == 'lambda':
        lams = np.linspace(lam_lo, lam_hi, frames)
        epss = np.full(frames, eps)
    else:
        lams = np.full(frames, lam if lam is not None else 0.5 * (lam_lo + lam_hi))
        epss = np.linspace(0, eps_max, frames)

    # ------------------------------------------------------------------ figure
    fig = plt.figure(figsize=(13, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.32, wspace=0.22)
    axR, axD = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axS, axF = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    def in_range(r):
        l = r['lam'] * scale
        return l.max() >= lam_lo and l.min() <= lam_hi

    shown = [r for r in res if in_range(r)]
    Rkey = 'ru'
    if map_relative:                       # R / R_flat, as in the paper's reflectivity maps
        for r in res:
            r['rrel'] = np.real(r['ru']) / np.real(r['ru_flat'])
        Rkey = 'rrel'
    Rall = np.concatenate([np.real(r[Rkey]).ravel() for r in shown])
    Rmax = min(1.0, np.nanpercentile(Rall, 99.9))
    Rmin = max(0.0, np.nanpercentile(Rall, 0.1))
    if lossless:
        Dall = np.concatenate([safe_log10(r['ee']).ravel() for r in shown])
        Dall = Dall[np.isfinite(Dall)]
        dlim = (np.percentile(Dall, 1), np.percentile(Dall, 99.5))
    else:
        Aall = np.concatenate([np.real(r['ee']).ravel() for r in shown])
        amax = max(1e-6, np.nanpercentile(Aall, 99.9))
    for r in shown:
        L = r['lam'] * scale
        axR.pcolormesh(L, r['Eps'], np.clip(np.real(r[Rkey]), Rmin, Rmax), cmap='hot', vmin=Rmin, vmax=Rmax,
                       shading='gouraud')
        if lossless:
            axD.pcolormesh(L, r['Eps'], safe_log10(r['ee']), cmap='hot', vmin=dlim[0], vmax=dlim[1],
                           shading='gouraud')
        else:
            axD.pcolormesh(L, r['Eps'], np.clip(np.real(r['ee']), 0, 1), cmap='magma', vmin=0,
                           vmax=amax, shading='gouraud')
    for ax, t in ((axR, r'Reflectivity $R/R_{flat}$' if map_relative else 'Reflectivity $R$'),
                  (axD, r'Energy defect $\log_{10}|D|$' if lossless else r'Absorptance $A = 1 - R - T$')):
        ax.set_xlim(lam_lo, lam_hi)
        ax.set_ylim(0, eps_max)
        ax.set_xlabel(r'$\lambda$' + unit)
        ax.set_ylabel(r'$\varepsilon$')
        ax.set_title(t, fontsize=11)
    fig.colorbar(axR.collections[0], ax=axR, pad=0.01)
    fig.colorbar(axD.collections[0], ax=axD, pad=0.01)
    cursors = [ax.axvline(lams[0], color='c', lw=1.2) for ax in (axR, axD)]
    hcur = [ax.axhline(epss[0], color='c', lw=0.8, ls='--') for ax in (axR, axD)]
    dots = [ax.plot([lams[0]], [epss[0]], 'o', color='c', ms=6, mec='k')[0] for ax in (axR, axD)]

    # spectrum / eps-profile panel
    def spectra(e):
        l, R = _spectrum(res, e, 'ru')
        _, T = _spectrum(res, e, 'rl')
        _, D = _spectrum(res, e, 'ee')
        return l * scale, R, T, D

    if sweep == 'lambda':
        l, R, T, D = spectra(eps)
        m = (l >= lam_lo) & (l <= lam_hi)
        weak = np.nanmax(R[m]) < 0.2              # weak reflector (dielectric): zoom on R
        lnR, = axS.plot(l[m], R[m], 'k-', lw=1.6, label='R (reflected)')
        if not np.allclose(T[m], 0):
            axT = axS.twinx() if weak else axS
            axT.plot(l[m], np.real(T[m]), 'b-', lw=1.2, label='T (transmitted)')
            if weak:
                axT.set_ylabel('T', color='b')
                axT.tick_params(axis='y', colors='b')
        if not lossless:
            axS.plot(l[m], np.real(D[m]), 'r-', lw=1.2, label='A (absorbed)')
        axS.set_xlim(lam_lo, lam_hi)
        axS.set_xlabel(r'$\lambda$' + unit)
        axS.set_title(rf'Spectrum at $\varepsilon$ = {eps:g}', fontsize=11)
        sx = l
        mark, = axS.plot([], [], 'o', color='c', mec='k', ms=8)
    else:
        weak = False
        rr, dl = _locate(res, 2 * np.pi / (lams[0] / scale))
        j = np.argmin(abs(rr['delta'] - dl))
        axS.plot(rr['Eps'], np.real(rr['ru'][:, j]), 'k-', lw=1.6, label='R (reflected)')
        if not np.allclose(rr['rl'][:, j], 0):
            axS.plot(rr['Eps'], np.real(rr['rl'][:, j]), 'b-', lw=1.2, label='T (transmitted)')
        if not lossless:
            axS.plot(rr['Eps'], np.real(rr['ee'][:, j]), 'r-', lw=1.2, label='A (absorbed)')
        axS.set_xlim(0, eps_max)
        axS.set_xlabel(r'$\varepsilon$')
        axS.set_title(rf'At $\lambda$ = {lams[0]:.4g}' + unit, fontsize=11)
        mark, = axS.plot([], [], 'o', color='c', mec='k', ms=8)
    if not (sweep == 'lambda' and weak):
        axS.set_ylim(-0.02, 1.02)
    else:
        axS.set_ylabel('R')
    axS.grid(alpha=0.3)
    axS.legend(loc='best', fontsize=8)

    # field panel: fix colour scale from a few sample frames
    def field_frame(k):
        omega = 2 * np.pi / (lams[k] / scale)
        r, dl = _locate(res, omega)
        X, Zu, Hu, Zw, Hw = FR.field(r, epss[k], dl)
        return X, Zu, Hu, Zw, Hw, r, dl

    samp = [field_frame(k) for k in np.linspace(0, frames - 1, 9).astype(int)]
    qfun = (lambda H: np.real(H)) if quantity == 're' else (lambda H: np.abs(H))
    vmax = np.percentile(np.concatenate([np.abs(qfun(s[2])).ravel() for s in samp] +
                                        [np.abs(qfun(s[4])).ravel() for s in samp]), 99)
    cmap, vmin = ('RdBu_r', -vmax) if quantity == 're' else ('inferno', 0)
    pol = 'H_y' if info.get('mode', 'TM') == 'TM' else 'E_y'
    axF.set_xlabel(r'$x$' + unit)
    axF.set_ylabel(r'$z$' + unit)
    axF.set_title((r'Re ' if quantity == 're' else '|') + f'${pol}$' + ('' if quantity == 're' else '|') +
                  '  (incident + reflected above, transmitted below)', fontsize=10)
    sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin, vmax), cmap=cmap)
    fig.colorbar(sm, ax=axF, pad=0.01)
    title = fig.suptitle('', fontsize=12)
    lims_fixed = (vmin, vmax)
    peak_txt = axF.text(0.01, 0.98, '', transform=axF.transAxes, va='top', fontsize=9,
                        bbox=dict(fc='w', alpha=0.7, lw=0))
    field_art = []

    def label(k, r):
        nu, nw = complex(r['n_u']), complex(r['n_w'])
        return (f"{info.get('desc', scenario)[:95]}\n"
                rf"$\lambda$ = {lams[k]:.4g}{unit.replace('(', '').replace(')', '')},  $\varepsilon$ = {epss[k]:.3f},"
                rf"  $n^u$ = {nu.real:.3g}{nu.imag:+.3g}i,  $n^w$ = {nw.real:.3g}{nw.imag:+.3g}i,  {info.get('mode', 'TM')}")

    def update(k):
        for a in field_art:
            a.remove()
        field_art.clear()
        X, Zu, Hu, Zw, Hw, r, dl = field_frame(k)
        Xs = np.broadcast_to(X * scale, Zu.shape)
        vmin, vmax = lims_fixed
        hmax = max(np.max(np.abs(Hu)), np.max(np.abs(Hw)))
        if field_scale == 'frame':            # rescale every frame; the peak |H| is printed
            vmax = np.percentile(np.abs(np.concatenate([qfun(Hu).ravel(), qfun(Hw).ravel()])), 99.5)
            vmin = -vmax if quantity == 're' else 0
        sm.set_clim(vmin, vmax)
        field_art.append(axF.pcolormesh(Xs, Zu * scale, qfun(Hu), cmap=cmap, vmin=vmin, vmax=vmax, shading='gouraud'))
        field_art.append(axF.pcolormesh(Xs, Zw * scale, qfun(Hw), cmap=cmap, vmin=vmin, vmax=vmax, shading='gouraud'))
        field_art.append(axF.plot(Xs[:, -1], Zu[:, -1] * scale, 'k-', lw=1.5)[0])      # the grating surface
        axF.set_xlim(Xs.min(), Xs.max())
        axF.set_ylim(-FR.H * scale, FR.H * scale)
        axF.set_aspect('equal')
        for c_ in cursors:
            c_.set_xdata([lams[k], lams[k]])
        for h_ in hcur:
            h_.set_ydata([epss[k], epss[k]])
        for d_ in dots:
            d_.set_data([lams[k]], [epss[k]])
        if sweep == 'lambda':
            l, R, T, D = spectra(epss[k])
            mark.set_data([lams[k]], [np.interp(lams[k], l, R)])
        else:
            j = np.argmin(abs(r['delta'] - dl))
            i = np.argmin(abs(r['Eps'] - epss[k]))
            mark.set_data([epss[k]], [np.real(r['ru'][i, j])])
        title.set_text(label(k, r))
        peak_txt.set_text(f'max |{pol}| / |incident| = {hmax:.2f}')
        return field_art

    anim = animation.FuncAnimation(fig, update, frames=frames, interval=1000 / fps, blit=False)
    tag = out or f'movie_{scenario}_{sweep}'
    os.makedirs(OUTDIR, exist_ok=True)
    use_mp4 = fmt in ('auto', 'mp4') and animation.writers.is_available('ffmpeg')
    if not use_mp4:
        try:                                  # ffmpeg shipped with imageio-ffmpeg (pip)
            import imageio_ffmpeg
            matplotlib.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
            use_mp4 = fmt in ('auto', 'mp4')
        except ImportError:
            pass
    paths = []
    if use_mp4:
        p = os.path.join(OUTDIR, tag + '.mp4')
        anim.save(p, writer=animation.FFMpegWriter(fps=fps, bitrate=2400), dpi=dpi)
        paths.append(p)
    if fmt in ('gif', 'both') or not use_mp4:
        p = os.path.join(OUTDIR, tag + '.gif')
        anim.save(p, writer=animation.PillowWriter(fps=fps), dpi=max(60, dpi - 25))
        paths.append(p)
    if verbose:
        print(f'saved {", ".join(paths)}  ({time.time() - t0:.0f} s)', flush=True)
    return paths, anim, fig


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--preset', choices=list(PRESETS))
    ap.add_argument('--list', action='store_true', help='list the movie presets')
    ap.add_argument('--scenario', help='any refl_map.py scenario (see refl_map.py --list-scenarios)')
    ap.add_argument('--q', type=int, nargs='*')
    ap.add_argument('--lam-range', type=float, nargs=2, help='wavelength range (um if the scenario has a period)')
    ap.add_argument('--eps', type=float, default=0.1, help='eps of the field picture (lambda sweep)')
    ap.add_argument('--lam', type=float, help='fixed wavelength for --sweep eps')
    ap.add_argument('--sweep', choices=['lambda', 'eps'], default='lambda')
    ap.add_argument('--M', type=int)
    ap.add_argument('--max-delta', type=float)
    ap.add_argument('--frames', type=int, default=150)
    ap.add_argument('--fps', type=int, default=20)
    ap.add_argument('--quantity', choices=['re', 'abs'], default='re')
    ap.add_argument('--field-scale', choices=['frame', 'fixed'], default='frame',
                    help="colour scale of the field panel: rescaled every frame (peak printed) or fixed")
    ap.add_argument('--format', choices=['auto', 'mp4', 'gif', 'both'], default='auto')
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--show', action='store_true', help='also open the animation window')
    a = ap.parse_args()
    if a.list:
        for k, v in PRESETS.items():
            print(f'{k:16s} scenario={v["scenario"]}, sweep={v["sweep"]}, lam_range={v.get("lam_range")}')
        raise SystemExit
    if not a.show:
        matplotlib.use('Agg')
    if a.preset:
        cfg = dict(PRESETS[a.preset])
        cfg['out'] = f'movie_{a.preset}'
    else:
        cfg = dict(scenario=a.scenario or 'Ag_disp', q=tuple(a.q) if a.q else None, lam_range=a.lam_range,
                   eps=a.eps, lam=a.lam, sweep=a.sweep, M=a.M, max_delta=a.max_delta)
    cfg.update(frames=a.frames, fps=a.fps, quantity=a.quantity, workers=a.workers, fmt=a.format,
               field_scale=a.field_scale)
    paths, anim, fig = make_movie(**cfg)
    if a.show:
        plt.show()
