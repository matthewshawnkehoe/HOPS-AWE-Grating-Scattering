"""refl_movie_3D.py -- movies of the 3D (crossed-grating) HOPS/AWE reflectivity map, energy defect and field.

3D analogue of refl_movie.py.  For a refl_map_3D.py scenario the HOPS/AWE coefficients are computed ONCE
per frequency window (with the volume fields kept on one vertical slice); every frame is then only a
Pade/Taylor summation.  Each frame shows

  top row     R(lambda, eps) (or R/R_flat) and the energy defect log10|D| / absorptance A = 1 - R - T,
              with a cursor at the current (lambda, eps); the spectrum R, T, A at the current eps
  bottom-left  Re u_total in the vertical plane y = y0 (x-z cut through two periods): incident + reflected
               above the grating, transmitted below; extended beyond z = a, -b by the exact Rayleigh series
  bottom-mid   |u| ON the crossed grating surface z = eps f(x, y) over 2 x 2 periods (top view) -- shows
               the 2D standing-wave / surface-plasmon hot-spot pattern that a 1D grating cannot produce
  bottom-right the diffraction orders in reciprocal space: every order (alpha + p, beta + q) as a dot; the
               light circle |k_par| = n^u omega (and n^w omega, dashed) moves outward as lambda decreases,
               orders enter at the Rayleigh anomalies; dot area ~ reflected efficiency R_pq (black) and
               transmitted efficiency T_pq (blue ring), area ~ sqrt(efficiency) so weak orders stay visible

  --sweep lambda (default) lambda runs through --lam-range at fixed eps = --eps
  --sweep eps    eps runs from 0 to eps_max at fixed lambda = --lam

Presets (--list):
  au_crossed_spp      gold cos(x)cos(y) grating, P = 0.8 um: the (+-1, +-1) orders turn evanescent at
                      0.566 um and couple to surface plasmons at 0.60 um (deep dip, hot spots at the crests)
  au_crossed_spp_eps  same at lambda = 0.60 um, eps grows 0 -> 0.2 (critical coupling, R -> small)
  silver_crossed      silver cos(4x)cos(4y) (paper-style nondimensional units): the (4, 4) plasmon at
                      lambda ~ 1.23 d/2pi
  egg_silver          silver "egg-crate" (cos x + cos y + cos x cos y)/3: several plasmon resonances
  dielectric_crossed  vacuum / n = 1.1, cos(x)cos(y), joint windows: Rayleigh/Wood anomalies as orders
                      (+-1, 0), (0, +-1), (+-1, +-1), ... enter the light circles
  conical             y-invariant cos(x) grating under CONICAL incidence (beta = 0.3): the orders shift off
                      the p-axis -- a 3D effect impossible in the 2D code
  oblique_gold        gold cos(x)cos(y) at theta = 20 deg, phi = 30 deg: asymmetric order pattern

Usage:  python refl_movie_3D.py --preset au_crossed_spp        (writes figures_3d/movies/movie3D_au_crossed_spp.mp4)
        python refl_movie_3D.py --scenario Ag_crossed --lam-range 0.35 0.6 --eps 0.1
Output: MP4 if ffmpeg is available (pip install imageio-ffmpeg), else GIF.
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

import refl_map_3D as rm
import hops3d as h3
from hops.summation import sum_series
from hops.plotting import safe_log10

OUTDIR = os.path.join(os.path.dirname(rm.OUTDIR), 'movies')

PRESETS = {
    'au_crossed_spp': dict(scenario='Au_crossed', q=(0, 1), lam_range=(0.50, 0.78), eps=0.12, M=10,
                           sweep='lambda', extra=dict(max_delta=0.01)),
    'au_crossed_spp_eps': dict(scenario='Au_crossed', q=(0, 1), lam=0.60, lam_range=(0.50, 0.78), M=10,
                               sweep='eps', extra=dict(max_delta=0.01)),
    'silver_crossed': dict(scenario='silver_lattice', q=(3, 4, 5, 6), lam_range=(1.0, 1.7), eps=0.12, M=12,
                           sweep='lambda', map_relative=True),
    'egg_silver': dict(scenario='egg_silver', q=(1, 2, 3), lam_range=None, eps=0.15, M=12, sweep='lambda',
                       map_relative=True),
    'dielectric_crossed': dict(scenario='dielectric_joint', q=(1, 2), lam_range=None, eps=0.15, M=12,
                               sweep='lambda', map_relative=True),
    'conical': dict(scenario='conical_1d', q=(1, 2), lam_range=None, eps=0.15, M=12, sweep='lambda',
                    map_relative=True),
    'oblique_gold': dict(scenario='oblique_gold', q=(1, 2), lam_range=None, eps=0.15, M=12, sweep='lambda',
                         map_relative=True),
}


# ----------------------------------------------------------------------------
# field reconstruction
# ----------------------------------------------------------------------------
def _cheb_interp_matrix(Nz, nfine):
    t = np.cos(np.pi * np.arange(Nz + 1) / Nz)
    tf = np.linspace(1, -1, nfine)
    V = np.polynomial.chebyshev.chebvander(t, Nz)
    Vf = np.polynomial.chebyshev.chebvander(tf, Nz)
    return Vf @ np.linalg.inv(V), tf


def _upsample(u, n_fine, axis=0):
    """Periodic spectral interpolation along one axis."""
    Nx = u.shape[axis]
    if Nx == 1:
        return np.repeat(u, n_fine, axis=axis)
    U = np.moveaxis(np.fft.fft(u, axis=axis), axis, 0)
    Uf = np.zeros((n_fine,) + U.shape[1:], dtype=complex)
    h = Nx // 2
    Uf[:h] = U[:h]
    Uf[-h:] = U[-h:]
    return np.moveaxis(np.fft.ifft(Uf, axis=0) * (n_fine / Nx), 0, axis)


def _outgoing(v):
    g = np.lib.scimath.sqrt(np.asarray(v, dtype=complex))
    return np.where(np.imag(g) < 0, -g, g)


class FieldRenderer3D:
    def __init__(self, info, sumtype, nx_fine=96, nz_fine=48, n_ext=40, height=None):
        self.info, self.st = info, sumtype
        self.N, self.M, self.Nx, self.Ny, self.Nz = info['N'], info['M'], info['Nx'], info['Ny'], info['Nz']
        self.a, self.b = info['a'], info['b']
        d = 2 * np.pi
        self.nx = nx_fine
        self.C, tf = _cheb_interp_matrix(self.Nz, nz_fine)
        self.zu = (self.a / 2) * (tf - 1) + self.a
        self.zw = (self.b / 2) * (tf - 1)
        self.x = d * np.arange(nx_fine) / nx_fine
        self.y = d * np.arange(nx_fine) / nx_fine
        self.H = max(height if height is not None else d / 2, self.a, self.b)
        self.ze_u = np.linspace(self.H, self.a, n_ext) if self.H > self.a else np.zeros(0)
        self.ze_w = np.linspace(-self.b, -self.H, n_ext) if self.H > self.b else np.zeros(0)
        f = info['f']
        self.f_fine = np.real(_upsample(_upsample(f.astype(complex), nx_fine, 0), nx_fine, 1))
        self.px = np.fft.fftfreq(nx_fine) * nx_fine                       # fine-grid wavenumbers (x)
        self.kx = h3.wavenumbers(self.Nx, d)
        self.ky = h3.wavenumbers(self.Ny, d)

    def _sum(self, c, eps, delta):
        """c (..., M+1, N+1) -> summed values (...)."""
        cc = np.swapaxes(c, -1, -2)                                         # (..., N+1, M+1)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return sum_series(self.st, cc, eps, delta, self.N, self.M)

    def xz_cut(self, r, eps, delta):
        """Total field in the plane y = yy[iy]: returns X (nx+1,1), Zu, Hu, Zw, Hw (two periods)."""
        iy = r['iy']
        y0 = 2 * np.pi * iy / self.Ny
        u = self._sum(r['u_xz'], eps, delta)                                # (Nx, Nz+1)
        w = self._sum(r['w_xz'], eps, delta)
        u = _upsample(u, self.nx, 0) @ self.C.T
        w = _upsample(w, self.nx, 0) @ self.C.T
        g = eps * np.real(_upsample(self.info['f'][:, iy:iy + 1].astype(complex), self.nx, 0))  # (nx, 1)
        Zu = g + self.zu[None, :] * (self.a - g) / self.a
        Zw = g + self.zw[None, :] * (self.b + g) / self.b
        s = 1 + delta
        alpha, beta = r['alpha'] * s, r['beta'] * s
        omega = r['omega_bar'] * s
        k_u, k_w = r['n_u'] * omega, r['n_w'] * omega
        gam_u = r['gamma_u_bar'] * s
        # beyond the artificial boundaries: exact Rayleigh expansion of the (x, y) traces, evaluated at y0
        for tr_key, zs, z0, k, up in (('ubar_n_m', self.ze_u, self.a, k_u, True),
                                      ('wbar_n_m', self.ze_w, -self.b, k_w, False)):
            if not zs.size:
                continue
            tr = self._sum(r[tr_key], eps, delta)                           # (Nx, Ny)
            ch = np.fft.fft2(tr) / (self.Nx * self.Ny)
            gp = _outgoing(k ** 2 - (alpha + self.kx[:, None]) ** 2 - (beta + self.ky[None, :]) ** 2)
            E = np.exp((1j if up else -1j) * gp[:, :, None] * (zs - z0)[None, None, :])    # (Nx, Ny, nz)
            cz = np.einsum('pqz,pq,q->pz', E, ch, np.exp(1j * self.ky * y0))               # (Nx, nz)
            full = np.zeros((self.nx, zs.size), dtype=complex)
            h = self.Nx // 2
            full[:h] = cz[:h]
            full[-h:] = cz[-h:]
            ext = np.fft.ifft(full, axis=0) * self.nx
            if up:
                u = np.concatenate([ext, u], axis=1)
                Zu = np.concatenate([np.broadcast_to(zs, (self.nx, zs.size)), Zu], axis=1)
            else:
                w = np.concatenate([w, ext], axis=1)
                Zw = np.concatenate([Zw, np.broadcast_to(zs, (self.nx, zs.size))], axis=1)
        phase = np.exp(1j * (alpha * self.x + beta * y0))[:, None]
        Hu = phase * (u + np.exp(-1j * gam_u * Zu))
        Hw = phase * w
        d = 2 * np.pi
        X = np.concatenate([self.x, self.x + d, [2 * d]])[:, None]
        tile = lambda A: np.concatenate([A, A * np.exp(1j * alpha * d), A[:1] * np.exp(2j * alpha * d)], axis=0)
        tileZ = lambda A: np.concatenate([A, A, A[:1]], axis=0)
        return X, tileZ(Zu), tile(Hu), tileZ(Zw), tile(Hw)

    def surface(self, r, eps, delta):
        """Total field on the grating surface z = eps f(x, y) (= W, the transmitted side), 2 x 2 periods."""
        s = 1 + delta
        alpha, beta = r['alpha'] * s, r['beta'] * s
        Wt = self._sum(r['W_n_m'], eps, delta)                               # (Nx, Ny)
        Wf = _upsample(_upsample(Wt, self.nx, 0), self.nx, 1)
        Wf = Wf * np.exp(1j * (alpha * self.x[:, None] + beta * self.y[None, :]))
        d = 2 * np.pi
        X = np.concatenate([self.x, self.x + d])
        Y = np.concatenate([self.y, self.y + d])
        big = np.block([[Wf, Wf * np.exp(1j * beta * d)], [Wf * np.exp(1j * alpha * d), Wf * np.exp(1j * (alpha + beta) * d)]])
        F = np.tile(self.f_fine, (2, 2))
        return X, Y, big, F

    def orders(self, r, eps, delta):
        """Per-order efficiencies at (eps, delta): (alpha_p, beta_q) grids, R_pq, T_pq, k_u, k_w."""
        _, _, _, Ru, Tw = h3.energy_defect_3d(r['tau2'], r['ubar_n_m'], r['wbar_n_m'], r['kx'], r['ky'],
                                              r['alpha'], r['beta'], r['gamma_u_bar'], r['gamma_w_bar'],
                                              np.array([eps]), np.array([delta]), self.N, self.M, self.st,
                                              return_modes=True)
        s = 1 + delta
        AP = r['alpha'] * s + r['kx'][:, None] + 0 * r['ky'][None, :]
        BQ = r['beta'] * s + r['ky'][None, :] + 0 * r['kx'][:, None]
        omega = r['omega_bar'] * s
        return AP, BQ, np.real(Ru[0, 0]), np.real(Tw[0, 0]), np.real(r['n_u']) * omega, np.real(r['n_w']) * omega


# ----------------------------------------------------------------------------
def _locate(results, omega):
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
        vals.append(np.real(r[key][j]))
    lam, vals = np.concatenate(lam), np.concatenate(vals)
    o = np.argsort(lam)
    return lam[o], vals[o]


def make_movie(scenario, q, lam_range=None, eps=0.1, lam=None, sweep='lambda', frames=120, fps=15, M=None,
               n_grid=50, out=None, workers=None, dpi=80, fmt='auto', verbose=True, extra=None,
               map_relative=False, quantity='re'):
    t0 = time.time()
    sc0 = rm.SCENARIOS[scenario]
    period = sc0.get('period')
    scale = (period / (2 * np.pi)) if period else 1.0
    over = dict(keep_fields=True, relative=False, solver='coupled')
    if M:
        over['M'] = M
    for k, v in (extra or {}).items():
        over[k] = v
    if lam_range is not None:
        lam_nd = np.array(lam_range) / scale
        over['omega_range'] = (2 * np.pi / lam_nd.max(), 2 * np.pi / lam_nd.min())
    solver = over.pop('solver')
    res, info = rm.run(scenario, qq=q, N_Eps=n_grid, N_delta=n_grid, verbose=False, workers=workers,
                       solver=solver, **over)
    res = sorted(res, key=lambda r: r['omega_bar'])
    if verbose:
        print(f'{scenario}: {len(res)} windows computed in {time.time() - t0:.0f} s', flush=True)
    st = 1 if info['Taylor'] else 2
    FR = FieldRenderer3D(info, st)
    lossless = info['lossless']
    lam_all = np.concatenate([r['lam'] for r in res]) * scale
    if lam_range is None:
        lam_range = (lam_all.min(), lam_all.max())
    lam_lo, lam_hi = lam_range
    unit = r' ($\mu$m)' if period else ''
    eps_max = info['eps_max']
    if sweep == 'lambda':
        lams = np.linspace(lam_hi, lam_lo, frames)          # decreasing lambda: orders enter the light circle
        epss = np.full(frames, eps)
    else:
        lams = np.full(frames, lam if lam is not None else 0.5 * (lam_lo + lam_hi))
        epss = np.linspace(0, eps_max, frames)

    # ------------------------------------------------------------------ static panels
    fig = plt.figure(figsize=(16, 9.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.25], hspace=0.30, wspace=0.28)
    axR, axD, axS = (fig.add_subplot(gs[0, i]) for i in range(3))
    axF, axT, axO = (fig.add_subplot(gs[1, i]) for i in range(3))
    shown = [r for r in res if (r['lam'] * scale).max() >= lam_lo and (r['lam'] * scale).min() <= lam_hi]
    Rkey = 'ru'
    if map_relative:
        for r in res:
            r['rrel'] = np.real(r['ru']) / np.real(r['ru_flat'])
        Rkey = 'rrel'
    Rall = np.concatenate([np.real(r[Rkey]).ravel() for r in shown])
    Rmax, Rmin = min(1.0, np.nanpercentile(Rall, 99.9)), max(0.0, np.nanpercentile(Rall, 0.1))
    if lossless:
        Dall = np.concatenate([safe_log10(r['ee']).ravel() for r in shown])
        Dall = Dall[np.isfinite(Dall)]
        dlim = (np.percentile(Dall, 1), np.percentile(Dall, 99.5))
    else:
        amax = max(1e-6, np.nanpercentile(np.concatenate([np.real(r['ee']).ravel() for r in shown]), 99.9))
    for r in shown:
        L = r['lam'] * scale
        axR.pcolormesh(L, r['Eps'], np.clip(np.real(r[Rkey]), Rmin, Rmax), cmap='hot', vmin=Rmin, vmax=Rmax,
                       shading='gouraud')
        if lossless:
            axD.pcolormesh(L, r['Eps'], safe_log10(r['ee']), cmap='hot', vmin=dlim[0], vmax=dlim[1], shading='gouraud')
        else:
            axD.pcolormesh(L, r['Eps'], np.clip(np.real(r['ee']), 0, 1), cmap='magma', vmin=0, vmax=amax,
                           shading='gouraud')
    for ax, t in ((axR, r'Reflectivity $R/R_{flat}$' if map_relative else 'Reflectivity $R$'),
                  (axD, r'Energy defect $\log_{10}|D|$' if lossless else r'Absorptance $A = 1 - R - T$')):
        ax.set_xlim(lam_lo, lam_hi)
        ax.set_ylim(0, eps_max)
        ax.set_xlabel(r'$\lambda$' + unit)
        ax.set_ylabel(r'$\varepsilon$')
        ax.set_title(t, fontsize=10)
    fig.colorbar(axR.collections[0], ax=axR, pad=0.01)
    fig.colorbar(axD.collections[0], ax=axD, pad=0.01)
    cursors = [ax.axvline(lams[0], color='c', lw=1.2) for ax in (axR, axD)]
    hcur = [ax.axhline(epss[0], color='c', lw=0.8, ls='--') for ax in (axR, axD)]
    dots = [ax.plot([lams[0]], [epss[0]], 'o', color='c', ms=6, mec='k')[0] for ax in (axR, axD)]

    if sweep == 'lambda':
        l, R = _spectrum(res, eps, 'ru')
        _, T = _spectrum(res, eps, 'rl')
        _, A = _spectrum(res, eps, 'ee')
        l = l * scale
        m = (l >= lam_lo) & (l <= lam_hi)
        weak = np.nanmax(R[m]) < 0.2              # weak reflector (dielectric): zoom on R, T on a twin axis
        axS.plot(l[m], R[m], 'k-', lw=1.6, label='R (reflected)')
        if not np.allclose(T[m], 0):
            axT2 = axS.twinx() if weak else axS
            axT2.plot(l[m], T[m], 'b-', lw=1.2, label='T (transmitted)')
            if weak:
                axT2.set_ylabel('T', color='b')
                axT2.tick_params(axis='y', colors='b')
        if not lossless:
            axS.plot(l[m], A[m], 'r-', lw=1.2, label='A (absorbed)')
        axS.set_xlim(lam_lo, lam_hi)
        axS.set_xlabel(r'$\lambda$' + unit)
        axS.set_title(rf'Spectrum at $\varepsilon$ = {eps:g}', fontsize=10)
        spec = (l, R)
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
        axS.set_title(rf'At $\lambda$ = {lams[0]:.4g}' + unit, fontsize=10)
    mark, = axS.plot([], [], 'o', color='c', mec='k', ms=8)
    if weak:
        axS.set_ylabel('R')
    else:
        axS.set_ylim(-0.02, 1.02)
    axS.grid(alpha=0.3)
    axS.legend(loc='best', fontsize=8)

    # dynamic panels
    qf = (lambda H: np.real(H)) if quantity == 're' else (lambda H: np.abs(H))
    cm_f, sym = ('RdBu_r', True) if quantity == 're' else ('inferno', False)
    smF = plt.cm.ScalarMappable(norm=plt.Normalize(-1, 1), cmap=cm_f)
    smT = plt.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap='inferno')
    fig.colorbar(smF, ax=axF, pad=0.01)
    fig.colorbar(smT, ax=axT, pad=0.01)
    iy0 = res[0].get('iy', 0)
    y0 = 2 * np.pi * iy0 / info['Ny'] * scale
    axF.set_title(('Re ' if quantity == 're' else '|') + r'$u_{tot}$' + ('' if quantity == 're' else '|') +
                  f' in the plane y = {y0:.3g}' + unit.replace('(', '').replace(')', ''), fontsize=10)
    axT.set_title(r'$|u_{tot}|$ on the grating surface $z = \varepsilon f(x,y)$ (top view)', fontsize=10)
    axO.set_title('Diffraction orders (reciprocal space)', fontsize=10)
    title = fig.suptitle('', fontsize=11)
    art = []
    ptxt = axF.text(0.01, 0.98, '', transform=axF.transAxes, va='top', fontsize=8, bbox=dict(fc='w', alpha=0.7, lw=0))

    def label(k, r):
        nu, nw = complex(r['n_u']), complex(r['n_w'])
        return (f"3D: {info.get('desc', scenario)[:110]}\n"
                rf"$\lambda$ = {lams[k]:.4g}{unit.replace('(', '').replace(')', '')},  $\varepsilon$ = {epss[k]:.3f},"
                rf"  $n^u$ = {nu.real:.3g}{nu.imag:+.3g}i,  $n^w$ = {nw.real:.3g}{nw.imag:+.3g}i,"
                rf"  $\alpha$ = {r['alpha']:.3g}, $\beta$ = {r['beta']:.3g}")

    th = np.linspace(0, 2 * np.pi, 200)

    def update(k):
        for a_ in art:
            a_.remove()
        art.clear()
        omega = 2 * np.pi / (lams[k] / scale)
        r, dl = _locate(res, omega)
        e = epss[k]
        # x-z cut
        X, Zu, Hu, Zw, Hw = FR.xz_cut(r, e, dl)
        Xs = np.broadcast_to(X * scale, Zu.shape)
        vmax = np.percentile(np.abs(np.concatenate([qf(Hu).ravel(), qf(Hw).ravel()])), 99.5)
        vmin = -vmax if sym else 0
        smF.set_clim(vmin, vmax)
        art.append(axF.pcolormesh(Xs, Zu * scale, qf(Hu), cmap=cm_f, vmin=vmin, vmax=vmax, shading='gouraud'))
        art.append(axF.pcolormesh(Xs, Zw * scale, qf(Hw), cmap=cm_f, vmin=vmin, vmax=vmax, shading='gouraud'))
        art.append(axF.plot(Xs[:, -1], Zu[:, -1] * scale, 'k-', lw=1.4)[0])
        axF.set_xlim(Xs.min(), Xs.max())
        axF.set_ylim(-FR.H * scale, FR.H * scale)
        axF.set_aspect('equal')
        axF.set_xlabel(r'$x$' + unit)
        axF.set_ylabel(r'$z$' + unit)
        hmax = max(np.abs(Hu).max(), np.abs(Hw).max())
        # surface top view
        Xg, Yg, S, F = FR.surface(r, e, dl)
        Sa = np.abs(S)
        smax, smin = np.percentile(Sa, 99.5), np.percentile(Sa, 0.5)
        if smax - smin < 1e-6 * smax:
            smin, smax = smin * 0.999, smax * 1.001
        smT.set_clim(smin, smax)
        art.append(axT.pcolormesh(Xg * scale, Yg * scale, Sa.T, cmap='inferno', vmin=smin, vmax=smax, shading='gouraud'))
        art.append(axT.contour(Xg * scale, Yg * scale, F.T, levels=[-0.5, 0.5], colors=['w', 'w'],
                               linewidths=0.6, linestyles=['--', '-']))
        axT.axhline(y0, color='c', lw=0.8, ls=':')
        axT.set_aspect('equal')
        axT.set_xlabel(r'$x$' + unit)
        axT.set_ylabel(r'$y$' + unit)
        # diffraction orders
        AP, BQ, Rpq, Tpq, ku, kw = FR.orders(r, e, dl)
        lim = 1.25 * max(ku, kw if lossless else ku, 1.0)
        sel = (np.abs(AP) < lim) & (np.abs(BQ) < lim)
        art.append(axO.scatter(AP[sel], BQ[sel], s=6, c='0.7', zorder=1))
        for vals, col in ((Rpq, 'k'), (Tpq, 'b')):     # marker area ~ sqrt(efficiency): weak orders stay visible
            mk = sel & (vals > 1e-5)
            if mk.any():
                art.append(axO.scatter(AP[mk], BQ[mk], s=700 * np.sqrt(vals[mk]), facecolors='none' if col == 'b' else col,
                                       edgecolors=col, lw=1.4, alpha=0.8, zorder=3))
        art.append(axO.plot(ku * np.cos(th), ku * np.sin(th), 'k-', lw=1)[0])
        if np.real(r['n_w']) > 0 and abs(np.imag(r['n_w'])) < 1e-3:
            art.append(axO.plot(kw * np.cos(th), kw * np.sin(th), 'b--', lw=1)[0])
        axO.set_xlim(-lim, lim)
        axO.set_ylim(-lim, lim)
        axO.set_aspect('equal')
        axO.set_xlabel(r'$\alpha + p$')
        axO.set_ylabel(r'$\beta + q$')
        lines = []
        for nm, vals in (('R', Rpq), ('T', Tpq)):
            big = [i for i in np.argsort(vals.ravel())[::-1][:3] if vals.ravel()[i] > 1e-4]
            if big:
                lines.append(f'{nm}: ' + ', '.join(
                    f'({int(round(r["kx"][i // vals.shape[1]]))},{int(round(r["ky"][i % vals.shape[1]]))}) '
                    f'{vals.ravel()[i]:.3f}' for i in big))
        art.append(axO.text(0.01, 0.01, '\n'.join(lines), transform=axO.transAxes, fontsize=7, va='bottom',
                            bbox=dict(fc='w', alpha=0.75, lw=0)))
        # cursors
        for c_ in cursors:
            c_.set_xdata([lams[k], lams[k]])
        for h_ in hcur:
            h_.set_ydata([e, e])
        for d_ in dots:
            d_.set_data([lams[k]], [e])
        if sweep == 'lambda':
            mark.set_data([lams[k]], [np.interp(lams[k], *spec)])
        else:
            j = np.argmin(abs(r['delta'] - dl))
            i = np.argmin(abs(r['Eps'] - e))
            mark.set_data([e], [np.real(r['ru'][i, j])])
        title.set_text(label(k, r))
        ptxt.set_text(f'max |u| / |incident| = {hmax:.2f}')
        return art

    anim = animation.FuncAnimation(fig, update, frames=frames, interval=1000 / fps, blit=False)
    tag = out or f'movie3D_{scenario}_{sweep}'
    os.makedirs(OUTDIR, exist_ok=True)
    use_mp4 = fmt in ('auto', 'mp4') and animation.writers.is_available('ffmpeg')
    if not use_mp4:
        try:
            import imageio_ffmpeg
            matplotlib.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
            use_mp4 = fmt in ('auto', 'mp4')
        except ImportError:
            pass
    paths = []
    if use_mp4:
        p = os.path.join(OUTDIR, tag + '.mp4')
        anim.save(p, writer=animation.FFMpegWriter(fps=fps, bitrate=2800), dpi=dpi)
        paths.append(p)
    if fmt in ('gif', 'both') or not use_mp4:
        p = os.path.join(OUTDIR, tag + '.gif')
        anim.save(p, writer=animation.PillowWriter(fps=fps), dpi=max(55, dpi - 20))
        paths.append(p)
    if verbose:
        print(f'saved {", ".join(paths)}  ({time.time() - t0:.0f} s)', flush=True)
    return paths, anim, fig


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--preset', choices=list(PRESETS))
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--scenario', help='any refl_map_3D.py scenario')
    ap.add_argument('--q', type=int, nargs='*')
    ap.add_argument('--lam-range', type=float, nargs=2)
    ap.add_argument('--eps', type=float, default=0.1)
    ap.add_argument('--lam', type=float)
    ap.add_argument('--sweep', choices=['lambda', 'eps'], default='lambda')
    ap.add_argument('--M', type=int)
    ap.add_argument('--frames', type=int, default=120)
    ap.add_argument('--fps', type=int, default=15)
    ap.add_argument('--quantity', choices=['re', 'abs'], default='re')
    ap.add_argument('--format', choices=['auto', 'mp4', 'gif', 'both'], default='auto')
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--show', action='store_true')
    a = ap.parse_args()
    if a.list:
        for k, v in PRESETS.items():
            print(f'{k:20s} scenario={v["scenario"]}, sweep={v["sweep"]}, lam_range={v.get("lam_range")}')
        raise SystemExit
    if not a.show:
        matplotlib.use('Agg')
    if a.preset:
        cfg = dict(PRESETS[a.preset])
        cfg['out'] = f'movie3D_{a.preset}'
    else:
        cfg = dict(scenario=a.scenario or 'Au_crossed', q=tuple(a.q) if a.q else None, lam_range=a.lam_range,
                   eps=a.eps, lam=a.lam, sweep=a.sweep, M=a.M)
    cfg.update(frames=a.frames, fps=a.fps, quantity=a.quantity, workers=a.workers, fmt=a.format)
    paths, anim, fig = make_movie(**cfg)
    if a.show:
        plt.show()
