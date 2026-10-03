"""refl_map.py -- Python port of refl_map.m (Reflectivity Map R and energy defect D),
extended with a materials library (refractiveindex.info), TE/TM, arbitrary superstrates,
the thesis examples (Kehoe PhD thesis, Ch. 6, Figs. 17-33) and many more profiles.

Quick start (PyCharm: set the script parameters in the run configuration)
  python refl_map.py                                   # refl_map.m as shipped (silver, Fig. 10a)
  python refl_map.py --scenario gold                   # paper Fig. 10b
  python refl_map.py --list-scenarios                  # every preset (paper, thesis, new materials)
  python refl_map.py --list-materials                  # the curated refractiveindex.info library
  python refl_map.py --nw Al --period 0.5              # any material below vacuum, grating period 0.5 um
  python refl_map.py --nu water --nw Au --period 1.0   # gold grating under water (SPR sensing)
  python refl_map.py --nu 1.5 --nw 1.0                 # glass over air: total internal reflection
  python refl_map.py --nw 3.8313+2.9043i --profile sin4x --M 15      # literal index (thesis Fig. 20a)
  python refl_map.py --nw main/Au/nk/Olmon-sc.yml --rii-db ~/rii/database/data --period 10

Physical units.  The computation is non-dimensional: period d = 2*pi, c0 = 1, and the six
bands q = 1..6 cover omega in [1, 7] (lambda = 2*pi/omega in [0.9, 6.3]).  For a grating of
physical period P (micrometres, --period) the vacuum wavelength is  lambda_phys = P / omega,
so band q is centred on lambda_q = P / (q + 1/2) (P = 1 um: q=1 -> 0.67 um ... q=6 -> 0.15 um).
Materials given by name are evaluated at these wavelengths:
  --index-mode band   (default) n evaluated at the centre of each band (piecewise-constant
                      dispersion: inside a band n is held fixed, which is what the HOPS/AWE
                      frequency expansion assumes);
  --index-mode fixed  one n for all bands, at --lambda-ref (the thesis' "representative value").
Literal indices (e.g. --nw 0.05+2.275i) are used as given.

Band centres.  The MATLAB code centres band q at omega_q = q + 1/2, which is the midpoint
between Rayleigh singularities only for n^u = 1 and alpha = 0.  --band-center rayleigh
(default 'auto' = rayleigh) uses paper eqs. (28)-(30):
    omega_q = (c0/n^u)(alpha + q + 1/2),   |delta| < sigma / (2 alpha + 2q + 1)   (d = 2 pi),
identical to MATLAB for n^u = 1, alpha = 0.
"""
import argparse
import os
import time
import warnings

# One BLAS thread per process: the bands run in parallel processes, and the matrices
# (33 x 33) are too small to benefit from threaded BLAS (avoids oversubscription).
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from hops import (cheb, setup_2d, csqrt, setup_zeta_psi_n_m, two_layer_solve_fast,
                  energy_defect, fourier_repr_lipschitz, fourier_repr_rough)
from hops.plotting import safe_log10
from hops.operators import two_layer_solve_operator, two_layer_solve_auto, two_layer_solve_lean
from hops.coupled import two_layer_solve_coupled
from hops import materials as mat

SOLVERS = {'coupled': two_layer_solve_coupled, 'operator': two_layer_solve_operator, 'fast': two_layer_solve_fast,
           'auto': two_layer_solve_auto, 'lean': two_layer_solve_lean}

# ----------------------------------------------------------------------------
# User settings (mirror the top of refl_map.m)
# ----------------------------------------------------------------------------
SCENARIO = 'silver'
PlotLambda = 1
PlotRelative = 1
RunNumber = 100
Mode = 2                  # 1 = TE, 2 = TM
N_delta = 100
N_Eps = 100
SHOW = True               # plt.show() at the end (set False for batch runs)
SOLVER = 'coupled'        # 'coupled' (fastest) | 'operator' | 'fast' (= two_layer_solve_fast.m) | 'auto' | 'lean'
WORKERS = os.cpu_count()  # frequency bands solved in parallel processes (1 = serial)
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'refl_map')

# ----------------------------------------------------------------------------
# Scenarios.  Keys: n_u, n_w (literal or material), profile, eps_max, alpha, M (=N unless N
# given), Nx, Nz, a, b, mode ('TM'/'TE'), taylor (True = Taylor, False = Pade), q (bands),
# period (um, for named materials), band_center, n_eps/n_delta (grid), desc.
# ----------------------------------------------------------------------------
_P = dict(n_u=1.0, profile='cosx', eps_max=0.2, alpha=0.0, M=15, Nx=32, Nz=32, a=1.0, b=1.0,
          mode='TM', taylor=False, q=(1, 2, 3, 4, 5, 6), period=None, band_center='auto')


def _s(desc, **kw):
    d = dict(_P)
    d.update(kw)
    d['desc'] = desc
    return d


SCENARIOS = {
    # ---- paper (J. Sci. Comput. 2024) / MATLAB --------------------------------------------
    'silver': _s('refl_map.m in src.zip, paper Fig. 10a: vacuum / silver, f = cos(4x), Pade',
                 n_w=0.05 + 2.275j, profile='cos4x'),
    'gold': _s('paper Fig. 10b: vacuum / gold (Johnson & Christy), f = cos(4x), Pade',
               n_w=1.48 + 1.883j, profile='cos4x'),
    'dielectric': _s('test_scenarios.m, GitHub README, paper Fig. 9 / thesis Fig. 17: vacuum / n=1.1, Taylor',
                     n_w=1.1, M=16, taylor=True),
    'dielectric_alpha': _s('paper Fig. 14: as dielectric with alpha = 0.01 (paper: 1000 x 1000 grid)',
                           n_w=1.1, M=16, taylor=True, alpha=0.01, band_center='matlab'),
    # ---- thesis, Chapter 6 (TM: 6.3-6.4, TE: 6.5) ------------------------------------------
    'thesis17': _s('thesis Fig. 17: (6.4)-(6.5), n_w = 1.1, cos(x), Taylor, N=M=16', n_w=1.1, M=16, taylor=True),
    'thesis18': _s('thesis Fig. 18: alpha = 1e-4, N_eps = N_delta = 1000, Taylor', n_w=1.1, M=16, taylor=True,
                   alpha=1e-4, n_eps=1000, n_delta=1000, band_center='matlab'),
    'thesis19a': _s('thesis Fig. 19a: silver (J&C), cos(4x), Pade, N=M=15', n_w=0.05 + 2.275j, profile='cos4x'),
    'thesis19b': _s('thesis Fig. 19b: gold (J&C), cos(4x), Pade', n_w=1.48 + 1.883j, profile='cos4x'),
    'thesis20a': _s('thesis Fig. 20a: tungsten (Ordal) 3.8313+2.9043i, sin(4x), Pade',
                    n_w=3.8313 + 2.9043j, profile='sin4x'),
    'thesis20b': _s('thesis Fig. 20b: iron (Ordal) 4.274+9.579i, sin(4x), Pade', n_w=4.274 + 9.579j, profile='sin4x'),
    'thesis21': _s('thesis Fig. 21: ZnGeP2 n=3.1874, cos(3x), alpha=0.01, Nx=Nz=64, N=M=13, q=1..3',
                   n_w=3.1874, profile='cos3x', alpha=0.01, Nx=64, Nz=64, M=13, q=(1, 2, 3), band_center='matlab'),
    'thesis22': _s('thesis Fig. 22: ZnO n=2.1054, sin(3x), alpha=0.01, Nx=Nz=64, N=M=13, q=1..3',
                   n_w=2.1054, profile='sin3x', alpha=0.01, Nx=64, Nz=64, M=13, q=(1, 2, 3), band_center='matlab'),
    'thesis23': _s('thesis Fig. 23: non-physical n_u=5, n_w=8.1, cos(x), eps_max=0.4, N=M=12, q=1',
                   n_u=5.0, n_w=8.1, eps_max=0.4, M=12, q=(1,), band_center='matlab'),
    'thesis24': _s('thesis Fig. 24: n_u=15, n_w=20i, alpha=0.1, cos(x), N=M=15, q=1, 1000x1000',
                   n_u=15.0, n_w=20j, alpha=0.1, q=(1,), n_eps=1000, n_delta=1000, band_center='matlab'),
    'thesis25': _s('thesis Fig. 25: n_u=10, n_w=40i, alpha=0.1, sin(x), N=M=20, q=1',
                   n_u=10.0, n_w=40j, alpha=0.1, profile='sinx', M=20, q=(1,), band_center='matlab'),
    'thesis26': _s('thesis Fig. 26: f_s1 = cos(4x)/4, eps_max=4, a=b=10, Nx=256, Nz=128, N=M=20, q=1',
                   n_w=1.1, profile='fs1', eps_max=4.0, a=10.0, b=10.0, Nx=256, Nz=128, M=20, q=(1,)),
    'thesis27': _s('thesis Fig. 27: f_s2 = exp(cos 3x)/3 - c0, eps_max=2, a=b=4, Nx=256, Nz=128, N=M=20, q=1',
                   n_w=1.1, profile='fs2', eps_max=2.0, a=4.0, b=4.0, Nx=256, Nz=128, M=20, q=(1,)),
    'thesis28a': _s('thesis Fig. 28a/b: rough f_r,P (P=120), eps_max=2, a=b=4, Nx=1024, Nz=128, N=M=20, q=1',
                    n_w=1.1, profile='rough120', eps_max=2.0, a=4.0, b=4.0, Nx=1024, Nz=128, M=20, q=(1,)),
    'thesis28c': _s('thesis Fig. 28c/d: Lipschitz f_L,P (P=120), as 28a',
                    n_w=1.1, profile='lipschitz120', eps_max=2.0, a=4.0, b=4.0, Nx=1024, Nz=128, M=20, q=(1,)),
    'thesis29': _s('thesis Fig. 29: TE, n_w=1.1, cos(x), Taylor, N=M=15', mode='TE', n_w=1.1, taylor=True),
    'thesis30': _s('thesis Fig. 30: TE, alpha=1e-4, 1000x1000, Taylor', mode='TE', n_w=1.1, taylor=True,
                   alpha=1e-4, n_eps=1000, n_delta=1000, band_center='matlab'),
    'thesis31a': _s('thesis Fig. 31a: TE, copper 0.94+1.337i, sin(5x), a=b=2/pi, Pade',
                    mode='TE', n_w=0.94 + 1.337j, profile='sin5x', a=2 / np.pi, b=2 / np.pi),
    'thesis31b': _s('thesis Fig. 31b: TE, cobalt 2.1396+3.9840i, sin(5x), a=b=2/pi, Pade',
                    mode='TE', n_w=2.1396 + 3.9840j, profile='sin5x', a=2 / np.pi, b=2 / np.pi),
    'thesis32': _s('thesis Fig. 32: TE, n_u=5, n_w=20i, alpha=0.001, cos(x), a=b=pi/2, q=1',
                   mode='TE', n_u=5.0, n_w=20j, alpha=0.001, a=np.pi / 2, b=np.pi / 2, q=(1,), band_center='matlab'),
    'thesis33a': _s('thesis Fig. 33a/b: TE, n_u=10, n_w=25i, alpha=0.1, sin(x), a=b=2/pi, q=1',
                    mode='TE', n_u=10.0, n_w=25j, alpha=0.1, profile='sinx', a=2 / np.pi, b=2 / np.pi,
                    q=(1,), band_center='matlab'),
    'thesis33c': _s('thesis Fig. 33c/d: as 33a with cos(x)', mode='TE', n_u=10.0, n_w=25j, alpha=0.1,
                    a=2 / np.pi, b=2 / np.pi, q=(1,), band_center='matlab'),
    # ---- new: dispersive materials from refractiveindex.info ------------------------------
    # Periods are chosen so that the first diffraction order of f = cos(x) excites surface
    # plasmons in the visible/NIR: grating-coupled SPP when n_spp * P / lambda = |alpha + m|,
    # n_spp = Re sqrt(eps_m eps_d / (eps_m + eps_d)).  For m = 1 this happens just beyond the
    # Rayleigh anomaly, i.e. for omega < 1: band q = 0 (omega in (0, 1), sub-wavelength grating),
    # which the MATLAB code never used.  max_delta = 0.05 re-evaluates the dispersive index
    # every 10 % in frequency.  Check: Ag, P = 0.5 um -> HOPS/AWE dip at 0.525 um vs. the
    # flat-surface SPP condition 0.524 um.
    'Ag_disp': _s('silver (J&C) with dispersion, period 0.5 um: grating-coupled surface plasmons at '
                  'lambda ~ 0.51 um (m=1), next to the Rayleigh anomaly at 0.5 um', n_w='Ag', period=0.5, q=(0, 1, 2), max_delta=0.05),
    'Al_uv': _s('aluminium (Rakic), period 0.25 um: UV plasmonics (SPP ~ 0.26 um)', n_w='Al', period=0.25, q=(0, 1, 2), max_delta=0.05),
    'Na': _s('sodium (Smith 1969): lowest-loss plasmonic metal, sharpest SPP dips; period 0.5 um',
             n_w='Na', period=0.5, q=(0, 1, 2), max_delta=0.05),
    'TiN': _s('titanium nitride (refractory plasmonics), period 0.6 um', n_w='TiN', period=0.6, q=(0, 1, 2), max_delta=0.05),
    'ITO_enz': _s('indium tin oxide, period 2 um: dielectric in VIS, epsilon-near-zero/metallic in NIR',
                  n_w='ITO', period=2.0, max_delta=0.05),
    'AZO_ir': _s('Al:ZnO, period 5 um: mid-IR plasmonic transparent conductor', n_w='AZO', period=5.0, max_delta=0.05),
    'VO2_cold': _s('VO2 at 25 C (insulating), period 5 um: compare with VO2_hot (switchable grating)',
                   n_w='VO2_cold', period=5.0, max_delta=0.05),
    'VO2_hot': _s('VO2 at 100 C (metallic), period 5 um', n_w='VO2_hot', period=5.0, max_delta=0.05),
    'Si': _s('crystalline silicon, period 1 um: high index, lossless in NIR, absorbing in VIS/UV',
             n_w='Si', period=1.0, windows='joint', max_delta=0.05),
    'TiO2': _s('rutile TiO2 (n ~ 2.5), period 1 um: dielectric metasurface material', n_w='TiO2', period=1.0,
               q=(1, 2), windows='joint', max_delta=0.05),
    'Ti_absorber': _s('titanium (n ~ k): lossy metal, broadband absorption', n_w='Ti', period=0.6, q=(0, 1, 2), max_delta=0.05),
    'SiC_reststrahlen': _s('crystalline SiC (Lorentz model), period 10.5 um, band q = 0: grating-coupled surface '
                           'phonon polaritons in the 10.3-12.5 um Reststrahlen band', n_w='SiC', period=10.5, q=(0,),
                           max_delta=0.02),
    'sapphire_reststrahlen': _s('sapphire (Querry), period 11 um, band q = 0: surface phonon polaritons in the '
                                '11-17 um Reststrahlen band', n_w='sapphire_IR', period=11.0, q=(0,), max_delta=0.02),
    'water_over_gold': _s('gold grating under water, period 0.6 um: grating-coupled SPR biosensor '
                          '(SPP dip ~ 0.83 um)', n_u='water', n_w='Au', period=0.6, q=(0, 1, 2), max_delta=0.05),
    'glass_over_silver': _s('silver grating under BK7 glass (n_u ~ 1.52), period 0.5 um', n_u='BK7', n_w='Ag',
                            period=0.5, q=(0, 1, 2), max_delta=0.05),
    'glass_over_air': _s('BK7 glass over air at 53 deg incidence (alpha = 6, q = 1): total internal reflection '
                         '(R_flat = 1), frustrated by diffraction into air', n_u='BK7', n_w='air', period=1.0,
                         alpha=6.0, q=(1,), max_delta=0.05),
    'high_contrast': _s('vacuum over silicon (lossless, n = 3.48 at 1.55 um), joint windows',
                        n_w='Si_IR', period=5.0, q=(1, 2), windows='joint', max_delta=0.05),
}

PROFILES = ('cosx', 'sinx', 'cos3x', 'sin3x', 'cos4x', 'sin4x', 'sin5x', 'fs1', 'fs2',
            'rough', 'lipschitz', 'rough120', 'lipschitz120', 'expr:<numpy expression in x>')


def profile_fn(name, xx):
    """Grating profiles f(x), f'(x) (thesis (6.4)-(6.40), (6.22)-(6.24))."""
    if name == 'cos4x':
        return np.cos(4 * xx), -4 * np.sin(4 * xx)
    if name == 'cosx':
        return np.cos(xx), -np.sin(xx)
    if name == 'sinx':
        return np.sin(xx), np.cos(xx)
    if name == 'sin4x':
        return np.sin(4 * xx), 4 * np.cos(4 * xx)
    if name == 'cos3x':
        return np.cos(3 * xx), -3 * np.sin(3 * xx)
    if name == 'sin3x':
        return np.sin(3 * xx), 3 * np.cos(3 * xx)
    if name == 'sin5x':
        return np.sin(5 * xx), 5 * np.cos(5 * xx)
    if name == 'fs1':
        return np.cos(4 * xx) / 4, -np.sin(4 * xx)
    if name == 'fs2':
        g = np.exp(np.cos(3 * xx)) / 3
        return g - g.mean(), -np.sin(3 * xx) * np.exp(np.cos(3 * xx))
    if name == 'lipschitz':
        return fourier_repr_lipschitz(40, xx)
    if name == 'rough':
        return fourier_repr_rough(40, xx)
    if name == 'lipschitz120':
        return fourier_repr_lipschitz(120, xx)
    if name == 'rough120':
        return fourier_repr_rough(120, xx)
    if name.startswith('expr:'):
        ns = {k: getattr(np, k) for k in ('sin', 'cos', 'exp', 'tanh', 'sinh', 'cosh', 'abs', 'pi', 'sqrt')}
        ns['x'] = xx
        f = np.asarray(eval(name[5:], {'__builtins__': {}}, ns), dtype=float) * np.ones_like(xx)
        f = f - f.mean()
        p = np.fft.fftfreq(xx.size, d=xx[1] - xx[0]) * 2 * np.pi
        return f, np.real(np.fft.ifft(1j * p * np.fft.fft(f)))
    raise ValueError(f'unknown profile {name}; choose one of {PROFILES}')


def _windows(qq, n_u_q, n_w_q, alpha, sigma, band_center, mode, min_width, verbose, max_delta=None):
    """Frequency windows (key, q, omega_bar, delta_max) on which separate HOPS/AWE expansions are made.

    'paper': one window per band q (paper eqs. (28)-(30), or MATLAB's centring).
    'joint': each band q is further split at the Rayleigh (Wood) frequencies of the LOWER layer,
             n_w omega = |alpha + p|, which lie inside it when n_w is real and larger than n_u.
             The delta-series of every quantity has a branch point there, so the AWE expansion
             of a single band cannot converge across them (energy defect ~1e-2...1e-4 for
             n_w = 1.5-3.5).  With 'joint' each sub-window is bounded by consecutive
             singularities of either layer and gets its own expansion (paper eq. (29) applied to
             both layers).  Costs one extra solve per lower-layer anomaly.
    max_delta: additionally split every window so that |delta| <= max_delta (finer piecewise-
             constant dispersion for named materials, and shorter AWE expansions).
    """
    out = []
    for q in qq:
        om_bar, dmax = _band_window(q, np.real(n_u_q[q]), alpha, sigma, band_center)
        if mode != 'joint' and not max_delta:
            out.append((f'q{q}', q, om_bar, dmax))
            continue
        lo, hi = om_bar * (1 - dmax / sigma), om_bar * (1 + dmax / sigma)     # the full band
        lo = max(lo, 0.02 * hi)                  # q = 0 starts at omega = 0 (static limit)
        n_w = n_w_q[q]
        cuts = []
        if mode == 'joint' and abs(np.imag(n_w)) < 1e-3 * abs(n_w) and np.real(n_w) > 0:
            cuts = [w for w in _wood_anomalies(np.real(n_w), alpha, np.array([lo, hi])) if lo < w < hi]
        edges = [lo] + cuts + [hi]
        segs = [(a_, b_) for a_, b_ in zip(edges[:-1], edges[1:])
                if (b_ - a_) >= min_width * (hi - lo)]          # drop slivers next to singularities
        j = 0
        for a_, b_ in segs:
            pieces = [a_, b_]
            if max_delta:
                md = max_delta / sigma          # geometric split: every piece has |delta| <= max_delta
                k = max(1, int(np.ceil(np.log(b_ / a_) / np.log((1 + md) / (1 - md)))))
                pieces = [a_] + list(a_ * (b_ / a_) ** (np.arange(1, k + 1) / k))
            for c_, e_ in zip(pieces[:-1], pieces[1:]):
                ob = 0.5 * (c_ + e_)
                out.append((f'q{q}w{j}', q, ob, sigma * (e_ - c_) / (2 * ob)))
                j += 1
        if verbose and cuts:
            print(f'band q={q}: split at lower-layer Wood anomalies omega = '
                  f'{", ".join(f"{w:.3f}" for w in cuts)} -> {j} windows')
    return out


def _band_window(q, n_u_re, alpha, sigma, band_center):
    """(omega_bar, delta_max) for band q."""
    if band_center == 'matlab':
        return (q + 0.5), sigma / (2 * q + 1)
    # paper eqs. (28)-(30) with d = 2 pi, c0 = 1
    return (alpha + q + 0.5) / n_u_re, sigma / (2 * alpha + 2 * q + 1)


def _wood_anomalies(n_re, alpha, omega):
    """Rayleigh/Wood frequencies n*omega = |alpha + p| inside [omega.min(), omega.max()]."""
    if n_re <= 0:
        return []
    lo, hi = omega.min() * n_re, omega.max() * n_re
    ps = np.arange(-int(hi) - 2, int(hi) + 3)
    return sorted({abs(alpha + p) / n_re for p in ps if lo < abs(alpha + p) < hi})


def _band(win, c):
    """Everything for one frequency window (module level so it can run in a worker process).
    win = (key, q, omega_bar, delta_max)."""
    key, q, omega_bar, dmax = win
    Nx, Nz, N, M, d, c_0, alpha_bar = c['Nx'], c['Nz'], c['N'], c['M'], c['d'], c['c_0'], c['alpha_bar']
    f, f_x, Eps, N_Eps, N_delta, sigma = (c['f'], c['f_x'], c['Eps'], c['N_Eps'], c['N_delta'], c['sigma'])
    n_u, n_w = c['n_u_w'][key], c['n_w_w'][key]
    identy = np.eye(Nz + 1)
    delta = np.array([0.0]) if N_delta == 1 else np.linspace(-dmax, dmax, N_delta)
    omega = (1 + delta) * omega_bar
    lam = 2 * np.pi * c_0 / omega
    k_u_bar = n_u * omega_bar / c_0
    gamma_u_bar = csqrt(k_u_bar ** 2 - alpha_bar ** 2)
    xx, pp, alpha_bar_p, gamma_u_bar_p, _, _ = setup_2d(Nx, d, alpha_bar, gamma_u_bar)
    k_w_bar = n_w * omega_bar / c_0
    gamma_w_bar = csqrt(k_w_bar ** 2 - alpha_bar ** 2)
    xx, pp, alpha_bar_p, gamma_w_bar_p, _, _ = setup_2d(Nx, d, alpha_bar, gamma_w_bar)
    zeta_n_m, psi_n_m = setup_zeta_psi_n_m(xx, pp, alpha_bar, gamma_u_bar, f, f_x, Nx, N, M)
    tau2 = 1.0 if c['Mode'] == 1 else (n_u / n_w) ** 2

    t0 = time.time()
    solve = SOLVERS[c['solver']]
    U_n_m, W_n_m, ubar_n_m, wbar_n_m = solve(
        tau2, zeta_n_m, psi_n_m, gamma_u_bar_p, gamma_w_bar_p, N, Nx, f, f_x, pp, alpha_bar,
        gamma_u_bar, gamma_w_bar, c['Dz'], c['a'], c['b'], Nz, M, identy, alpha_bar_p)
    t1 = time.time()
    ub = np.transpose(ubar_n_m, (2, 1, 0))    # permute(.,[3 2 1]) -> (N+1, M+1, Nx)
    wb = np.transpose(wbar_n_m, (2, 1, 0))
    ee_flat, ru_flat, rl_flat = energy_defect(tau2, ub, wb, d, alpha_bar, gamma_u_bar, gamma_w_bar,
                                              Eps, delta, Nx, 0, 0, N_Eps, N_delta, 1)
    st = 1 if c['Taylor'] else 2
    ee, ru, rl = energy_defect(tau2, ub, wb, d, alpha_bar, gamma_u_bar, gamma_w_bar,
                               Eps, delta, Nx, N, M, N_Eps, N_delta, st,
                               taylor_full_order=c.get('taylor_full_order', False))
    wood_w = _wood_anomalies(np.real(n_w), alpha_bar, omega) if abs(np.imag(n_w)) < 1e-3 * abs(n_w) else []
    if c['verbose']:
        extra = (f'; lower-layer Wood anomalies inside the band at omega = '
                 f'{", ".join(f"{w:.3f}" for w in wood_w[:6])}' if wood_w else '')
        print(f'{key}: omega in [{omega.min():.3f}, {omega.max():.3f}], n_u = {complex(n_u):.4g}, n_w = {complex(n_w):.4g}; two_layer_solve ({c["solver"]}) '
              f'{t1 - t0:.1f} s, energy_defect {time.time() - t1:.1f} s{extra}', flush=True)
    RR = ru / ru_flat if c['relative'] else ru
    extra_fields = {}
    if c.get('keep_fields'):
        # volume fields u_{n,m}(x,z'), w_{n,m}(x,z') driven by the interface data U_{n,m}, W_{n,m}
        # (for field pictures/movies, see refl_movie.py)
        from hops import field_tfe_helmholtz_m_and_n, field_tfe_helmholtz_m_and_n_lf
        extra_fields = dict(
            u_n_m=field_tfe_helmholtz_m_and_n(U_n_m, f, pp, gamma_u_bar_p, alpha_bar, gamma_u_bar, c['Dz'],
                                              c['a'], Nx, Nz, N, M, identy, alpha_bar_p),
            w_n_m=field_tfe_helmholtz_m_and_n_lf(W_n_m, f, pp, gamma_w_bar_p, alpha_bar, gamma_w_bar, c['Dz'],
                                                 c['b'], Nx, Nz, N, M, identy, alpha_bar_p),
            gamma_u_bar=gamma_u_bar, gamma_w_bar=gamma_w_bar, tau2=tau2, dmax=dmax)
    return dict(**extra_fields, q=q, key=key, omega_bar=omega_bar, delta=delta, omega=omega, lam=lam, Eps=Eps, ee=ee, ru=ru, rl=rl,
                ru_flat=ru_flat, RR=RR, U_n_m=U_n_m, W_n_m=W_n_m, n_u=n_u, n_w=n_w,
                ubar_n_m=ubar_n_m, wbar_n_m=wbar_n_m, wood_w=np.array(wood_w))


def _resolve_indices(spec, qq, alpha, sigma, band_center, period, index_mode, lambda_ref, db, n_u_spec=None):
    """n for every band: literal -> constant; material -> evaluated at band-centre wavelengths."""
    lit = mat.parse_index(spec)
    if lit is not None:
        return {q: lit for q in qq}, None
    if period is None and lambda_ref is None:
        info = mat.MATERIALS.get(spec, {})
        period = info.get('period') or 1.0
        print(f"note: '{spec}' is a named material; using its suggested grating period {period} um "
              f"(set --period to change)")
    out, lams = {}, {}
    for q in qq:
        if index_mode == 'fixed':
            L = lambda_ref if lambda_ref is not None else period / (qq[0] + 0.5)
        else:
            nu_re = 1.0 if n_u_spec is None else np.real(n_u_spec[q])
            omega_bar, _ = _band_window(q, nu_re, alpha, sigma, band_center)
            L = period / omega_bar
        out[q] = mat.refractive_index(spec, L, db=db)
        lams[q] = L
    return out, lams


def run(scenario=SCENARIO, qq=None, N_Eps=None, N_delta=None, RunNumber=RunNumber, Mode=None,
        Taylor=None, verbose=True, solver=None, workers=None, **over):
    """Compute the reflectivity map for a scenario (+ keyword overrides of any scenario key:
    n_u, n_w, profile, eps_max, alpha, M, N, Nx, Nz, a, b, mode, taylor, period, band_center,
    index_mode, lambda_ref, rii_db, relative, sigma, windows)."""
    solver = solver or SOLVER
    workers = WORKERS if workers is None else workers
    sc = dict(SCENARIOS[scenario])
    sc.update({k: v for k, v in over.items() if v is not None})
    if Taylor is not None:
        sc['taylor'] = Taylor
    if Mode is not None:
        sc['mode'] = 'TE' if Mode == 1 else 'TM'
    qq = tuple(qq) if qq is not None else tuple(sc['q'])
    N_Eps = N_Eps or sc.get('n_eps', 100)
    N_delta = N_delta or sc.get('n_delta', 100)
    sigma = sc.get('sigma', 0.99)
    M, Nx, Eps_Max = sc['M'], sc['Nx'], sc['eps_max']
    if RunNumber == 1:
        M, Nx, Eps_Max, sigma = 4, 16, 1e-2, 1e-2
    elif RunNumber == 2:
        M, Nx, Eps_Max, sigma = 6, 24, 0.1, 0.1
    elif RunNumber == 3:
        M, Nx, Eps_Max, sigma = 8, 32, 0.1, 0.5
    N = sc.get('N') or M
    Nz = sc['Nz']
    alpha_bar = sc['alpha']
    d, c_0 = 2 * np.pi, 1.0
    band_center = sc.get('band_center', 'auto')
    band_center = 'rayleigh' if band_center == 'auto' else band_center
    index_mode = sc.get('index_mode', 'band')
    period, lambda_ref, db = sc.get('period'), sc.get('lambda_ref'), sc.get('rii_db')

    # superstrate first (its real part sets the band centres), then the substrate
    n_u_q, lam_u = _resolve_indices(sc['n_u'], qq, alpha_bar, sigma, band_center, period, index_mode,
                                    lambda_ref, db)
    if lam_u is not None and period is None:
        period = mat.MATERIALS.get(sc['n_u'], {}).get('period', 1.0)
    n_w_q, lam_w = _resolve_indices(sc['n_w'], qq, alpha_bar, sigma, band_center, period, index_mode,
                                    lambda_ref, db, n_u_spec=n_u_q)
    if lam_w is not None and period is None:
        period = mat.MATERIALS.get(sc['n_w'], {}).get('period', 1.0)
    for q in qq:
        if abs(np.imag(n_u_q[q])) > 1e-3:
            warnings.warn(f'band {q}: absorbing superstrate n_u = {n_u_q[q]:.4g}; the incident wave decays')

    Dz, _ = cheb(Nz)
    Eps = np.linspace(0, Eps_Max, N_Eps)
    xx = (d / Nx) * np.arange(Nx)
    f, f_x = profile_fn(sc['profile'], xx)
    common = dict(N_delta=N_delta, sigma=sigma, c_0=c_0, d=d, n_u_q=n_u_q, n_w_q=n_w_q, alpha_bar=alpha_bar,
                  Nx=Nx, Nz=Nz, N=N, M=M, f=f, f_x=f_x, Mode=1 if sc['mode'] == 'TE' else 2, Dz=Dz,
                  a=sc['a'], b=sc['b'], Eps=Eps, N_Eps=N_Eps, Taylor=sc['taylor'], solver=solver,
                  verbose=verbose, band_center=band_center, relative=sc.get('relative', PlotRelative),
                  keep_fields=sc.get('keep_fields', False),
                  taylor_full_order=sc.get('taylor_full_order', False))
    wins = _windows(qq, n_u_q, n_w_q, alpha_bar, sigma, band_center, sc.get('windows', 'paper'),
                    sc.get('min_window', 0.02), verbose, sc.get('max_delta'))
    if sc.get('omega_range'):             # only the windows that overlap a frequency range
        o_lo, o_hi = sc['omega_range']
        wins = [w for w in wins if w[2] * (1 + w[3]) >= o_lo and w[2] * (1 - w[3]) <= o_hi]
    # indices per window (named materials: at the window-centre wavelength)
    def per_window(spec, nq):
        if mat.parse_index(spec) is not None or index_mode == 'fixed':
            return {w[0]: nq[w[1]] for w in wins}
        return {w[0]: mat.refractive_index(spec, period / w[2], db=db) for w in wins}
    common['n_u_w'] = per_window(sc['n_u'], n_u_q)
    common['n_w_w'] = per_window(sc['n_w'], n_w_q)
    workers = min(workers or 1, len(wins))
    if workers > 1:
        # the frequency windows are independent -> solve them in parallel processes
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(_band, wins, [common] * len(wins)))
    else:
        results = [_band(w, common) for w in wins]
    lossless = all(abs(np.imag(results[i]['n_w'])) < 1e-6 and abs(np.imag(results[i]['n_u'])) < 1e-6
                   for i in range(len(results)))
    info = dict(scenario=scenario, N=N, M=M, Nx=Nx, Nz=Nz, n_w=n_w_q[qq[0]], n_u=n_u_q[qq[0]],
                n_w_spec=sc['n_w'], n_u_spec=sc['n_u'], Taylor=sc['taylor'], mode=sc['mode'],
                period=period, profile=sc['profile'], eps_max=Eps_Max, alpha=alpha_bar,
                lossless=lossless, relative=common['relative'], desc=sc.get('desc', ''),
                a=sc['a'], b=sc['b'], f=f, f_x=f_x)
    return results, info


def plot(results, info, outdir=OUTDIR, tag=None):
    """MATLAB-like rendering: each band q is its own contourf call with automatic
    levels, all sharing one colour scale ('colormap hot'), as in refl_map.m.
    For R, isolated Pade spikes above 1 (a handful of points out of 60,000 next to
    the plasmon resonance) are drawn in the top colour, so the colour scale ends at 1."""
    from hops.plotting import matlab_contourf, shared_colorbar
    from matplotlib.colors import Normalize
    os.makedirs(outdir, exist_ok=True)
    tag = tag or info['scenario']
    period = info.get('period')
    if PlotLambda:
        if period:
            xfun = lambda r: r['lam'] * period / (2 * np.pi)
            xlabel = r'$\lambda$ ($\mu$m)' + f'   (period {period:g} $\\mu$m)'
        else:
            xfun = lambda r: r['lam']
            xlabel = r'$\lambda$'
    else:
        xfun = lambda r: r['omega']
        xlabel = r'$\omega$'
    relative = info.get('relative', PlotRelative)
    sub = ''
    if 'n_w_spec' in info:
        def _fmt(v):
            vv = mat.parse_index(v)
            return f'{vv:.4g}' if vv is not None else str(v)
        sub = (f"  {info.get('mode', 'TM')}, $n^u$={_fmt(info['n_u_spec'])}, $n^w$={_fmt(info['n_w_spec'])}, "
               f"{'Taylor' if info.get('Taylor') else 'Pade'}")
    figs = []
    for num, (title, getZ) in enumerate([('$D$', lambda r: safe_log10(r['ee'])),
                                         ('$R/R_{flat}$' if relative else '$R$', lambda r: np.real(r['RR']))],
                                        start=1):
        Zs = [getZ(r) for r in results]
        allz = np.concatenate([z[np.isfinite(z)].ravel() for z in Zs])
        vmin, vmax = allz.min(), allz.max()
        if num == 2 and np.mean(allz > 1.0 + 1e-9) < 1e-3 and vmax > 1.0:
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
        ax.set_title(title + sub + note, fontsize=11)
        fig.tight_layout()
        fname = os.path.join(outdir, f'refl_map_{tag}_{"D" if num == 1 else "R"}.png')
        fig.savefig(fname, dpi=130)
        figs.append(fname)
    return figs


def load_results(scenario, outdir=OUTDIR):
    """Re-load a saved run (figures/refl_map_<scenario>.npz) for re-plotting."""
    d = np.load(os.path.join(outdir, f'refl_map_{scenario}.npz'))
    keys = sorted({k.split('_')[0] for k in d.files if k.startswith('q')})
    res = [dict(key=kk, **{k: d[f'{kk}_{k}'] for k in ('lam', 'omega', 'Eps', 'ee', 'ru', 'RR')}) for kk in keys]
    info = dict(scenario=scenario)
    if 'period' in d.files and np.isfinite(d['period']):
        info['period'] = float(d['period'])
    return res, info


def print_materials():
    print(f"{'category':22s} {'key':13s} {'data range (um)':19s} {'P(um)':>5s}   n+ik at "
          f"0.4 / 0.6 / 1.0 / 1.55 / 10 um")
    for cat, key, rng, per, vals, note in mat.list_materials():
        print(f'{cat:22s} {key:13s} {rng:19s} {per:5g}   ' + ' '.join(f'{v:>14s}' for v in vals))
        print(f'{"":22s} {"":13s} {note}')
    print("\nAny literal index also works (e.g. --nw 0.05+2.275i, --nw 20i), and any page of the full "
          "refractiveindex.info database via --rii-db <path>/database/data --nw main/Au/nk/Olmon-sc.yml")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scenario', default=SCENARIO, choices=list(SCENARIOS))
    ap.add_argument('--list-scenarios', action='store_true')
    ap.add_argument('--list-materials', action='store_true')
    g = ap.add_argument_group('physics (override the scenario)')
    g.add_argument('--nu', help='upper index n^u: number (1.5, 0.05+2.275i, 20i) or material key/path')
    g.add_argument('--nw', help='lower index n^w: number or material key/path')
    g.add_argument('--mode', choices=['TE', 'TM'])
    g.add_argument('--alpha', type=float, help='Bloch/incidence parameter alpha')
    g.add_argument('--period', type=float, help='grating period in micrometres (maps bands to wavelengths)')
    g.add_argument('--index-mode', choices=['band', 'fixed'], help='dispersion handling for named materials')
    g.add_argument('--lambda-ref', type=float, help='wavelength (um) for --index-mode fixed')
    g.add_argument('--rii-db', help='path to refractiveindex.info database/data (for arbitrary pages)')
    g = ap.add_argument_group('geometry / numerics')
    g.add_argument('--profile', help=f'grating profile: {", ".join(PROFILES)}')
    g.add_argument('--eps-max', type=float)
    g.add_argument('--a', type=float, help='upper artificial boundary z = a')
    g.add_argument('--b', type=float, help='lower artificial boundary z = -b')
    g.add_argument('--M', type=int, help='frequency order M (and N = M unless --N)')
    g.add_argument('--N', type=int, help='height order N')
    g.add_argument('--Nx', type=int)
    g.add_argument('--Nz', type=int)
    g.add_argument('--summation', choices=['taylor', 'pade'])
    g.add_argument('--band-center', choices=['auto', 'matlab', 'rayleigh'])
    g.add_argument('--taylor-full-order', action='store_true',
                   help='Taylor: sum all min(N,M)+1 polar orders (MATLAB taylorsum_2_coeff sums only half of them)')
    g.add_argument('--max-delta', type=float,
                   help='split windows so that |delta| <= this (finer dispersion sampling), e.g. 0.05')
    g.add_argument('--windows', choices=['paper', 'joint'],
                   help="'joint': also split each band at the lower layer's Wood anomalies (lossless high-index n_w)")
    g.add_argument('--sigma', type=float, help='filling fraction of the delta window (0.99)')
    g.add_argument('--q', type=int, nargs='*', help='frequency bands (default from the scenario)')
    g.add_argument('--neps', type=int)
    g.add_argument('--ndelta', type=int)
    g.add_argument('--absolute', action='store_true', help='plot R instead of R/R_flat')
    g = ap.add_argument_group('run')
    g.add_argument('--solver', default=SOLVER, choices=list(SOLVERS))
    g.add_argument('--workers', type=int, default=WORKERS)
    g.add_argument('--tag', help='name used for the output files (default: scenario)')
    g.add_argument('--no-show', action='store_true')
    args = ap.parse_args(argv)

    if args.list_scenarios:
        w = max(map(len, SCENARIOS))
        for k, v in SCENARIOS.items():
            print(f'{k:{w}s}  {v["desc"]}')
        return
    if args.list_materials:
        print_materials()
        return
    if args.no_show:
        matplotlib.use('Agg')
    over = dict(n_u=args.nu, n_w=args.nw, mode=args.mode, alpha=args.alpha, period=args.period,
                index_mode=args.index_mode, lambda_ref=args.lambda_ref, rii_db=args.rii_db,
                profile=args.profile, eps_max=args.eps_max, a=args.a, b=args.b, M=args.M, N=args.N,
                Nx=args.Nx, Nz=args.Nz, band_center=args.band_center, sigma=args.sigma, windows=args.windows, max_delta=args.max_delta,
                taylor=None if args.summation is None else args.summation == 'taylor',
                relative=False if args.absolute else None,
                taylor_full_order=True if args.taylor_full_order else None)
    t0 = time.time()
    res, info = run(args.scenario, args.q, args.neps, args.ndelta, solver=args.solver,
                    workers=args.workers, **over)
    custom = any(v is not None for k, v in over.items() if k != 'relative')
    clean = lambda v: os.path.splitext(os.path.basename(str(v)))[0].replace('+', 'p').replace('.', 'd')
    tag = args.tag or (args.scenario if not custom else
                       f"{clean(info['n_u_spec'])}_over_{clean(info['n_w_spec'])}_{info['mode']}_{info['profile']}")
    for fn in plot(res, info, tag=tag):
        print('saved', fn)
    np.savez_compressed(os.path.join(OUTDIR, f'refl_map_{tag}.npz'),
                        period=np.nan if info['period'] is None else info['period'],
                        **{f'{r["key"]}_{k}': r[k] for r in res for k in ('lam', 'omega', 'Eps', 'ee', 'ru', 'RR')})
    print(f'total {time.time() - t0:.1f} s')
    if SHOW and not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
