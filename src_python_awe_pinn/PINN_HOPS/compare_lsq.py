"""compare_lsq.py -- least-squares interface PINN (random-feature / ELM PINN) vs HOPS/AWE.

Parts (each writes into results/lsq/):
  points   the five single-point cases of compare_point.py: R, T, D, interface field, PINN loss, time
  conv     error in R and D vs number of network features (dielectric, gold), for tanh / sin / gauss and
           the I-PINN choice (tanh above, sin below); and the generic 'rfm' network for contrast
  coords   features in physical (x, z) vs in the flattened TFE coordinate of each layer (gold, cos 4x)
  map      reflectivity map and energy defect on the refl_map band q = 1 (one LSQ solve per grid point)
           vs refl_map.py (HOPS/AWE)

    python compare_lsq.py                 # all parts
    python compare_lsq.py --parts points conv
"""
import argparse
import json
import os
import time
from dataclasses import replace

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from pinn_hops import Grating2D
from pinn_hops.lsq_pinn import LSQPinn
from pinn_hops.hops_reference import hops_point, hops_refl_map
from compare_point import CASES

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'lsq')

# network size per profile: Fourier modes K and z-features per mode (separable basis)
# network size per profile: Fourier modes K, z-features per mode nz (separable basis), collocation Mx x Ms;
# coords='tfe': features in the flattened coordinate of each layer (see lsq_pinn.MappedFeatures and 'coords')
CFG = {'cosx': dict(K=10, nz=12, Mx=64, Ms=32), 'cos4x': dict(K=32, nz=12, Mx=128, Ms=24, coords='tfe')}
ACT = {'cosx': 'sin', 'cos4x': 'sin'}           # fastest-converging activation (see the conv part)


def make(G, K, nz, Mx, Ms, activation='tanh', basis='separable', rz=3.0, **kw):
    n = nz * (2 * K + 1) if basis == 'separable' else nz
    return LSQPinn(G, K=K, n_feat=n, Mx=Mx, Ms=Ms, activation=activation, basis=basis, rz=rz, **kw)


def run_points(cases=None):
    rows = []
    for name in cases or CASES:
        G = Grating2D(**CASES[name]['G'])
        ref = hops_point(G, N=16, Nx=128, Nz=48, summation='pade')
        t0 = time.time()
        hops_point(G, N=16, Nx=128, Nz=48, summation='pade', fields=False)
        t_h = time.time() - t0
        t0 = time.time()
        P = make(G, activation=ACT[G.profile], **CFG[G.profile]).solve()
        R, T, D = P.energy()
        t_l = time.time() - t0
        U = P.evaluate('u', ref['x'], G.g(ref['x']))
        W = P.evaluate('w', ref['x'], G.g(ref['x']))
        rows.append(dict(case=name, R_hops=float(ref['R']), R_lsq=float(R), R_relerr=float(abs(R - ref['R']) / ref['R']),
                         D_hops=float(ref['D']), D_lsq=float(D), D_absdiff=float(abs(D - ref['D'])),
                         U_relerr=float(np.abs(U - ref['U']).max() / np.abs(ref['U']).max()),
                         W_relerr=float(np.abs(W - ref['W']).max() / np.abs(ref['W']).max()),
                         loss=P.loss, unknowns=P.shape[1], rows=P.shape[0], rank=int(P.rank),
                         t_lsq_s=t_l, t_hops_s=t_h))
        print(rows[-1], flush=True)
    os.makedirs(OUT, exist_ok=True)
    json.dump(rows, open(os.path.join(OUT, 'points.json'), 'w'), indent=1)
    return rows


def run_conv():
    """convergence in the number of features, for several activations"""
    out = {}
    for name, nzs in (('dielectric', (2, 4, 8, 12, 16, 24, 32)), ('gold', (2, 4, 8, 12, 16))):
        G = Grating2D(**CASES[name]['G'])
        ref = hops_point(G, N=16, Nx=128, Nz=48, summation='pade', fields=False)
        cfg = dict(CFG[G.profile], coords='physical')
        cfg['K'] = 40 if G.profile == 'cos4x' else cfg['K']
        cfg['Mx'] = 160 if G.profile == 'cos4x' else cfg['Mx']
        for act in ('tanh', 'sin', 'gauss', ('tanh', 'sin')):
            key = f"{name}|{act if isinstance(act, str) else 'ipinn_' + '_'.join(act)}"
            out[key] = []
            for nz in nzs:
                cfg['nz'] = nz
                P = make(G, activation=act, **cfg).solve()
                R, T, D = P.energy()
                out[key].append(dict(nz=nz, unknowns=P.shape[1], R_relerr=abs(R - ref['R']) / ref['R'],
                                     D_absdiff=abs(D - ref['D']), loss=P.loss, t=P.solve_time))
                print(key, out[key][-1], flush=True)
        # generic (non-separable) random-feature network, for contrast
        key = f'{name}|rfm_tanh'
        out[key] = []
        for n in (100, 200, 400, 800):
            P = LSQPinn(G, K=cfg['K'], n_feat=n, basis='rfm', Mx=cfg['Mx'], Ms=32, rz=3.0).solve()
            R, T, D = P.energy()
            out[key].append(dict(nz=n, unknowns=P.shape[1], R_relerr=abs(R - ref['R']) / ref['R'],
                                 D_absdiff=abs(D - ref['D']), loss=P.loss, t=P.solve_time))
            print(key, out[key][-1], flush=True)
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(os.path.join(OUT, 'conv.json'), 'w'), indent=1, default=float)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, name in zip(axs, ('dielectric', 'gold')):
        for key, v in out.items():
            if key.startswith(name):
                ax.loglog([r['unknowns'] for r in v], [max(r['R_relerr'], 1e-16) for r in v], 'o-', label=key.split('|')[1])
        pinn = _gradient_pinn_err(name)
        if pinn:
            ax.axhline(pinn, color='k', ls='--', label='gradient-trained PINN (Adam+L-BFGS)')
        ax.set_xlabel('unknowns (complex output weights, both subdomains)')
        ax.set_ylabel('relative error in R vs HOPS/AWE')
        ax.set_title(name)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    fig.suptitle('Least-squares interface PINN: convergence in the network size')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'conv.png'), dpi=110)
    plt.close(fig)
    return out


def _gradient_pinn_err(name):
    fn = os.path.join(os.path.dirname(OUT), 'point', 'summary.json')
    if os.path.exists(fn):
        for r in json.load(open(fn)):
            if r['case'] == name:
                return r['R_relerr']
    return None


def run_map(scenario='dielectric', n_grid=21, q=1):
    sc = {'dielectric': dict(n_w=1.1, profile='cosx'), 'gold': dict(n_w=1.48 + 1.883j, profile='cos4x'),
          'silver': dict(n_w=0.05 + 2.275j, profile='cos4x')}[scenario]
    t0 = time.time()
    res, info = hops_refl_map(scenario, q=(q,), N_Eps=n_grid, N_delta=n_grid)
    t_h = time.time() - t0
    r = res[0]
    Eps, omega = r['Eps'], r['omega']
    Rh, Dh, Rflat = np.real(r['ru']), np.real(r['ee']), np.real(r['ru_flat'])
    cfg = dict(CFG[sc['profile']])
    Rl, Dl = np.zeros_like(Rh), np.zeros_like(Rh)
    t0 = time.time()
    for i, e in enumerate(Eps):
        for j, o in enumerate(omega):
            G = Grating2D(eps=float(e), omega=float(o), **sc)
            P = make(G, activation=ACT[sc['profile']], **cfg).solve()
            Rl[i, j], _, Dl[i, j] = P.energy()
    t_l = time.time() - t0
    summ = dict(scenario=scenario, grid=f'{n_grid}x{n_grid}', t_hops_s=t_h, t_lsq_s=t_l,
                R_absdiff_max=float(np.abs(Rl - Rh).max()), R_absdiff_median=float(np.median(np.abs(Rl - Rh))),
                R_reldiff_median=float(np.median(np.abs(Rl - Rh) / np.abs(Rh))),
                D_hops_median_log10=float(np.log10(np.median(np.abs(Dh)))),
                D_lsq_median_log10=float(np.log10(np.median(np.abs(Dl)))),
                D_lsq_max_log10=float(np.log10(np.abs(Dl).max())), lossless=bool(info['lossless']))
    print(summ, flush=True)
    os.makedirs(OUT, exist_ok=True)
    tag = f'map_{scenario}_q{q}'
    json.dump(summ, open(os.path.join(OUT, tag + '.json'), 'w'), indent=1)
    np.savez_compressed(os.path.join(OUT, tag + '.npz'), Eps=Eps, omega=omega, R_hops=Rh, R_lsq=Rl, D_hops=Dh,
                        D_lsq=Dl, R_flat=Rflat)
    lam = 2 * np.pi / omega
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.5))
    RR = [Rh / Rflat, Rl / Rflat]
    lo, hi = min(z.min() for z in RR), max(z.max() for z in RR)
    for ax, Z, t in ((axs[0, 0], RR[0], 'HOPS/AWE (refl_map.py)'), (axs[0, 1], RR[1], 'least-squares PINN')):
        m = ax.contourf(lam, Eps, Z, levels=np.linspace(lo, hi, 15), cmap='hot')
        ax.contour(lam, Eps, Z, levels=np.linspace(lo, hi, 15)[1:-1], colors='k', linewidths=0.4)
        fig.colorbar(m, ax=ax)
        ax.set_title(r'$R/R_{flat}$: ' + t)
    m = axs[0, 2].pcolormesh(lam, Eps, np.log10(np.abs(Rl - Rh) + 1e-17), cmap='viridis', shading='auto')
    fig.colorbar(m, ax=axs[0, 2])
    axs[0, 2].set_title(r'$\log_{10}|R_{LSQ} - R_{HOPS}|$')
    Ds = [np.log10(np.abs(Dh) + 1e-17), np.log10(np.abs(Dl) + 1e-17)]
    lo, hi = min(z.min() for z in Ds), max(z.max() for z in Ds)
    lab = 'log10|D|' if info['lossless'] else 'log10 A (= D)'
    for ax, Z, t in ((axs[1, 0], Ds[0], 'HOPS/AWE'), (axs[1, 1], Ds[1], 'least-squares PINN')):
        m = ax.contourf(lam, Eps, Z, levels=np.linspace(lo, hi, 15), cmap='hot')
        fig.colorbar(m, ax=ax)
        ax.set_title(f'{lab}: {t}')
    axs[1, 2].axis('off')
    axs[1, 2].text(0, 1, '\n'.join(f'{k}: {v:.4g}' if isinstance(v, float) else f'{k}: {v}' for k, v in summ.items()),
                   va='top', family='monospace', fontsize=9)
    for ax in axs.ravel()[:5]:
        ax.set_xlabel(r'$\lambda$')
        ax.set_ylabel(r'$\varepsilon$')
    fig.suptitle(f'{scenario}, band q={q}: HOPS/AWE vs least-squares interface PINN (same equations (6))')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, tag + '.png'), dpi=100)
    plt.close(fig)
    return summ


def run_coords():
    """physical vs flattened (TFE) coordinates in the features, gold f = cos 4x, eps = 0.1 and 0.2"""
    out = []
    for eps in (0.1, 0.2):
        G = Grating2D(eps=eps, n_w=1.48 + 1.883j, profile='cos4x')
        ref = hops_point(G, N=24, Nx=128, Nz=48, summation='pade', fields=False)
        for coords in ('physical', 'tfe'):
            for K in (16, 24, 32, 40):
                P = make(G, K=K, nz=12, Mx=max(128, 4 * K), Ms=24, activation='sin', coords=coords).solve()
                R, T, D = P.energy()
                out.append(dict(eps=eps, coords=coords, K=K, unknowns=P.shape[1],
                                R_relerr=float(abs(R - ref['R']) / ref['R']), loss=P.loss, t=P.solve_time))
                print(out[-1], flush=True)
    os.makedirs(OUT, exist_ok=True)
    json.dump(out, open(os.path.join(OUT, 'coords.json'), 'w'), indent=1)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for eps, ls in ((0.1, '-'), (0.2, '--')):
        for coords, c in (('physical', 'C0'), ('tfe', 'C3')):
            v = [r for r in out if r['eps'] == eps and r['coords'] == coords]
            ax.semilogy([r['K'] for r in v], [r['R_relerr'] for r in v], 'o' + ls, color=c,
                        label=f'{coords} (x, z), eps = {eps}' if coords == 'physical' else f'TFE (x, zeta), eps = {eps}')
    ax.set_xlabel('Fourier modes K in the network')
    ax.set_ylabel('relative error in R vs HOPS')
    ax.set_title('gold, f = cos 4x: network coordinates')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'coords.png'), dpi=110)
    plt.close(fig)
    return out


def check_map_edges(scenario='dielectric', q=1):
    """Where LSQ-PINN and the HOPS/AWE map disagree: re-solve those (eps, omega) points with HOPS at
    delta = 0 (no frequency expansion, N = 24).  Shows which of the two maps carries the error."""
    tag = f'map_{scenario}_q{q}'
    d = np.load(os.path.join(OUT, tag + '.npz'))
    sc = {'dielectric': dict(n_w=1.1, profile='cosx'), 'gold': dict(n_w=1.48 + 1.883j, profile='cos4x')}[scenario]
    Nx = 32 if sc['profile'] == 'cosx' else 128
    rows = []
    n = len(d['Eps']) - 1
    m = len(d['omega']) - 1
    for i, j in ((n, 0), (n, m), (n // 2, 0), (n // 2, m), (n, m // 2)):
        e, o = float(d['Eps'][i]), float(d['omega'][j])
        r = hops_point(Grating2D(eps=e, omega=o, **sc), N=24, Nx=Nx, Nz=32, summation='pade', fields=False)
        rows.append(dict(eps=e, omega=o, R_hops_point=float(r['R']), R_lsq=float(d['R_lsq'][i, j]),
                         R_awe_map=float(d['R_hops'][i, j]), diff_point_lsq=float(abs(r['R'] - d['R_lsq'][i, j])),
                         diff_point_awe=float(abs(r['R'] - d['R_hops'][i, j]))))
        print(rows[-1], flush=True)
    json.dump(rows, open(os.path.join(OUT, tag + '_edges.json'), 'w'), indent=1)
    return rows


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--parts', nargs='*', default=['points', 'conv', 'coords', 'map'])
    ap.add_argument('--map-scenarios', nargs='*', default=['dielectric', 'gold'])
    ap.add_argument('--grid', type=int, default=21)
    a = ap.parse_args()
    matplotlib.use('Agg')
    if 'points' in a.parts:
        run_points()
    if 'conv' in a.parts:
        run_conv()
    if 'coords' in a.parts:
        run_coords()
    if 'map' in a.parts:
        for s in a.map_scenarios:
            run_map(s, a.grid if s == 'dielectric' else min(a.grid, 11))
            check_map_edges(s)
