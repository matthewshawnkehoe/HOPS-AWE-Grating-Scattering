"""inverse_demo.py -- an inverse problem: recover the grating amplitude eps and the substrate index n_w
from measured diffraction efficiencies, with HOPS and with the least-squares interface PINN as the
forward model.

Data: the reflected efficiencies e_p = (gamma^u_p / gamma^u) |xi_p|^2 of the propagating orders
(p = -1, 0, 1) at 5 frequencies across the refl_map band q = 1 (15 numbers), computed by HOPS/AWE with a
fine discretisation ("measurements"), optionally with relative Gaussian noise.
Unknowns: eps (true 0.13) and n_w (true 1.1).
Solver: scipy.optimize.least_squares (trust-region reflective, finite-difference Jacobian); the forward
model is either HOPS (hops_point, N = 16, Pade) or the LSQ-PINN (compare_lsq.CFG).

    python inverse_demo.py            # noise-free and 1 % noise, both forward models
Writes results/inverse/inverse.json and inverse.png.
"""
import json
import os
import time

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

from pinn_hops import Grating2D
from pinn_hops.hops_reference import hops_point
from compare_lsq import make, CFG

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'inverse')
OMEGAS = 1.5 * (1 + np.linspace(-0.3, 0.3, 5))
TRUE = dict(eps=0.13, n_w=1.1)
START = dict(eps=0.05, n_w=1.3)


def efficiencies(G, u_top):
    Nx = u_top.size
    uh = np.fft.fft(u_top) / Nx
    p = G.wavenumbers(Nx)
    e = np.real(G.gamma_p(Nx, 'u') / G.gamma_u) * np.abs(uh) ** 2
    return np.array([e[p == q][0] for q in (-1, 0, 1)])


def forward_hops(eps, n_w, Nx=32):
    out = []
    for om in OMEGAS:
        G = Grating2D(eps=eps, n_w=n_w, omega=om)
        out.append(efficiencies(G, hops_point(G, N=16, Nx=Nx, Nz=32, summation='pade', fields=False)['ubar']))
    return np.concatenate(out)


def forward_lsq(eps, n_w):
    out = []
    for om in OMEGAS:
        G = Grating2D(eps=eps, n_w=n_w, omega=om)
        P = make(G, activation='sin', **CFG['cosx']).solve()      # ~1e-14 accurate, ~1 s
        x = 2 * np.pi * np.arange(64) / 64
        out.append(efficiencies(G, P.evaluate('u', x, np.full(64, G.a))))
    return np.concatenate(out)


def invert(forward, data, x0):
    t0 = time.time()
    n = [0]

    def res(p):
        n[0] += 1
        return (forward(p[0], p[1]) - data) / np.maximum(np.abs(data), 1e-12)   # relative misfit
    sol = least_squares(res, x0, bounds=([0.0, 1.0], [0.3, 2.0]), x_scale=[0.05, 0.1], xtol=1e-14,
                        ftol=1e-14, gtol=1e-14, diff_step=1e-6)
    return dict(eps=float(sol.x[0]), n_w=float(sol.x[1]), cost=float(sol.cost), n_forward=n[0],
                time_s=time.time() - t0)


def main():
    matplotlib.use('Agg')
    os.makedirs(OUT, exist_ok=True)
    clean = forward_hops(TRUE['eps'], TRUE['n_w'], Nx=64)
    rng = np.random.default_rng(1)
    rows = []
    x0 = [START['eps'], START['n_w']]
    for noise in (0.0, 0.01):
        data = clean * (1 + noise * rng.standard_normal(clean.size))
        for name, fwd in (('HOPS/AWE', forward_hops), ('LSQ-PINN', forward_lsq)):
            r = invert(fwd, data, x0)
            r.update(model=name, noise=noise, eps_err=abs(r['eps'] - TRUE['eps']), n_w_err=abs(r['n_w'] - TRUE['n_w']))
            rows.append(r)
            print(r, flush=True)
    json.dump(dict(true=TRUE, start=START, omegas=OMEGAS.tolist(), data=clean.tolist(), results=rows),
              open(os.path.join(OUT, 'inverse.json'), 'w'), indent=1)
    # misfit landscape (HOPS, noise-free) with the recovered points
    E = np.linspace(0.02, 0.22, 21)
    Nw = np.linspace(1.02, 1.3, 21)
    J = np.zeros((E.size, Nw.size))
    for i, e in enumerate(E):
        for j, nw in enumerate(Nw):
            J[i, j] = np.sum(((forward_hops(e, nw) - clean) / clean) ** 2)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    m = ax.contourf(Nw, E, np.log10(J + 1e-30), 30, cmap='viridis')
    fig.colorbar(m, ax=ax, label='log10 misfit')
    ax.plot(TRUE['n_w'], TRUE['eps'], 'r*', ms=14, label='true')
    ax.plot(START['n_w'], START['eps'], 'wo', label='start')
    for r, mk in zip(rows, ('s', 'x', 'D', '+')):
        ax.plot(r['n_w'], r['eps'], mk, color='w' if r['noise'] == 0 else 'orange', ms=9,
                label=f"{r['model']}, noise {r['noise']:.0%}", mfc='none')
    ax.set_xlabel('$n^w$')
    ax.set_ylabel(r'$\varepsilon$')
    ax.set_title('inverse problem: recover $(\\varepsilon, n^w)$ from 15 efficiencies')
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'inverse.png'), dpi=110)
    return rows


if __name__ == '__main__':
    main()
