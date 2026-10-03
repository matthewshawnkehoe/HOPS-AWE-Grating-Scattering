"""timing_benchmark.py -- cost per step of each solver, single thread (torch and BLAS), same machine.

  HOPS/AWE point solve (delta = 0, N = 16, Pade)      one solve
  LSQ interface PINN (compare_lsq.CFG)                 one assemble + least-squares solve
  gradient PINN (pinn.py, Taylor-mode derivatives)     one Adam step (2000 interior points)
  DeepXDE PINN (deepxde_solver.py, reverse-mode AD)    one Adam step (1000 interior points)
Writes results/variants/timing.json.
"""
import json
import os
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('DDE_BACKEND', 'pytorch')
import numpy as np
import torch

torch.set_num_threads(1)
from pinn_hops import Grating2D, GratingPINN
from pinn_hops.hops_reference import hops_point
from compare_lsq import make, CFG, ACT

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'variants')


def tm(f, n):
    f()
    t0 = time.time()
    for _ in range(n):
        f()
    return (time.time() - t0) / n


def main():
    rows = {}
    for name, G, K in (('dielectric', Grating2D(eps=0.1), 6),
                       ('gold', Grating2D(eps=0.1, n_w=1.48 + 1.883j, profile='cos4x'), 12)):
        Nx = 32 if G.profile == 'cosx' else 128
        r = dict(hops_s=tm(lambda: hops_point(G, N=16, Nx=Nx, Nz=48, summation='pade', fields=False), 3),
                 lsq_s=tm(lambda: make(G, activation=ACT[G.profile], **CFG[G.profile]).solve(), 2))
        P = GratingPINN(G, K=K, width=48, depth=4)
        opt = torch.optim.Adam(P.parameters(), lr=1e-3)
        P.sample()

        def step():
            opt.zero_grad()
            L, _ = P.loss()
            L.backward()
            opt.step()
        r['pinn_adam_step_s'] = tm(step, 10)
        import deepxde_solver as ds
        m, w = ds.build(G, K=K, n_domain=1000, n_if=128)
        m.compile('adam', lr=1e-3, loss_weights=w)
        m.train(iterations=1, display_every=10 ** 6, verbose=0)
        t0 = time.time()
        m.train(iterations=10, display_every=10 ** 6, verbose=0)
        r['deepxde_adam_step_s'] = (time.time() - t0) / 10
        rows[name] = r
        print(name, r, flush=True)
    json.dump(rows, open(os.path.join(OUT, 'timing.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
