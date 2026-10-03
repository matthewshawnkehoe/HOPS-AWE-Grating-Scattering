"""study_activations.py -- activation functions and the I-PINN choice for the gradient-trained PINN.

Same PINN as compare_point.py (flat ansatz, exact FFT TBC, Adam + L-BFGS), only the activation changes:
  tanh             (default of pinn.py)
  sin              (SIREN-like; periodic in the pre-activation)
  laaf             (layer-wise locally adaptive tanh, Jagtap, Kawaguchi & Karniadakis 2020)
  ipinn_tanh_sin   (I-PINN, Sarma et al. CMAME 2024: a different activation in each subdomain,
                    tanh above the interface, sin below)
Each run is a separate process-safe unit: results/variants/act_<case>_<variant>.json (skipped if present).

    python study_activations.py --cases dielectric gold --budget quick
"""
import argparse
import json
import os
import time

import numpy as np

from pinn_hops import Grating2D, GratingPINN
from pinn_hops.hops_reference import hops_point
from pinn_hops.hybrid import lsq_refine, lsgd_train
from compare_point import CASES, BUDGETS

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'variants')
VARIANTS = {'tanh': 'tanh', 'sin': 'sin', 'laaf': 'laaf', 'ipinn_tanh_sin': ('tanh', 'sin'),
            'tanh+lsq': 'tanh',        # the tanh run, then an exact least-squares output layer (hybrid.lsq_refine)
            'lsgd': 'tanh'}            # LSGD / variable projection training (hybrid.lsgd_train), same Adam budget


def run(case, variant, budget='quick', seed=0):
    os.makedirs(OUT, exist_ok=True)
    fn = os.path.join(OUT, f'act_{case}_{variant}.json')
    if os.path.exists(fn):
        return json.load(open(fn))
    cfg = CASES[case]
    G = Grating2D(**cfg['G'])
    ref = hops_point(G, N=16, Nx=128, Nz=48, summation='pade', fields=True)
    P = GratingPINN(G, width=48, depth=4, seed=seed, activation=VARIANTS[variant], **cfg['net'])
    t0 = time.time()
    if variant == 'tanh+lsq':
        P.load(os.path.join(OUT, f'act_{case}_tanh.pt'))
        prev = json.load(open(os.path.join(OUT, f'act_{case}_tanh.json')))
        lsq_refine(P)
        P.train_time = prev['t_train'] + (time.time() - t0)
    elif variant == 'lsgd':
        b = BUDGETS[budget]
        lsgd_train(P, outer=b['adam_iters'] // 50, inner=50)
    else:
        P.train(**BUDGETS[budget], verbose=True, log_every=500)
    R, T, D = P.energy()
    x = ref['x']
    U = P.evaluate('u', x, G.g(x))
    P.n_int, P.n_if = 4000, 512
    P.sample()
    L, _ = P.loss()
    row = dict(case=case, variant=variant, budget=budget, R_hops=ref['R'], R=R, R_relerr=abs(R - ref['R']) / ref['R'],
               D_hops=ref['D'], D=D, D_absdiff=abs(D - ref['D']),
               U_relerr=float(np.abs(U - ref['U']).max() / np.abs(ref['U']).max()),
               loss=float(L), t_train=P.train_time)
    json.dump(row, open(fn, 'w'), indent=1)
    P.save(fn.replace('.json', '.pt'))
    print(row, flush=True)
    return row


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--cases', nargs='*', default=['dielectric', 'gold'])
    ap.add_argument('--variants', nargs='*', default=list(VARIANTS))
    ap.add_argument('--budget', default='quick')
    a = ap.parse_args()
    for c in a.cases:
        for v in a.variants:
            run(c, v, a.budget)
