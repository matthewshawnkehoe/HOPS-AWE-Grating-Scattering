"""inverse_paper.py -- periodic grating identification by nonlinear least squares with HOPS
(Kaplan & Nicholls, Appl. Numer. Math. 143 (2019) 20-34), and where the PINN hybrid does / does not help.

Paper setup (Sect. 6): unknown g(x_j) at Nx = 32 equispaced points, data u_a(x_j) = u(x_j, a) (TE, b = 1),
Gauss-Newton or Levenberg-Marquardt with a first-order finite-difference Jacobian (Remark 5.1), initial
guess g = 0, stop when the update is below 1e-7, HOPS with N = 10 Taylor orders:
  (6.1) g = eps exp(cos 2x)                          alpha = 0,   k_u = 1.1, k_w = 5.5, a = 1   (6.2)
  (6.3) g = eps sech(2x)                             alpha = 0.2, k_u = 1.3, k_w = 6.8, a = 1   (6.5)
  (6.4) g = eps [tanh(2(x + 3pi/5)) - tanh(2(x - 3pi/5))]                                   (6.5)

Part A  forward accuracy at fixed frequency: u_a from the eps-series (N = 10, 16) summed by Taylor,
        Pade, or the physics-informed (hybrid) sum, vs a converged solve (N = 30, Nx = 96, Nz = 64).
Part B  the paper's experiment (data made with the same N = 10 HOPS model, as in the paper): GN and LM
        iterations and L-infinity errors, cf. the paper's Tables 1-12.
Part C  the same inversion with data from the converged model (no "inverse crime"): the recovery error
        is then set by the forward-model error, amplified by the ill-conditioning of u_a(g).

    python inverse_paper.py            # results/inverse/{forward.json, paper.json, crime_free.json, inverse.md}
"""
import json
import os
import time
import warnings

import numpy as np

from hybrid.series import forward_taylor, forward_pade, forward_hybrid

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'inverse')
P62 = dict(n_u=1.1, n_w=5.5, omega=1.0, alpha=0.0, a=1.0, b=1.0, mode='TE')
P65 = dict(n_u=1.3, n_w=6.8, omega=1.0, alpha=0.2, a=1.0, b=1.0, mode='TE')
PROFILES = {
    '6.1': (lambda x: np.exp(np.cos(2 * x)), P62),
    '6.3': (lambda x: 1 / np.cosh(2 * (x - np.pi)), P65),
    '6.4': (lambda x: np.tanh(2 * (x - np.pi + 3 * np.pi / 5)) - np.tanh(2 * (x - np.pi - 3 * np.pi / 5)), P65),
}
MODELS = {'taylor': forward_taylor, 'pade': forward_pade, 'hybrid': forward_hybrid}


def converged(prof, eps, phys):
    x96 = 2 * np.pi * np.arange(96) / 96
    u, s, _ = forward_hybrid(eps * prof(x96), phys, N=30, Nz=64, return_all=True)
    return u[::3], s['loss']                                  # on the 32-point grid


def invert(F, d, Nx=32, method='GN', tol=1e-7, max_it=20, h=1e-7, g_true=None):
    res = lambda u: np.concatenate([(u - d).real, (u - d).imag])
    g = np.zeros(Nx)
    lam = 1e-3
    r = res(F(g))
    hist = []
    t0 = time.time()
    for it in range(1, max_it + 1):
        J = np.empty((r.size, Nx))
        for j in range(Nx):
            e = np.zeros(Nx)
            e[j] = h
            J[:, j] = (res(F(g + e)) - r) / h
        if method == 'GN':
            v = np.linalg.lstsq(J, -r, rcond=None)[0]
            g, r = g + v, res(F(g + v))
        else:
            JTJ, JTr = J.T @ J, J.T @ r
            while True:
                v = np.linalg.solve(JTJ + lam * np.diag(np.diag(JTJ)), -JTr)
                rn = res(F(g + v))
                if rn @ rn <= r @ r or lam > 1e10:
                    g, r, lam = g + v, rn, max(lam / 3, 1e-14)
                    break
                lam *= 4
        hist.append((float(np.max(np.abs(v))), float(np.sqrt(r @ r)),
                     float(np.max(np.abs(g - g_true))) if g_true is not None else np.nan))
        if hist[-1][0] < tol:
            break
    return g, it, hist, time.time() - t0


def part_a():
    rows = []
    x64 = 2 * np.pi * np.arange(64) / 64
    for p in ('6.1', '6.3'):
        prof, phys = PROFILES[p]
        for eps in (0.05, 0.1, 0.2, 0.3):
            ref, loss = converged(prof, eps, phys)
            row = dict(profile=p, eps=eps, ref_loss=loss)
            for m, fn in MODELS.items():
                for N in (10, 16):
                    t0 = time.time()
                    u = fn(eps * prof(x64), phys, N=N)[::2]
                    row[f'{m}_N{N}'] = float(np.max(np.abs(u - ref)))
                    row[f'{m}_N{N}_t'] = time.time() - t0
            rows.append(row)
            print('A', {k: (f'{v:.1e}' if isinstance(v, float) else v) for k, v in row.items() if not k.endswith('_t')},
                  flush=True)
    return rows


def part_bc(crime=True):
    rows = []
    x = 2 * np.pi * np.arange(32) / 32
    for p, (prof, phys) in PROFILES.items():
        for eps in ((0.001, 0.01, 0.05, 0.1) if crime else (0.01, 0.05, 0.1)):
            gt = eps * prof(x)
            if crime:
                d = forward_taylor(gt, phys, N=10)
            else:
                d, _ = converged(prof, eps, phys)
            for method in ('GN', 'LM'):
                for m in (('taylor',) if crime else ('taylor', 'pade')):
                    F = lambda g, fn=MODELS[m]: fn(g, phys, N=10)
                    g, it, hist, t = invert(F, d, method=method, g_true=gt)
                    err = float(np.max(np.abs(g - gt)))
                    # the paper's count: iterations until the absolute L-inf error is below 1e-7
                    it_tol = next((k + 1 for k, h_ in enumerate(hist) if h_[2] < 1e-7), None)
                    rows.append(dict(profile=p, eps=eps, method=method, model=m, iterations=it, it_to_1e7=it_tol,
                                     abs_err=err, rel_err=err / float(np.max(np.abs(gt))), time_s=t,
                                     err_history=[h_[2] for h_ in hist]))
                    print('B' if crime else 'C', p, eps, method, m, f'it {it} (to 1e-7: {it_tol})  abs {err:.2e}  rel {rows[-1]["rel_err"]:.2e}',
                          flush=True)
    return rows


def main():
    warnings.filterwarnings('ignore')
    os.makedirs(OUT, exist_ok=True)
    parts = {}
    for name, fn in (('forward', part_a), ('paper', lambda: part_bc(True)), ('crime_free', lambda: part_bc(False))):
        f = os.path.join(OUT, name + '.json')
        if os.path.exists(f):
            parts[name] = json.load(open(f))
            continue
        parts[name] = fn()
        json.dump(parts[name], open(f, 'w'), indent=1)
    write_md(parts)


def write_md(parts):
    L = ['### A. Forward accuracy at fixed frequency (max |u_a - u_a,converged|, Nx = 64)', '',
         '| profile | eps | Taylor N=10 | Pade N=10 | hybrid N=10 | Taylor N=16 | Pade N=16 | hybrid N=16 |',
         '|---|---|---|---|---|---|---|---|']
    for r in parts['forward']:
        L.append(f"| ({r['profile']}) | {r['eps']} | " + ' | '.join(f"{r[f'{m}_N{N}']:.1e}" for N in (10, 16)
                                                           for m in ('taylor', 'pade', 'hybrid')) + ' |')
    L += ['', '### B. The paper\'s experiment (data from the same N = 10 model)', '',
          'iterations to reach absolute L-inf error < 1e-7 (the paper\'s count) / final absolute error', '',
          '| profile | eps | GN | LM |', '|---|---|---|---|']
    for p in PROFILES:
        for eps in sorted({r['eps'] for r in parts['paper'] if r['profile'] == p}):
            c = {r['method']: r for r in parts['paper'] if r['profile'] == p and r['eps'] == eps}
            L.append(f"| ({p}) | {eps} | {c['GN']['it_to_1e7'] or '-'} / {c['GN']['abs_err']:.1e} | "
                     f"{c['LM']['it_to_1e7'] or '-'} / {c['LM']['abs_err']:.1e} |")
    L += ['', '### C. Data from the converged model: relative L-inf error of the recovered g', '',
          '| profile | eps | GN Taylor | GN Pade | LM Taylor | LM Pade |', '|---|---|---|---|---|---|']
    for p in PROFILES:
        for eps in sorted({r['eps'] for r in parts['crime_free'] if r['profile'] == p}):
            c = {(r['method'], r['model']): r for r in parts['crime_free'] if r['profile'] == p and r['eps'] == eps}
            L.append(f"| ({p}) | {eps} | " + ' | '.join(f"{c[k]['rel_err']:.1e}" for k in
                                                        (('GN', 'taylor'), ('GN', 'pade'), ('LM', 'taylor'), ('LM', 'pade'))) + ' |')
    open(os.path.join(OUT, 'inverse.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
