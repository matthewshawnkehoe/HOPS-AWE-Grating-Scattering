"""deepxde_solver.py -- the same grating problem (6a)-(6h) solved with the DeepXDE library (Lu et al. 2021).

This is a standard DeepXDE PINN, built from DeepXDE's own components:
  * geometry: dde.geometry.Rectangle [0, 2 pi] x [-b, a]; interior points from DeepXDE's sampler
    (uniform in the rectangle, resampled by dde.callbacks.PDEPointResampler);
  * network: dde.nn.PFNN -- parallel sub-networks, one per output (Re u, Im u, Re w, Im w): the
    "one network per subdomain" choice of XPINN / I-PINN, made with DeepXDE's own class;
    activation 'tanh', 'sin' or DeepXDE's adaptive 'LAAF-10 tanh';
  * exact periodicity through net.apply_feature_transform (cos kx, sin kx, z);
  * flat-interface ansatz through net.apply_output_transform (u = u_flat + eps N_u, ...);
  * PDE (6a, 6b) with dde.grad.hessian / jacobian, masked to the layer each point lies in
    (z > g(x): Helmholtz for u; z < g(x): for w);
  * interface (6c, 6d) and the exact FFT-based DtN conditions (6e, 6f) as dde.icbc.PointSetOperatorBC
    on fixed point sets (the uniform boundary grid makes the FFT possible);
  * Adam, then DeepXDE's L-BFGS; float64.

    python deepxde_solver.py --case dielectric --iters 1000 --lbfgs 300
    python deepxde_solver.py --case gold --activation "LAAF-10 tanh"
Writes results/variants/deepxde_<case>_<activation>.json
"""
import argparse
import json
import os
import time

os.environ.setdefault('DDE_BACKEND', 'pytorch')
import numpy as np
import deepxde as dde
import torch

from pinn_hops import Grating2D
from pinn_hops.hops_reference import hops_point

dde.config.set_default_float('float64')
torch.set_num_threads(int(os.environ.get('TORCH_THREADS', '1')))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results', 'variants')
CASES = {'dielectric': (dict(n_w=1.1, profile='cosx', eps=0.1, omega=1.5), 6),
         'gold': (dict(n_w=1.48 + 1.883j, profile='cos4x', eps=0.1, omega=1.5), 12),
         'silver': (dict(n_w=0.05 + 2.275j, profile='cos4x', eps=0.1, omega=1.5), 12)}


def build(G, K=6, width=48, depth=4, activation='tanh', n_domain=2000, n_if=256, Nx_bc=64, lam=10.0):
    gu, gw, t2, al = complex(G.gamma_u), complex(G.gamma_w), complex(G.tau2), G.alpha
    k2 = {'u': complex(G.k_u) ** 2, 'w': complex(G.k_w) ** 2}
    fname = {'cosx': (torch.cos, lambda x: -torch.sin(x)), 'cos4x': (lambda x: torch.cos(4 * x),
                                                                        lambda x: -4 * torch.sin(4 * x))}[G.profile]
    g_t = lambda x: G.eps * fname[0](x)
    gx_t = lambda x: G.eps * fname[1](x)
    r, t = G.flat_solution()

    geom = dde.geometry.Rectangle([0, -G.b], [2 * np.pi, G.a])

    # ---- PDE (6a), (6b), masked by layer
    def pde(X, Y):
        x, z = X[:, 0:1], X[:, 1:2]
        up = (z > g_t(x)).to(X.dtype)
        res = []
        for off, layer, mask in ((0, 'u', up), (2, 'w', 1 - up)):
            c = []
            for j in (0, 1):
                c.append(dde.grad.hessian(Y, X, component=off + j, i=0, j=0) +
                         dde.grad.hessian(Y, X, component=off + j, i=1, j=1))
            lap = torch.complex(c[0], c[1])
            ux = torch.complex(dde.grad.jacobian(Y, X, i=off, j=0), dde.grad.jacobian(Y, X, i=off + 1, j=0))
            u = torch.complex(Y[:, off:off + 1], Y[:, off + 1:off + 2])
            rr = (lap + 2j * al * ux + (k2[layer] - al ** 2) * u) / (1 + abs(k2[layer])) * mask
            res += [rr.real, rr.imag]
        return res

    # ---- boundary / interface point sets (order fixes their offsets in the batch)
    xi = np.linspace(0, 2 * np.pi, n_if, endpoint=False) + np.pi / n_if
    P_if = np.stack([xi, G.g(xi)], 1)
    xb = 2 * np.pi * np.arange(Nx_bc) / Nx_bc
    P_top = np.stack([xb, np.full(Nx_bc, G.a)], 1)
    P_bot = np.stack([xb, np.full(Nx_bc, -G.b)], 1)

    def cplx(Y, off):
        return torch.complex(Y[:, off], Y[:, off + 1])

    def dz(Y, X, off):
        return torch.complex(dde.grad.jacobian(Y, X, i=off, j=1)[:, 0], dde.grad.jacobian(Y, X, i=off + 1, j=1)[:, 0])

    def dxx(Y, X, off):
        return torch.complex(dde.grad.jacobian(Y, X, i=off, j=0)[:, 0], dde.grad.jacobian(Y, X, i=off + 1, j=0)[:, 0])

    def make(name, part, b, n):
        """PointSetOperatorBC function.  DeepXDE puts the BC point sets first in the batch, in the order
        the BCs are listed, so this BC's points are rows b .. b+n-1 of X (needed for the FFT)."""
        sl = slice(b, b + n)

        def func(X, Y, Xnp):
            x, z = X[sl, 0], X[sl, 1]
            if name in ('dir', 'neu'):
                gx = gx_t(x)
                inc = torch.exp(-1j * gu * z)
                u, w = cplx(Y, 0)[sl], cplx(Y, 2)[sl]
                if name == 'dir':                                    # (6c): u - w - zeta, zeta = -inc
                    v = u - w + inc
                else:                                                # (6d)
                    nu = dz(Y, X, 0)[sl] - gx * dxx(Y, X, 0)[sl] - 1j * al * gx * u
                    nw = dz(Y, X, 2)[sl] - gx * dxx(Y, X, 2)[sl] - 1j * al * gx * w
                    v = (nu - t2 * nw - (1j * gu + 1j * al * gx) * inc) / (1 + abs(gu))
            else:                                                    # (6e), (6f): exact DtN by FFT
                off, sgn, layer = (0, 1, 'u') if name == 'top' else (2, -1, 'w')
                f = cplx(Y, off)[sl]
                gp = torch.tensor(G.gamma_p(n, layer))
                v = (dz(Y, X, off)[sl] - torch.fft.ifft(sgn * 1j * gp * torch.fft.fft(f))) / (1 + abs(k2[layer]) ** 0.5)
            v = v.real if part == 0 else v.imag
            full = torch.zeros(X.shape[0], 1, dtype=X.dtype)
            full[sl, 0] = v
            return full
        return func

    bcs, beg = [], 0
    for name, pts in (('dir', P_if), ('neu', P_if), ('top', P_top), ('bot', P_bot)):
        for part in (0, 1):
            bcs.append(dde.icbc.PointSetOperatorBC(pts, np.zeros((len(pts), 1)), make(name, part, beg, len(pts))))
            beg += len(pts)

    data = dde.data.PDE(geom, pde, bcs, num_domain=n_domain, num_boundary=0, num_test=500)
    kk = torch.arange(1, K + 1, dtype=torch.float64)
    net = dde.nn.PFNN([2 * K + 1] + [[width] * 4] * depth + [4], activation, 'Glorot normal')

    def feat(X):
        kx = X[:, 0:1] * kk
        return torch.cat([torch.cos(kx), torch.sin(kx), X[:, 1:2]], 1)

    def out(X, Y):
        z = X[:, 1:2]
        u0 = r * torch.exp(1j * gu * z)
        w0 = t * torch.exp(-1j * gw * z)
        e = G.eps
        return torch.cat([u0.real + e * Y[:, 0:1], u0.imag + e * Y[:, 1:2],
                          w0.real + e * Y[:, 2:3], w0.imag + e * Y[:, 3:4]], 1)
    net.apply_feature_transform(feat)
    net.apply_output_transform(out)
    model = dde.Model(data, net)
    weights = [1.0] * 4 + [lam] * 8
    return model, weights


def run(case='dielectric', activation='tanh', iters=1000, lbfgs=300, seed=0, n_domain=1000):
    dde.config.set_random_seed(seed)
    gk, K = CASES[case]
    G = Grating2D(**gk)
    model, weights = build(G, K=K, activation=activation, n_domain=n_domain, n_if=128)
    t0 = time.time()
    model.compile('adam', lr=2e-3, loss_weights=weights,
                  decay=('inverse time', max(1, iters // 4), 1.0))
    model.train(iterations=iters, display_every=100,
                callbacks=[dde.callbacks.PDEPointResampler(period=250)])
    if lbfgs:
        dde.optimizers.config.set_LBFGS_options(maxiter=lbfgs)
        model.compile('L-BFGS', loss_weights=weights)
        model.train(display_every=max(lbfgs // 2, 1))
    t_train = time.time() - t0
    # R, T, D from the traces, as for the other solvers
    Nx = 64
    x = 2 * np.pi * np.arange(Nx) / Nx
    Yt = model.predict(np.stack([x, np.full(Nx, G.a)], 1))
    Yb = model.predict(np.stack([x, np.full(Nx, -G.b)], 1))
    R, T, D = G.energy(Yt[:, 0] + 1j * Yt[:, 1], Yb[:, 2] + 1j * Yb[:, 3])
    ref = hops_point(G, N=16, Nx=128, Nz=48, summation='pade')
    Yi = model.predict(np.stack([ref['x'], G.g(ref['x'])], 1))
    U = Yi[:, 0] + 1j * Yi[:, 1]
    row = dict(case=case, variant='deepxde_' + activation.replace(' ', '_'), R_hops=float(ref['R']), R=float(R),
               R_relerr=float(abs(R - ref['R']) / ref['R']), D_hops=float(ref['D']), D=float(D),
               D_absdiff=float(abs(D - ref['D'])),
               U_relerr=float(np.abs(U - ref['U']).max() / np.abs(ref['U']).max()),
               loss=float(np.sum(model.losshistory.loss_train[-1])), t_train=t_train, iters=iters, lbfgs=lbfgs)
    print(row)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f"deepxde_{case}_{activation.replace(' ', '_')}.json"), 'w') as fh:
        json.dump(row, fh, indent=1)
    return row


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='dielectric', choices=list(CASES))
    ap.add_argument('--activation', default='tanh')
    ap.add_argument('--iters', type=int, default=1000)
    ap.add_argument('--lbfgs', type=int, default=300)
    ap.add_argument('--n-domain', type=int, default=1000)
    a = ap.parse_args()
    run(a.case, a.activation, a.iters, a.lbfgs, n_domain=a.n_domain)
