"""Hybrid PINN training: gradient-trained hidden layers + least-squares output layer.

The last layer of pinn.FieldNet is Linear(width, 2) -> (Re, Im).  With real weights W (2 x width) and
bias b this is the complex combination  u = sum_j (W_0j + i W_1j) h_j + (b_0 + i b_1),  i.e. a complex
coefficient vector c over the real last-hidden-layer features h_j and a constant.  Because problem (6) is
linear, the PINN loss is a linear least-squares problem in c for fixed hidden layers, so:

  * `lsq_refine(P)`   : after ordinary Adam/L-BFGS training, replace the output layer by its exact
                        least-squares optimum (one solve, ~1 s);
  * `lsgd_train(P)`   : the LSGD / variable-projection scheme of Cyr, Gulian, Patel, Perego & Trask (2020):
                        alternate (i) an exact least-squares solve for the output layer and (ii) Adam steps on
                        the hidden layers with the output layer frozen.  Each Adam step then optimises the
                        reduced (projected) loss, which removes the worst-conditioned directions of the
                        PINN landscape.

The least-squares rows are those of lsq_pinn.LSQPinn (identical weighting to the PINN loss), with the
network's hidden features in place of random features.
"""
import time

import numpy as np
import torch

from .lsq_pinn import LSQPinn


class NetFeatures:
    """Last-hidden-layer features (and x/z derivatives) of a FieldNet, times eps (flat ansatz scaling)."""

    def __init__(self, net, eps):
        self.net, self.eps = net, eps
        self.n = net.mlp[-1].in_features + 1

    @torch.no_grad()
    def __call__(self, x, z, order=2):
        net = self.net
        last = net.mlp[-1]
        # run derivs() with an identity output layer to get the hidden features
        W, b = last.weight.data.clone(), last.bias.data.clone()
        width = last.in_features
        tmp = torch.nn.Linear(width, width)
        tmp.weight.data = torch.eye(width)
        tmp.bias.data.zero_()
        net.mlp[-1] = tmp
        try:
            D = net.derivs(torch.as_tensor(x), torch.as_tensor(z), (), order=max(order, 1))
        finally:
            net.mlp[-1] = last
            last.weight.data, last.bias.data = W, b
        out = {}
        for k, v in D.items():
            v = v.numpy()
            const = np.ones((v.shape[0], 1)) if k == 'f' else np.zeros((v.shape[0], 1))
            out[k] = self.eps * np.hstack([v, const])
        return out


def _set_output(net, c):
    last = net.mlp[-1]
    last.weight.data = torch.tensor(np.vstack([c[:-1].real, c[:-1].imag]))
    last.bias.data = torch.tensor([c[-1].real, c[-1].imag])


def lsq_refine(P, **kw):
    """Exact least-squares output layer for a trained GratingPINN P (fixed eps, omega, flat ansatz)."""
    G = P.G
    fe = (NetFeatures(P.net_u, G.eps), NetFeatures(P.net_w, G.eps))
    L = LSQPinn(G, K=P.net_u.K, features=fe, rcond=kw.pop('rcond', 1e-12), **kw).solve()
    _set_output(P.net_u, L.c_u)
    _set_output(P.net_w, L.c_w)
    return L


def lsgd_train(P, outer=40, inner=50, lr=2e-3, lbfgs_iters=0, verbose=True, **kw):
    """LSGD (Cyr et al. 2020): alternate exact LSQ output layer and Adam on the hidden layers."""
    t0 = time.time()
    hidden = [p for n_, p in list(P.net_u.named_parameters()) + list(P.net_w.named_parameters())
              if not n_.startswith(f'mlp.{len(P.net_u.mlp) - 1}.')]
    opt = torch.optim.Adam(hidden, lr=lr)
    sched = torch.optim.lr_scheduler.ExponentialLR(opt, gamma=0.05 ** (1.0 / max(1, outer * inner)))
    for o in range(outer):
        L = lsq_refine(P, **kw)
        P.sample()
        for _ in range(inner):
            opt.zero_grad()
            loss, parts = P.loss()
            loss.backward()
            opt.step()
            sched.step()
        P.history.append(('lsgd', o, L.loss, parts, time.time() - t0))
        if verbose:
            print(f'  lsgd {o:3d}  lsq-loss {L.loss:.3e}  adam-loss {float(loss):.3e}  ({time.time() - t0:.0f} s)',
                  flush=True)
    L = lsq_refine(P, **kw)
    P.train_time = time.time() - t0
    return L
