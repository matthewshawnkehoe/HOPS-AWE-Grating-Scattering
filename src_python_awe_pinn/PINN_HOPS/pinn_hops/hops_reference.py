"""Bridge to the 2D HOPS/AWE code (the Python port in ../HOPS_Python, package `hops`) for comparisons.

Set the environment variable HOPS_PYTHON to the HOPS_Python directory if it is not a sibling of this
project.  Two kinds of references are provided:

* hops_point(G)      -- HOPS at ONE (eps, omega): the HOPS/AWE solver is run with omega_bar = omega
                        (delta = 0) and the eps-series is summed (Taylor or Pade).  Gives interface data
                        U, W, the traces at z = a, -b, the volume fields in physical coordinates and R, T, D.
* hops_refl_map(...) -- the actual reflectivity map / energy defect of refl_map.py (joint (eps, delta)
                        expansion, Taylor or Pade), for map/spectrum comparisons.
"""
import os
import sys
import warnings
import numpy as np

_HOPS = os.environ.get('HOPS_PYTHON', os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), 'HOPS_Python'))
if _HOPS not in sys.path:
    sys.path.insert(0, _HOPS)


def _hops():
    import hops
    return hops


def _sum_eps(c, eps, kind):
    """c (..., N+1): sum_n c_n eps^n by Taylor (all orders) or [N/2, N/2] Pade (padesum.m)."""
    h = _hops()
    if kind == 'taylor':
        return np.polynomial.polynomial.polyval(eps, np.moveaxis(c, -1, 0))
    N = c.shape[-1] - 1
    out = np.empty(c.shape[:-1], dtype=complex)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for idx in np.ndindex(out.shape):
            out[idx] = h.padesum(c[idx], eps, N // 2)[0]
    return out


def hops_point(G, N=16, Nx=32, Nz=32, summation='pade', fields=True, nz_out=33):
    """HOPS solution of the grating G (problem.Grating2D) at its (eps, omega), delta = 0."""
    h = _hops()
    from hops.coupled import two_layer_solve_coupled
    M = 2                                     # delta = 0: only the m = 0 coefficients are used
    d = 2 * np.pi
    x = (d / Nx) * np.arange(Nx)
    f, f_x = G.f(x), G.f_x(x)
    gu, gw = G.gamma_u, G.gamma_w
    xx, pp, abp, gup, _, _ = h.setup_2d(Nx, d, G.alpha, gu)
    _, _, _, gwp, _, _ = h.setup_2d(Nx, d, G.alpha, gw)
    zeta, psi = h.setup_zeta_psi_n_m(xx, pp, G.alpha, gu, f, f_x, Nx, N, M)
    Dz, _ = h.cheb(Nz)
    tau2 = G.tau2
    args = (tau2, zeta, psi, gup, gwp, N, Nx, f, f_x, pp, G.alpha, gu, gw, Dz, G.a, G.b, Nz, M, np.eye(Nz + 1), abp)
    U, W, ub, wb = two_layer_solve_coupled(*args)
    s = lambda A: _sum_eps(A[:, 0, :], G.eps, summation)            # m = 0 column, series in eps
    out = dict(x=x, U=s(U), W=s(W), ubar=s(ub), wbar=s(wb))
    out['R'], out['T'], out['D'] = G.energy(out['ubar'], out['wbar'])
    if fields:
        un = h.field_tfe_helmholtz_m_and_n(U, f, pp, gup, G.alpha, gu, Dz, G.a, Nx, Nz, N, M, np.eye(Nz + 1), abp)
        wn = h.field_tfe_helmholtz_m_and_n_lf(W, f, pp, gwp, G.alpha, gw, Dz, G.b, Nx, Nz, N, M, np.eye(Nz + 1), abp)
        tz = np.cos(np.pi * np.arange(Nz + 1) / Nz)
        g = G.eps * f[:, None]
        zu = (G.a / 2) * (tz - 1) + G.a                       # transformed coordinates z'
        zw = (G.b / 2) * (tz - 1)
        out['Zu'] = g + zu[None, :] * (G.a - g) / G.a           # physical z of the collocation nodes
        out['Zw'] = g + zw[None, :] * (G.b + g) / G.b
        out['u'] = _sum_eps(un[:, :, 0, :], G.eps, summation)     # (Nx, Nz+1) at the nodes (x, Zu)
        out['w'] = _sum_eps(wn[:, :, 0, :], G.eps, summation)
    return out


def hops_refl_map(scenario='dielectric', q=(1,), N_Eps=41, N_delta=41, **over):
    """refl_map.py results (list of windows with omega, lam, Eps, ru, rl, ee, RR)."""
    import refl_map as rm
    res, info = rm.run(scenario, qq=q, N_Eps=N_Eps, N_delta=N_delta, verbose=False, workers=1, **over)
    return res, info


class HOPSFieldTorch:
    """The HOPS solution as a differentiable torch function of PHYSICAL (x, z).

    The HOPS volume fields are known at the Fourier-Chebyshev nodes (x_j, z'_l) of the transformed
    coordinates; here they are turned into their spectral interpolant
        u(x, z) = sum_p sum_k C_pk e^{i p x} T_k(t),  t = 2 z'/a - 1,  z' = a (z - g)/(a - g)   (upper)
    (lower: t = 2 z'/b + 1, z' = b (z - g)/(b + g)), evaluated with torch so that autograd gives exact
    derivatives.  Plugging this into GratingPINN.residuals() evaluates the PINN loss of the HOPS solution:
    it must be ~0 if the PINN loss and the HOPS solver discretise the same equations (6a)-(6h)."""

    def __init__(self, ref, G):
        import torch
        self.torch = torch
        self.G = G
        Nx, Nz1 = ref['u'].shape
        Nz = Nz1 - 1
        tz = np.cos(np.pi * np.arange(Nz + 1) / Nz)
        Vinv = np.linalg.inv(np.polynomial.chebyshev.chebvander(tz, Nz))
        self.C = {}
        for layer in ('u', 'w'):
            F = np.fft.fft(ref[layer], axis=0) / Nx                   # (Nx, Nz+1) Fourier in x
            self.C[layer] = torch.tensor(F @ Vinv.T)                  # (Nx, Nz+1) Chebyshev coefficients
        self.p = torch.tensor(np.fft.fftfreq(Nx, 1.0 / Nx))
        self.Nz = Nz

    def _cheb(self, t):
        T = [self.torch.ones_like(t), t]
        for k in range(2, self.Nz + 1):
            T.append(2 * t * T[-1] - T[-2])
        return self.torch.stack(T, dim=1)                              # (n, Nz+1)

    def __call__(self, layer, x, z, eps, omega):
        torch = self.torch
        G = self.G
        f = {'cosx': torch.cos(x), 'cos4x': torch.cos(4 * x), 'cos4x_over4': 0.25 * torch.cos(4 * x),
             'sinx': torch.sin(x), 'cos2x': torch.cos(2 * x)}[G.profile]
        g = G.eps * f
        if layer == 'u':
            zp = G.a * (z - g) / (G.a - g)
            t = 2 * zp / G.a - 1
        else:
            zp = G.b * (z - g) / (G.b + g)
            t = 2 * zp / G.b + 1
        E = torch.exp(1j * x[:, None] * self.p[None, :])               # (n, Nx)
        Tk = self._cheb(t).to(torch.complex128)                        # (n, Nz+1)
        val = torch.sum((E @ self.C[layer]) * Tk, dim=1)
        return val.real, val.imag
