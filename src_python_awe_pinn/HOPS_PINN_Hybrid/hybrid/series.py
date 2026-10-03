"""HOPS eps-series for ONE interface shape at ONE frequency (delta = 0), for arbitrary sampled profiles,
and the three forward maps used by the inverse problem (Kaplan & Nicholls, Appl. Numer. Math. 2019):

  forward_taylor   HOPS (TFE, N orders) summed by Taylor at the actual size -- the paper's forward model
  forward_pade     same series, Pade-summed
  forward_hybrid   same series, physics-informed summation (least-squares PINN over the HOPS basis,
                   hybrid.core.PISum), optionally enriched with random features

The profile is given by its samples g(x_j) on the uniform grid; it is split as g = eps_s * f with
eps_s = max|g| (so |f| <= 1), f_x and f_xx by FFT.  The data is the far-field (near-field) pattern
u_a(x) = u(x, a), cf. the paper's L[U].
"""
import time
import warnings

import numpy as np

from .core import spectral_dx, PISum, TWO_PI, Grating2D, RandomFeatures     # noqa: F401  (sets sys.path)
from hops import cheb, setup_2d, setup_zeta_psi_n_m, padesum        # noqa: E402
from hops.coupled import two_layer_solve_coupled                     # noqa: E402
from hops import field_tfe_helmholtz_m_and_n, field_tfe_helmholtz_m_and_n_lf  # noqa: E402
from pinn_hops.problem import outgoing_sqrt                          # noqa: E402


class ShapeSeries:
    """eps-series (M = 0) of HOPS for a sampled profile g; exposes the attributes PISum needs."""

    def __init__(self, g, n_u, n_w, omega, alpha=0.0, a=1.0, b=1.0, mode='TE', N=10, Nz=None, fields=True):
        g = np.asarray(g, float)
        self.Nx = Nx = g.size
        self.Nz = Nz = Nz or Nx
        self.N, self.M = N, 0
        self.a, self.b, self.mode = a, b, mode
        self.n_u, self.n_w = complex(n_u), complex(n_w)
        self.omega_bar = omega
        self.alpha_bar = alpha
        self.eps = float(np.max(np.abs(g))) or 1.0
        self.x = TWO_PI * np.arange(Nx) / Nx
        self.f = g / self.eps
        self.fx = np.real(spectral_dx(self.f.astype(complex)))
        self.fxx = np.real(spectral_dx(self.f.astype(complex), 2))
        D, t = cheb(Nz)
        self.zeta = 0.5 * (t + 1)
        self.Dzeta = 2 * D
        self.tau2 = 1.0 + 0j if mode == 'TE' else (self.n_u / self.n_w) ** 2
        gu = outgoing_sqrt(self.n_u ** 2 * omega ** 2 - alpha ** 2)
        gw = outgoing_sqrt(self.n_w ** 2 * omega ** 2 - alpha ** 2)
        xx, pp, abp, gup, _, _ = setup_2d(Nx, TWO_PI, alpha, gu)
        _, _, _, gwp, _, _ = setup_2d(Nx, TWO_PI, alpha, gw)
        M2 = 2                                   # solver needs M >= 1; only m = 0 is used (delta = 0)
        zeta, psi = setup_zeta_psi_n_m(xx, pp, alpha, gu, self.f, self.fx, Nx, N, M2)
        I = np.eye(Nz + 1)
        U, W, ub, wb = two_layer_solve_coupled(self.tau2, zeta, psi, gup, gwp, N, Nx, self.f, self.fx, pp, alpha,
                                               gu, gw, D, a, b, Nz, M2, I, abp)
        self.ubar_n = ub[:, 0, :]                # (Nx, N+1)
        self.wbar_n = wb[:, 0, :]
        self.nm = [(n, 0) for n in range(N + 1)]
        if fields:
            un = field_tfe_helmholtz_m_and_n(U, self.f, pp, gup, alpha, gu, D, a, Nx, Nz, N, M2, I, abp)
            wn = field_tfe_helmholtz_m_and_n_lf(W, self.f, pp, gwp, alpha, gw, D, b, Nx, Nz, N, M2, I, abp)
            self.basis = {'u': self._derivs(un[:, :, 0, :]), 'w': self._derivs(wn[:, :, 0, :])}

    def _derivs(self, B):
        Bx = spectral_dx(B)
        Dz = self.Dzeta
        return dict(f=B, x=Bx, z=np.einsum('lk,jkq->jlq', Dz, B), xx=spectral_dx(B, 2),
                    zz=np.einsum('lk,jkq->jlq', Dz @ Dz, B), xz=np.einsum('lk,jkq->jlq', Dz, Bx))

    def taylor_weights(self, eps, delta=0.0):
        return np.array([eps ** n for (n, _) in self.nm], dtype=complex)

    def grating(self, eps, omega):
        return Grating2D(n_u=self.n_u, n_w=self.n_w, omega=omega, eps=eps, alpha=self.alpha_bar, a=self.a,
                         b=self.b, mode=self.mode)


def forward_taylor(g, phys, N=10, Nz=None):
    S = ShapeSeries(g, N=N, Nz=Nz, fields=False, **phys)
    return np.polynomial.polynomial.polyval(S.eps, S.ubar_n.T)


def forward_pade(g, phys, N=10, Nz=None):
    S = ShapeSeries(g, N=N, Nz=Nz, fields=False, **phys)
    out = np.empty(S.Nx, complex)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for j in range(S.Nx):
            out[j] = padesum(S.ubar_n[j], S.eps, N // 2)[0]
    return out


def forward_hybrid(g, phys, N=10, Nz=None, enrich=None, return_all=False, ridge=None):
    S = ShapeSeries(g, N=N, Nz=Nz, fields=True, **phys)
    P = PISum(S, enrich=RandomFeatures(S, **enrich) if enrich else None, ridge=ridge)
    s = P.solve(S.eps, S.omega_bar)
    nu = P.nfeat['u']
    ua = P.feat['u']['f'][:, 0, :] @ s['c'][:nu]
    return (ua, s, S) if return_all else ua
