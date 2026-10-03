"""The 2D two-layer grating scattering problem of Kehoe & Nicholls (J. Sci. Comput. 2024), eqs. (6a)-(6h).

Unknowns (Bloch phase e^{i alpha x} factored out, both d-periodic, d = 2 pi):
    u(x, z)  scattered (reflected) field in  S^u = {g(x) < z < a}
    w(x, z)  transmitted field in            S^w = {-b < z < g(x)}
Equations
    (6a)  Delta u + 2 i alpha u_x + (gamma^u)^2 u = 0                    g < z < a
    (6b)  Delta w + 2 i alpha w_x + (gamma^w)^2 w = 0                    -b < z < g
    (6c)  u - w = zeta,        zeta = -exp(-i gamma^u g)                 z = g
    (6d)  [d_N u - i alpha g_x u] - tau^2 [d_N w - i alpha g_x w] = psi,  N = (-g_x, 1),
          psi = (i gamma^u + i alpha g_x) exp(-i gamma^u g)               z = g
    (6e)  d_z u - T^u[u] = 0,  T^u = Fourier multiplier  i gamma^u_p      z = a
    (6f)  d_z w - T^w[w] = 0,  T^w = Fourier multiplier -i gamma^w_p      z = -b
    (6g,h) u, w periodic in x   (built into the network: Fourier features in x)
with alpha = k^u sin(theta), gamma^q = sqrt((k^q)^2 - alpha^2), k^q = n^q omega / c0, tau^2 = 1 (TE) or
(n^u/n^w)^2 (TM), gamma^q_p = sqrt((k^q)^2 - (alpha + p)^2) with Im >= 0, and g(x) = eps f(x).

Reflectivity map and energy (Sect. 2.2): with u(x, a) = sum_p xi_p e^{ipx},
    R = sum_{p in U^u} (gamma^u_p / gamma^u) |xi_p|^2,   T = tau^2 sum_{p in U^w} (gamma^w_p / gamma^u) |w^_p(-b)|^2,
    D = 1 - R - T      (propagating sets: Re (alpha+p)^2 < Re k^2, the rule of energy_defect.m).
This module is independent of the HOPS code.
"""
from dataclasses import dataclass
import numpy as np


PROFILES = {
    'cosx': (lambda x: np.cos(x), lambda x: -np.sin(x)),
    'cos4x': (lambda x: np.cos(4 * x), lambda x: -4 * np.sin(4 * x)),
    'cos4x_over4': (lambda x: 0.25 * np.cos(4 * x), lambda x: -np.sin(4 * x)),
    'sinx': (lambda x: np.sin(x), lambda x: np.cos(x)),
    'cos2x': (lambda x: np.cos(2 * x), lambda x: -2 * np.sin(2 * x)),
    'cos3x': (lambda x: np.cos(3 * x), lambda x: -3 * np.sin(3 * x)),
    'sin3x': (lambda x: np.sin(3 * x), lambda x: 3 * np.cos(3 * x)),
    'sin4x': (lambda x: np.sin(4 * x), lambda x: 4 * np.cos(4 * x)),
    'sin5x': (lambda x: np.sin(5 * x), lambda x: 5 * np.cos(5 * x)),
}


def outgoing_sqrt(v):
    """sqrt with Im >= 0 (outgoing / decaying branch; setup_2d.m convention)."""
    v = np.asarray(v, dtype=complex)
    th = np.angle(v)
    th = np.where(th < 0, th + 2 * np.pi, th)
    return np.sqrt(np.abs(v)) * np.exp(0.5j * th)


@dataclass
class Grating2D:
    n_u: complex = 1.0
    n_w: complex = 1.1
    omega: float = 1.5
    eps: float = 0.1
    profile: str = 'cosx'
    alpha: float = 0.0            # lateral Bloch wavenumber (0 = normal incidence, as the paper's examples)
    a: float = 1.0
    b: float = 1.0
    mode: str = 'TM'
    c0: float = 1.0

    @property
    def k_u(self):
        return self.n_u * self.omega / self.c0

    @property
    def k_w(self):
        return self.n_w * self.omega / self.c0

    @property
    def gamma_u(self):
        return outgoing_sqrt(self.k_u ** 2 - self.alpha ** 2)

    @property
    def gamma_w(self):
        return outgoing_sqrt(self.k_w ** 2 - self.alpha ** 2)

    @property
    def tau2(self):
        return 1.0 + 0j if self.mode == 'TE' else (self.n_u / self.n_w) ** 2

    def f(self, x):
        return PROFILES[self.profile][0](x)

    def f_x(self, x):
        return PROFILES[self.profile][1](x)

    def g(self, x):
        return self.eps * self.f(x)

    def g_x(self, x):
        return self.eps * self.f_x(x)

    def zeta(self, x):
        return -np.exp(-1j * self.gamma_u * self.g(x))

    def psi(self, x):
        return (1j * self.gamma_u + 1j * self.alpha * self.g_x(x)) * np.exp(-1j * self.gamma_u * self.g(x))

    def wavenumbers(self, Nx):
        return np.fft.fftfreq(Nx, 1.0 / Nx)                      # integers p (d = 2 pi)

    def gamma_p(self, Nx, layer):
        k = self.k_u if layer == 'u' else self.k_w
        return outgoing_sqrt(k ** 2 - (self.alpha + self.wavenumbers(Nx)) ** 2)

    def propagating(self, Nx, layer):
        k = self.k_u if layer == 'u' else self.k_w
        return np.real((self.alpha + self.wavenumbers(Nx)) ** 2) < np.real(k ** 2)

    def energy(self, u_top, w_bot):
        """R, T, D from the traces u(x_j, a), w(x_j, -b) on a uniform grid x_j = 2 pi j / Nx."""
        Nx = u_top.size
        uh = np.fft.fft(u_top) / Nx
        wh = np.fft.fft(w_bot) / Nx
        gu, gw = self.gamma_p(Nx, 'u'), self.gamma_p(Nx, 'w')
        g0 = self.gamma_u
        R = np.sum(np.where(self.propagating(Nx, 'u'), gu / g0 * np.abs(uh) ** 2, 0.0))
        T = self.tau2 * np.sum(np.where(self.propagating(Nx, 'w'), gw / g0 * np.abs(wh) ** 2, 0.0))
        return np.real(R), np.real(T), np.real(1.0 - R - T)

    def flat_solution(self):
        """Exact solution for eps = 0 (flat interface): u = r e^{i gamma^u z}, w = t e^{-i gamma^w z}."""
        gu, gw, t2 = self.gamma_u, self.gamma_w, self.tau2
        r = (gu - t2 * gw) / (gu + t2 * gw)
        t = 2 * gu / (gu + t2 * gw)
        return r, t
