"""Least-squares interface PINN (random-feature / extreme-learning-machine PINN) for problem (6a)-(6h).

Why a different PINN
--------------------
Problem (6) is LINEAR in (u, w).  A network whose last layer is linear,
    u_theta(x, z) = u_flat(z) + sum_j c^u_j phi^u_j(x, z),     w_theta = w_flat + sum_j c^w_j phi^w_j,
has residuals of (6a)-(6f) that are affine in the output weights c.  The PINN loss (sum of squared residuals
at collocation points) is then a linear least-squares problem in c: its global minimiser comes from ONE
complex least-squares solve, instead of thousands of Adam / L-BFGS steps on a non-convex landscape.  This
is the physics-informed extreme learning machine (PIELM, Dwivedi & Srinivasan 2020) / random feature method
(RFM, Chen, Chi, E & Yang 2022) idea.  The hidden layer is random (RFM) or trained (see hybrid_lsq in
`GratingPINN`-based helpers below).

Everything else is exactly the problem that `pinn.py` and HOPS solve:
  * two subdomain networks (upper u, lower w) with separate features -- and, as in the I-PINN of
    Sarma et al. (CMAME 2024), a possibly different activation per subdomain;
  * exact periodicity: x enters only through cos(kx), sin(kx), k <= K;
  * the exact nonlocal DtN transparent conditions (6e, 6f) on a uniform x-grid (FFT, multiplier +-i gamma_p);
  * the flat-interface ansatz: u_flat, w_flat (Fresnel) are the eps = 0 solution; the features carry the
    correction.  (Here the ansatz only moves known terms to the right-hand side.)

Two hidden-layer families (`basis`):
  'rfm'        phi_j(x, z) = act(sum_k A_jk cos kx + B_jk sin kx + a_j zs + b_j), random A, B, a, b:
               a generic one-hidden-layer network in (x, z), the standard PIELM/RFM.
  'separable'  phi_{k,j}(x, z) = e^{ikx} act(a_j zs + b_j), |k| <= K: a Fourier-feature network whose
               hidden units are products of one x-feature and one z-feature (separable PINN, Cho et al. 2023).

The residual rows are nondimensionalised exactly as in pinn.py (PDE / (1 + |k^2|), fluxes / (1 + |gamma|),
weights lambda on (6c-f)) and scaled by 1/sqrt(#points) so the least-squares objective IS the PINN loss.
"""
import time
from dataclasses import replace

import numpy as np
import scipy.linalg as sla

TWO_PI = 2 * np.pi


def _act(name, h):
    """activation and its first two derivatives"""
    if name == 'tanh':
        t = np.tanh(h)
        d = 1 - t * t
        return t, d, -2 * t * d
    if name == 'sin':
        s = np.sin(h)
        return s, np.cos(h), -s
    if name == 'gauss':
        e = np.exp(-h * h)
        return e, -2 * h * e, (4 * h * h - 2) * e
    raise ValueError(name)


def _cheb_s(M):
    """Chebyshev-Gauss-Lobatto points on [0, 1] (clusters rows at the interface and the boundary)"""
    return 0.5 * (1 - np.cos(np.pi * np.arange(M) / (M - 1)))


class Features:
    """Random hidden layer of one subdomain network; returns values and x/z derivatives."""

    def __init__(self, n, K, activation='tanh', basis='rfm', z_range=(-1.0, 1.0), rx=1.0, rz=3.0, rb=None,
                 rng=None):
        rng = rng or np.random.default_rng(0)
        self.K, self.act, self.basis = K, activation, basis
        self.z0, self.z1 = z_range
        self.dzs = 2.0 / (self.z1 - self.z0)
        rb = rz if rb is None else rb
        if basis == 'rfm':
            self.A = rng.uniform(-rx, rx, (K, n))
            self.B = rng.uniform(-rx, rx, (K, n))
            self.k = np.arange(1, K + 1)
        else:
            nz = max(1, n // (2 * K + 1))
            self.k = np.arange(-K, K + 1)
            n = nz * (2 * K + 1)
        self.a = rng.uniform(-rz, rz, n if basis == 'rfm' else nz)
        self.b = rng.uniform(-rb, rb, n if basis == 'rfm' else nz)
        self.n = n

    def __call__(self, x, z, order=2):
        """dict f, x, z, xx, zz of shape (len(x), n) (complex for 'separable')"""
        zs = (2 * z - (self.z0 + self.z1)) / (self.z1 - self.z0)
        if self.basis == 'rfm':
            kx = x[:, None] * self.k[None, :]
            c, s = np.cos(kx), np.sin(kx)
            h = c @ self.A + s @ self.B + zs[:, None] * self.a + self.b
            hx = (-self.k * s) @ self.A + (self.k * c) @ self.B
            hxx = (-self.k ** 2 * c) @ self.A + (-self.k ** 2 * s) @ self.B
            hz = self.a * self.dzs
            a0, a1, a2 = _act(self.act, h)
            out = dict(f=a0, x=a1 * hx, z=a1 * hz)
            if order >= 2:
                out['xx'] = a2 * hx * hx + a1 * hxx
                out['zz'] = a2 * hz * hz + 0 * h
                out['xz'] = a2 * hx * hz
            return out
        # separable: e^{ikx} * psi_j(z), flattened over (k, j)
        E = np.exp(1j * x[:, None] * self.k[None, :])                           # (m, 2K+1)
        a0, a1, a2 = _act(self.act, zs[:, None] * self.a + self.b)              # (m, nz)
        hz = self.a * self.dzs
        prod = lambda P, Q: (P[:, :, None] * Q[:, None, :]).reshape(len(x), -1)
        out = dict(f=prod(E, a0), x=prod(1j * self.k * E, a0), z=prod(E, a1 * hz))
        if order >= 2:
            out['xx'] = prod(-(self.k ** 2) * E, a0)
            out['zz'] = prod(E, a2 * hz * hz)
            out['xz'] = prod(1j * self.k * E, a1 * hz)
        return out


G_XX = {'cosx': lambda x: -np.cos(x), 'cos4x': lambda x: -16 * np.cos(4 * x), 'cos4x_over4': lambda x: -4 * np.cos(4 * x),
        'sinx': lambda x: -np.sin(x), 'cos2x': lambda x: -4 * np.cos(2 * x)}


class MappedFeatures:
    """Features that live in the flattened (TFE) coordinate of their own layer,
        upper:  zeta = (z - g(x)) / (a - g(x)),     lower:  zeta = (z + b) / (g(x) + b),     zeta in [0, 1],
    i.e. phi(x, z) = Phi(x, zeta(x, z)), with physical derivatives by the chain rule.  This is the change
    of variables of the paper's TFE method: the network only has to represent u, w ON their own layer,
    whereas features in physical (x, z) must also continue u below the crests of the interface (and w above
    its troughs) -- the Rayleigh-hypothesis difficulty, which grows with the slope of g."""

    def __init__(self, F, layer, G):
        self.F, self.layer, self.G = F, layer, G
        self.n = F.n

    def __call__(self, x, z, order=2):
        G = self.G
        g, gx = G.g(x), G.g_x(x)
        gxx = G.eps * G_XX[G.profile](x)
        if self.layer == 'u':
            L = G.a - g
            zeta = (z - g) / L
            zz_ = 1 / L
            zx = gx * (z - G.a) / L ** 2
            zxx = gxx * (z - G.a) / L ** 2 + 2 * gx ** 2 * (z - G.a) / L ** 3
            zxz = gx / L ** 2
        else:
            L = g + G.b
            zeta = (z + G.b) / L
            zz_ = 1 / L
            zx = -(z + G.b) * gx / L ** 2
            zxx = -(z + G.b) * (gxx / L ** 2 - 2 * gx ** 2 / L ** 3)
            zxz = -gx / L ** 2
        D = self.F(x, zeta, max(order, 1) if order < 2 else 2)
        c = lambda v: v[:, None]
        out = dict(f=D['f'])
        if order >= 1:
            out['x'] = D['x'] + D['z'] * c(zx)
            out['z'] = D['z'] * c(zz_)
        if order >= 2:
            out['xx'] = D['xx'] + 2 * D['xz'] * c(zx) + D['zz'] * c(zx ** 2) + D['z'] * c(zxx)
            out['zz'] = D['zz'] * c(zz_ ** 2)
        return out


class LSQPinn:
    """Linear-least-squares (random-feature) interface PINN for one grating (eps, omega fixed)."""

    def __init__(self, grating, K=6, n_feat=1200, activation='tanh', basis='rfm', Mx=None, Ms=28,
                 Nx_bc=64, lam_bc=10.0, lam_if=10.0, rx=1.0, rz=3.0, rb=None, seed=0, rcond=1e-13, features=None,
                 coords='physical'):
        self.G = grating
        self.K = K
        act_u, act_w = (activation, activation) if isinstance(activation, str) else activation
        nu, nw = (n_feat, n_feat) if np.isscalar(n_feat) else n_feat
        rng = np.random.default_rng(seed)
        zr = (-grating.b, grating.a)
        if features is None:
            if coords == 'tfe':      # features in the flattened coordinate zeta in [0, 1] of each layer
                self.Fu = MappedFeatures(Features(nu, K, act_u, basis, (0.0, 1.0), rx, rz, rb, rng), 'u', grating)
                self.Fw = MappedFeatures(Features(nw, K, act_w, basis, (0.0, 1.0), rx, rz, rb, rng), 'w', grating)
            else:
                self.Fu = Features(nu, K, act_u, basis, zr, rx, rz, rb, rng)
                self.Fw = Features(nw, K, act_w, basis, zr, rx, rz, rb, rng)
        else:                                   # any callables (x, z, order) -> dict (e.g. trained hidden layers)
            self.Fu, self.Fw = features
        self.Mx = Mx or max(64, 8 * K)
        self.Ms, self.Nx_bc = Ms, Nx_bc
        self.lam_bc, self.lam_if, self.rcond = lam_bc, lam_if, rcond

    # ---- exact flat solution (moved to the right-hand side) --------------
    def flat(self, layer, z):
        G = self.G
        gu, gw, t2 = complex(G.gamma_u), complex(G.gamma_w), complex(G.tau2)
        if layer == 'u':
            amp, lam = (gu - t2 * gw) / (gu + t2 * gw), 1j * gu
        else:
            amp, lam = 2 * gu / (gu + t2 * gw), -1j * gw
        f = amp * np.exp(lam * z)
        return f, lam * f, lam * lam * f

    def points(self):
        G = self.G
        x = TWO_PI * (np.arange(self.Mx) + 0.5) / self.Mx
        s = _cheb_s(self.Ms)[1:-1]                 # interior rows; the edges carry the (6c-f) rows
        X, S = np.meshgrid(x, s, indexing='ij')
        g = G.g(X)
        pts = dict(int_u=(X.ravel(), (g + S * (G.a - g)).ravel()), int_w=(X.ravel(), (-G.b + S * (g + G.b)).ravel()))
        xi = TWO_PI * np.arange(2 * self.Mx) / (2 * self.Mx)
        pts['if'] = (xi, G.g(xi))
        xb = TWO_PI * np.arange(self.Nx_bc) / self.Nx_bc
        pts['top'] = (xb, np.full_like(xb, G.a))
        pts['bot'] = (xb, np.full_like(xb, -G.b))
        return pts

    # ---- assemble A c = rhs --------------------------------------------------
    def assemble(self):
        G = self.G
        k2 = {'u': complex(G.k_u) ** 2, 'w': complex(G.k_w) ** 2}
        al, gu, t2 = G.alpha, complex(G.gamma_u), complex(G.tau2)
        pts = self.points()
        F = {'u': self.Fu, 'w': self.Fw}
        nu, nw = self.Fu.n, self.Fw.n
        blocks, rhs = [], []

        def add(Au, Aw, r, wgt):
            m = len(r)
            Au = np.zeros((m, nu)) if Au is None else Au
            Aw = np.zeros((m, nw)) if Aw is None else Aw
            blocks.append(wgt * np.hstack([Au, Aw]))
            rhs.append(wgt * r)

        for layer in ('u', 'w'):
            x, z = pts['int_' + layer]
            d = F[layer](x, z, 2)
            A = d['xx'] + d['zz'] + 2j * al * d['x'] + (k2[layer] - al ** 2) * d['f']
            w = 1 / (1 + abs(k2[layer])) / np.sqrt(len(x))
            add(A if layer == 'u' else None, A if layer == 'w' else None, np.zeros(len(x), complex), w)
        # interface
        x, z = pts['if']
        du, dw = F['u'](x, z, 1), F['w'](x, z, 1)
        u0, u0z, _ = self.flat('u', z)
        w0, w0z, _ = self.flat('w', z)
        gx = G.g_x(x)
        inc = np.exp(-1j * gu * z)
        zeta, psi = -inc, (1j * gu + 1j * al * gx) * inc
        m = len(x)
        add(du['f'], -dw['f'], zeta - (u0 - w0), np.sqrt(self.lam_bc / m))
        Nu = du['z'] - gx[:, None] * du['x'] - 1j * al * gx[:, None] * du['f']
        Nw = dw['z'] - gx[:, None] * dw['x'] - 1j * al * gx[:, None] * dw['f']
        flux0 = (u0z - 1j * al * gx * u0) - t2 * (w0z - 1j * al * gx * w0)
        add(Nu, -t2 * Nw, psi - flux0, np.sqrt(self.lam_if / m) / (1 + abs(gu)))
        # transparent boundary conditions (flat part satisfies them exactly -> rhs 0)
        for key, layer, sgn in (('top', 'u', 1), ('bot', 'w', -1)):
            x, z = pts[key]
            d = F[layer](x, z, 1)
            gp = G.gamma_p(len(x), layer)
            T = np.fft.ifft(sgn * 1j * gp[:, None] * np.fft.fft(d['f'], axis=0), axis=0)
            A = d['z'] - T
            w = np.sqrt(self.lam_bc / len(x)) / (1 + abs(k2[layer]) ** 0.5)
            add(A if layer == 'u' else None, A if layer == 'w' else None, np.zeros(len(x), complex), w)
        return np.vstack(blocks).astype(complex), np.concatenate(rhs).astype(complex)

    def solve(self):
        t0 = time.time()
        A, r = self.assemble()
        cn = np.linalg.norm(A, axis=0)
        cn[cn == 0] = 1
        c, _, rank, sv = sla.lstsq(A / cn, r, cond=self.rcond, lapack_driver='gelsd')
        c = c / cn
        self.c_u, self.c_w = c[:self.Fu.n], c[self.Fu.n:]
        self.loss = float(np.sum(np.abs(A @ c - r) ** 2))   # = PINN loss (same weights as pinn.py)
        self.rank, self.shape, self.cond = rank, A.shape, sv[0] / sv[rank - 1]
        self.solve_time = time.time() - t0
        return self

    # ---- evaluation ------------------------------------------------------
    def evaluate(self, layer, x, z):
        x, z = np.atleast_1d(np.asarray(x, float)), np.atleast_1d(np.asarray(z, float))
        F, c = (self.Fu, self.c_u) if layer == 'u' else (self.Fw, self.c_w)
        return self.flat(layer, z)[0] + F(x, z, 0)['f'] @ c

    def energy(self, Nx=64):
        x = TWO_PI * np.arange(Nx) / Nx
        return self.G.energy(self.evaluate('u', x, np.full(Nx, self.G.a)), self.evaluate('w', x, np.full(Nx, -self.G.b)))


def lsq_energy(G, **kw):
    """R, T, D of the LSQ-PINN at the grating G (convenience)."""
    return LSQPinn(G, **kw).solve().energy()


def lsq_map(G0, Eps, omega, **kw):
    """LSQ-PINN reflectivity map: one least-squares solve per (eps, omega) grid point."""
    R, T, D = (np.zeros((len(Eps), len(omega))) for _ in range(3))
    for i, e in enumerate(Eps):
        for j, o in enumerate(omega):
            G = replace(G0, eps=float(e), omega=float(o), alpha=G0.alpha * o / G0.omega)
            R[i, j], T[i, j], D[i, j] = lsq_energy(G, **kw)
    return R, T, D
