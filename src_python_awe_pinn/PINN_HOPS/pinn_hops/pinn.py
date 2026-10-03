"""Physics-Informed Neural Network for the two-layer grating problem (6a)-(6h).

Design (and why)
----------------
* Two networks, u_theta (upper layer) and w_theta (lower layer).  The normal derivative jumps across
  z = g(x) (tau^2 != 1 in TM, and the incident field is not part of u, w), so one smooth network for both
  layers would fight the interface conditions.  Each outputs (Re, Im) of its complex field.
* Periodicity (6g, 6h) is imposed EXACTLY: x enters only through Fourier features cos(kx), sin(kx),
  k = 1..K, so every network output is 2 pi-periodic by construction (no periodicity loss).
* The transparent boundary conditions (6e, 6f) are the exact nonlocal DtN operators of the paper, not
  approximations: on z = a (z = -b) the network is evaluated on a uniform x-grid, FFT'ed, multiplied by
  i gamma^u_p (-i gamma^w_p) and inverse-FFT'ed -- all inside torch, so the loss is differentiable.
  (A local absorbing condition would add a modelling error HOPS does not have.)
* Collocation in the physical domain: z = g(x) + s (a - g(x)), s ~ U(0, 1) (upper) and
  z = -b + s (g(x) + b) (lower), i.e. uniform in the TFE variable z' of the paper; resampled regularly.
* Residuals are nondimensionalised: PDE residual / (1 + |k|^2), derivative conditions / (1 + |k|);
  (6c) and the TBCs get weight lambda_bc, (6d) weight lambda_if (defaults 10: they carry the data).
* Optional parametric inputs (eps, omega): u_theta(x, z; eps, omega).  One network then covers a whole
  (eps, omega) window -- the PINN counterpart of the HOPS/AWE joint expansion in (eps, delta).  The
  domain moves with eps (z = eps f(x)) and every PDE/TBC coefficient is evaluated per sample.
* float64 throughout; Adam (with decay) followed by L-BFGS (strong Wolfe), the standard recipe for
  reaching PINN accuracy limits.
"""
import math
import time
import numpy as np
import torch
import torch.nn as nn

from .problem import outgoing_sqrt

torch.set_default_dtype(torch.float64)
TWO_PI = 2 * math.pi


# ----------------------------------------------------------------------------
# networks
# ----------------------------------------------------------------------------
class FieldNet(nn.Module):
    """(x, z[, eps, omega]) -> (Re, Im) of a d-periodic complex field."""

    def __init__(self, K=8, width=64, depth=4, z_range=(-1.0, 1.0), param_ranges=None, activation='tanh'):
        super().__init__()
        self.K = K
        # activation: 'tanh' | 'sin' | 'laaf' (layer-wise adaptive tanh(n a h), Jagtap et al. 2020, n = 10)
        self.activation = activation
        if activation == 'laaf':
            self.laaf_n = 10.0
            self.laaf_a = nn.Parameter(torch.full((depth,), 0.1))
        self.register_buffer('kk', torch.arange(1, K + 1, dtype=torch.float64))
        self.z0, self.z1 = z_range
        self.param_ranges = list(param_ranges or [])
        n_in = 2 * K + 1 + len(self.param_ranges)
        layers, d = [], n_in
        for _ in range(depth):
            layers += [nn.Linear(d, width), _Act(self, len(layers) // 2)]
            d = width
        layers += [nn.Linear(d, 2)]
        self.mlp = nn.Sequential(*layers)
        for m in self.mlp:
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)

    def derivs(self, x, z, params=(), order=2):
        """Values and x/z derivatives (up to second order) by forward-mode propagation through the
        tanh MLP (Taylor-mode AD): ~5x cheaper than nested reverse-mode autograd for the Laplacian.
        Returns dict with keys 'f', 'x', 'z' (and 'xx', 'zz'), each (n, 2) = (Re, Im) parts."""
        n = x.shape[0]
        kk = self.kk
        kx = x[:, None] * kk[None, :]
        c, s_ = torch.cos(kx), torch.sin(kx)
        zs = (2 * z - (self.z0 + self.z1)) / (self.z1 - self.z0)
        dz_s = 2.0 / (self.z1 - self.z0)
        pf = [((2 * p - (lo + hi)) / (hi - lo))[:, None] for p, (lo, hi) in zip(params, self.param_ranges)]
        zero = torch.zeros(n, 1 + len(pf))
        h = torch.cat([c, s_, zs[:, None]] + pf, dim=1)
        hx = torch.cat([-kk * s_, kk * c, zero], dim=1)
        hz = torch.cat([torch.zeros(n, 2 * self.K), torch.full((n, 1), dz_s), torch.zeros(n, len(pf))], dim=1)
        hxx = torch.cat([-kk ** 2 * c, -kk ** 2 * s_, zero], dim=1) if order >= 2 else None
        hzz = torch.zeros_like(h) if order >= 2 else None
        layers = [m for m in self.mlp if isinstance(m, nn.Linear)]
        for i, lin in enumerate(layers):
            W = lin.weight.T
            h, hx, hz = h @ W + lin.bias, hx @ W, hz @ W
            if order >= 2:
                hxx, hzz = hxx @ W, hzz @ W
            if i < len(layers) - 1:
                a, d1, d2 = self.act(h, i, order)
                if order >= 2:
                    hxx = d2 * hx * hx + d1 * hxx
                    hzz = d2 * hz * hz + d1 * hzz
                h, hx, hz = a, d1 * hx, d1 * hz
        out = dict(f=h, x=hx, z=hz)
        if order >= 2:
            out.update(xx=hxx, zz=hzz)
        return out

    def act(self, h, i, order=2):
        """activation value and its first two derivatives"""
        if self.activation == 'sin':
            a = torch.sin(h)
            return a, torch.cos(h), (-a if order >= 2 else None)
        s = self.laaf_n * self.laaf_a[i] if self.activation == 'laaf' else 1.0
        t = torch.tanh(s * h)
        d = 1 - t * t
        return t, s * d, (-2 * s * s * t * d if order >= 2 else None)

    def forward(self, x, z, params=()):
        kx = x[:, None] * self.kk[None, :]
        zs = (2 * z - (self.z0 + self.z1)) / (self.z1 - self.z0)
        feats = [torch.cos(kx), torch.sin(kx), zs[:, None]]
        for p, (lo, hi) in zip(params, self.param_ranges):
            feats.append(((2 * p - (lo + hi)) / (hi - lo))[:, None])
        out = self.mlp(torch.cat(feats, dim=1))
        return out[:, 0], out[:, 1]


# ----------------------------------------------------------------------------
# problem coefficients per sample (supports parametric eps, omega)
# ----------------------------------------------------------------------------
class Coeffs:
    """Physical coefficients as torch tensors, broadcast over samples.

    The grating (n_u, n_w, profile, alpha_ref at omega_ref, a, b, mode) is fixed; eps and omega may vary
    per sample.  alpha(omega) = alpha_ref * omega / omega_ref (fixed incidence angle; alpha_ref = 0 is
    normal incidence as in the paper's examples)."""

    def __init__(self, grating, omega_ref=None):
        self.G = grating
        self.omega_ref = omega_ref or grating.omega
        self.tau2 = complex(grating.tau2)

    def k2(self, omega, layer):
        n = self.G.n_u if layer == 'u' else self.G.n_w
        return (n * omega / self.G.c0) ** 2

    def alpha(self, omega):
        return self.G.alpha * omega / self.omega_ref

    def gamma_u(self, omega):
        return torch.sqrt((self.k2(omega, 'u') - self.alpha(omega) ** 2).to(torch.complex128))

    def f(self, x):
        name = self.G.profile
        return {'cosx': torch.cos(x), 'cos4x': torch.cos(4 * x), 'cos4x_over4': 0.25 * torch.cos(4 * x),
                'sinx': torch.sin(x), 'cos2x': torch.cos(2 * x)}[name]

    def f_x(self, x):
        name = self.G.profile
        return {'cosx': -torch.sin(x), 'cos4x': -4 * torch.sin(4 * x), 'cos4x_over4': -torch.sin(4 * x),
                'sinx': torch.cos(x), 'cos2x': -2 * torch.sin(2 * x)}[name]


class _Act(nn.Module):
    def __init__(self, net, i):
        super().__init__()
        self.i = i
        object.__setattr__(self, 'net', net)          # no submodule cycle

    def forward(self, h):
        return self.net.act(h, self.i, order=0)[0]


def _grad(f, v):
    g = torch.autograd.grad(f.sum(), v, create_graph=True, allow_unused=True)[0]
    return torch.zeros_like(v) if g is None else g


def _grads(fr, fi, x, z, second=True):
    """First (and second) derivatives of the real and imaginary parts."""
    out = {}
    for name, f in (('r', fr), ('i', fi)):
        fx, fz = _grad(f, x), _grad(f, z)
        out['x' + name], out['z' + name] = fx, fz
        if second:
            out['xx' + name] = _grad(fx, x) if fx.requires_grad else torch.zeros_like(x)
            out['zz' + name] = _grad(fz, z) if fz.requires_grad else torch.zeros_like(z)
    c = lambda k: torch.complex(out[k + 'r'], out[k + 'i'])
    return {k: c(k) for k in (('x', 'z', 'xx', 'zz') if second else ('x', 'z'))}


# ----------------------------------------------------------------------------
# the PINN
# ----------------------------------------------------------------------------
class GratingPINN:
    def __init__(self, grating, K=8, width=64, depth=4, eps_range=None, omega_range=None,
                 lam_bc=10.0, lam_if=10.0, Nx_bc=64, n_int=2000, n_if=256, seed=0, ansatz='flat',
                 corr_scale=1.0, activation='tanh'):
        torch.manual_seed(seed)
        np.random.seed(seed)
        self.G = grating
        self.C = Coeffs(grating, omega_ref=grating.omega)
        self.eps_range = eps_range            # None -> fixed eps = grating.eps
        self.omega_range = omega_range        # None -> fixed omega = grating.omega
        pr = [r for r in (eps_range, omega_range) if r is not None]
        zr = (-grating.b, grating.a)
        # activation='tanh' | 'sin' | 'laaf', or a pair (upper, lower): the I-PINN choice of Sarma et al.
        # (CMAME 2024) -- a different activation function in each subdomain
        act_u, act_w = (activation, activation) if isinstance(activation, str) else activation
        self.net_u = FieldNet(K, width, depth, zr, pr, activation=act_u)
        self.net_w = FieldNet(K, width, depth, zr, pr, activation=act_w)
        self.lam_bc, self.lam_if = lam_bc, lam_if
        self.Nx_bc, self.n_int, self.n_if = Nx_bc, n_int, n_if
        self.n_param = 8 if pr else 1          # parameter samples per batch (parametric mode)
        self.history = []
        self.ansatz, self.corr_scale = ansatz, corr_scale
        self.s_layer = {'u': 1.0, 'w': 1.0}
        self.override = None

    # ---- parameters -------------------------------------------------------
    def parameters(self):
        return list(self.net_u.parameters()) + list(self.net_w.parameters())

    def _params(self, eps, omega):
        p = []
        if self.eps_range is not None:
            p.append(eps)
        if self.omega_range is not None:
            p.append(omega)
        return p

    def _draw_params(self, n):
        e = (np.random.uniform(*self.eps_range, n) if self.eps_range is not None else np.full(n, self.G.eps))
        o = (np.random.uniform(*self.omega_range, n) if self.omega_range is not None else np.full(n, self.G.omega))
        return e, o

    # ---- field evaluation -------------------------------------------------
    def flat(self, layer, z, omega):
        """Exact flat-interface solution and its z-derivatives: u0 = r e^{i gamma^u z}, w0 = t e^{-i gamma^w z}."""
        C = self.C
        al = C.alpha(omega)
        gu = _outgoing_sqrt_t(C.k2(omega, 'u') - al ** 2)
        gw = _outgoing_sqrt_t(C.k2(omega, 'w') - al ** 2)
        t2 = C.tau2
        if layer == 'u':
            amp, lam = (gu - t2 * gw) / (gu + t2 * gw), 1j * gu
        else:
            amp, lam = 2 * gu / (gu + t2 * gw), -1j * gw
        f = amp * torch.exp(lam * z)
        return f, lam * f, lam * lam * f

    def fields(self, layer, x, z, e, o, order=2):
        """Complex field and derivatives {'f','x','z'[,'xx','zz']} at the points.

        Ansatz (default, ansatz='flat'):  u = u_flat + eps * N_u,  w = w_flat + eps * N_w,
        where u_flat, w_flat is the EXACT solution for the flat interface (Fresnel), which satisfies
        (6a,b,e,f) identically.  The networks learn only the O(eps) correction (the analogue of HOPS'
        expansion about the flat interface), eps = 0 is reproduced exactly, and the loss is not dominated
        by the O(1) specular part.  ansatz='plain' uses the bare networks."""
        if getattr(self, 'override', None) is not None:     # any other torch field (tests: exact, HOPS)
            fr, fi = self.override(layer, x, z, e, o)
            d = _grads(fr, fi, x, z, second=order >= 2) if order >= 1 else {}
            d['f'] = torch.complex(fr, fi)
            return d
        net = self.net_u if layer == 'u' else self.net_w
        D = net.derivs(x, z, self._params(e, o), order=max(order, 1))
        cx = lambda k: torch.complex(D[k][:, 0], D[k][:, 1])
        keys = ['f', 'x', 'z'] + (['xx', 'zz'] if order >= 2 else [])
        out = {k: cx(k) for k in keys}
        if self.ansatz == 'flat':
            sc = e * self.corr_scale * self.s_layer[layer]
            out = {k: sc * v for k, v in out.items()}
            f0, f0z, f0zz = self.flat(layer, z, o)
            out['f'] = out['f'] + f0
            out['z'] = out['z'] + f0z
            if order >= 2:
                out['zz'] = out['zz'] + f0zz
        return out

    def field(self, layer, x, z, eps, omega):
        if getattr(self, 'override', None) is not None:     # plug in any other torch field (tests, HOPS)
            return self.override(layer, x, z, eps, omega)
        net = self.net_u if layer == 'u' else self.net_w
        fr, fi = net(x, z, self._params(eps, omega))
        return fr, fi

    # ---- sampling ---------------------------------------------------------
    def sample(self):
        """Collocation points: interior (both layers), interface, and uniform-x boundary grids."""
        G, C = self.G, self.C
        P = self.n_param
        e, o = self._draw_params(P)
        S = {}
        # interior points, n_int per layer spread over the P parameter samples
        m = max(1, self.n_int // P)
        x = np.random.uniform(0, TWO_PI, (P, m))
        s = np.random.uniform(0, 1, (P, m))
        g = e[:, None] * G.f(x)
        S['int_u'] = (x, g + s * (G.a - g), e[:, None] + 0 * x, o[:, None] + 0 * x)
        s = np.random.uniform(0, 1, (P, m))
        S['int_w'] = (x, -G.b + s * (g + G.b), e[:, None] + 0 * x, o[:, None] + 0 * x)
        # interface
        mi = max(16, self.n_if // P)
        x = np.random.uniform(0, TWO_PI, (P, mi))
        S['if'] = (x, e[:, None] * G.f(x), e[:, None] + 0 * x, o[:, None] + 0 * x)
        # boundary grids (uniform in x with a random shift -> FFT-based DtN operators)
        x0 = np.random.uniform(0, TWO_PI / self.Nx_bc, (P, 1))
        x = x0 + TWO_PI * np.arange(self.Nx_bc)[None, :] / self.Nx_bc
        S['top'] = (x, np.full_like(x, G.a), e[:, None] + 0 * x, o[:, None] + 0 * x)
        S['bot'] = (x, np.full_like(x, -G.b), e[:, None] + 0 * x, o[:, None] + 0 * x)
        S['shift'] = x0
        self.S = {k: tuple(torch.tensor(a_, requires_grad=(k != 'shift' and i < 2)) for i, a_ in enumerate(v))
                  if k != 'shift' else torch.tensor(v) for k, v in S.items()}
        return self.S

    # ---- residuals --------------------------------------------------------
    def residuals(self, S=None):
        S = S or self.S
        C, G = self.C, self.G
        out = {}
        for layer in ('u', 'w'):
            x, z, e, o = (t.reshape(-1) for t in S['int_' + layer])
            d = self.fields(layer, x, z, e, o, order=2)
            k2 = C.k2(o, layer)
            al = C.alpha(o)
            res = d['xx'] + d['zz'] + 2j * al * d['x'] + (k2 - al ** 2) * d['f']
            out['pde_' + layer] = res / (1 + abs(k2))
        # interface (6c), (6d)
        x, z, e, o = (t.reshape(-1) for t in S['if'])
        du = self.fields('u', x, z, e, o, order=1)
        dw = self.fields('w', x, z, e, o, order=1)
        u, w = du['f'], dw['f']
        gx = e * C.f_x(x)
        al = C.alpha(o)
        gu = C.gamma_u(o)
        inc = torch.exp(-1j * gu * z)                                # z = g on the interface
        zeta = -inc
        psi = (1j * gu + 1j * al * gx) * inc
        dN = lambda dd: dd['z'] - gx * dd['x']                        # N = (-g_x, 1)
        out['if_dir'] = u - w - zeta
        flux = (dN(du) - 1j * al * gx * u) - C.tau2 * (dN(dw) - 1j * al * gx * w) - psi
        out['if_neu'] = flux / (1 + abs(gu))
        # transparent boundary conditions (6e), (6f): exact DtN via FFT on the uniform x-grids
        for key, layer, sgn in (('top', 'u', +1), ('bot', 'w', -1)):
            x, z, e, o = S[key]
            P, Nx = x.shape
            d = self.fields(layer, x.reshape(-1), z.reshape(-1), e.reshape(-1), o.reshape(-1), order=1)
            f, fz = d['f'].reshape(P, Nx), d['z'].reshape(P, Nx)
            p = torch.tensor(np.fft.fftfreq(Nx, 1.0 / Nx))
            om = o[:, :1]
            k2 = C.k2(om, layer)
            ap = C.alpha(om) + p[None, :]
            gp = _outgoing_sqrt_t(k2 - ap ** 2)
            Tf = torch.fft.ifft(sgn * 1j * gp * torch.fft.fft(f, dim=1), dim=1)
            out['tbc_' + layer] = (fz - Tf).reshape(-1) / (1 + abs(k2).mean() ** 0.5)
        return out

    def loss(self, S=None):
        r = self.residuals(S)
        ms = lambda v: torch.mean(torch.abs(v) ** 2)
        parts = dict(pde=ms(r['pde_u']) + ms(r['pde_w']),
                     if_dir=ms(r['if_dir']), if_neu=ms(r['if_neu']),
                     tbc=ms(r['tbc_u']) + ms(r['tbc_w']))
        L = parts['pde'] + self.lam_bc * (parts['if_dir'] + parts['tbc']) + self.lam_if * parts['if_neu']
        return L, {k: float(v.detach()) for k, v in parts.items()}

    # ---- training ---------------------------------------------------------
    def train(self, adam_iters=3000, lbfgs_iters=1000, lr=2e-3, resample=250, verbose=True, log_every=500,
              time_limit=None):
        t0 = time.time()
        opt = torch.optim.Adam(self.parameters(), lr=lr)
        sched = torch.optim.lr_scheduler.ExponentialLR(opt, gamma=(0.05) ** (1.0 / max(1, adam_iters)))
        self.sample()
        for it in range(adam_iters):
            if it % resample == 0 and it > 0:
                self.sample()
            opt.zero_grad()
            L, parts = self.loss()
            L.backward()
            opt.step()
            sched.step()
            if it % log_every == 0 or it == adam_iters - 1:
                self.history.append(('adam', it, float(L), parts, time.time() - t0))
                if verbose:
                    print(f'  adam {it:5d}  loss {float(L):.3e}  ' +
                          ' '.join(f'{k} {v:.1e}' for k, v in parts.items()) + f'  ({time.time() - t0:.0f} s)',
                          flush=True)
            if time_limit and time.time() - t0 > time_limit:
                break
        if lbfgs_iters:
            # L-BFGS on a larger, fixed point set (deterministic objective)
            self.n_int, self.n_if = 2 * self.n_int, 2 * self.n_if
            n_param0 = self.n_param
            if self.eps_range is not None or self.omega_range is not None:
                self.n_param = 2 * self.n_param
            chunks = max(1, lbfgs_iters // 250)
            for c in range(chunks):
                self.sample()
                opt = torch.optim.LBFGS(self.parameters(), lr=1.0, max_iter=lbfgs_iters // chunks,
                                        history_size=50, tolerance_grad=1e-12, tolerance_change=1e-15,
                                        line_search_fn='strong_wolfe')

                def closure():
                    opt.zero_grad()
                    L, _ = self.loss()
                    L.backward()
                    return L
                opt.step(closure)
                L, parts = self.loss()
                self.history.append(('lbfgs', c, float(L), parts, time.time() - t0))
                if verbose:
                    print(f'  lbfgs chunk {c + 1}/{chunks}  loss {float(L):.3e}  ' +
                          ' '.join(f'{k} {v:.1e}' for k, v in parts.items()) + f'  ({time.time() - t0:.0f} s)',
                          flush=True)
                if time_limit and time.time() - t0 > 2 * time_limit:
                    break
            self.n_int, self.n_if, self.n_param = self.n_int // 2, self.n_if // 2, n_param0
        self.train_time = time.time() - t0
        return self

    # ---- post-processing --------------------------------------------------
    @torch.no_grad()
    def evaluate(self, layer, x, z, eps=None, omega=None):
        eps = self.G.eps if eps is None else eps
        omega = self.G.omega if omega is None else omega
        x = torch.as_tensor(np.asarray(x, dtype=float))
        z = torch.as_tensor(np.asarray(z, dtype=float))
        e = torch.full_like(x, eps)
        o = torch.full_like(x, omega)
        return self.fields(layer, x, z, e, o, order=0)['f'].numpy()

    def energy(self, eps=None, omega=None, Nx=64):
        """R, T, D of the PINN solution (same formula as the paper / energy_defect.m)."""
        G = self.G
        eps = G.eps if eps is None else eps
        omega = G.omega if omega is None else omega
        x = TWO_PI * np.arange(Nx) / Nx
        ut = self.evaluate('u', x, np.full(Nx, G.a), eps, omega)
        wb = self.evaluate('w', x, np.full(Nx, -G.b), eps, omega)
        from dataclasses import replace
        Gp = replace(G, eps=eps, omega=omega, alpha=G.alpha * omega / self.C.omega_ref)
        return Gp.energy(ut, wb)

    def save(self, path):
        torch.save(dict(u=self.net_u.state_dict(), w=self.net_w.state_dict(), history=self.history), path)

    def load(self, path):
        d = torch.load(path)
        self.net_u.load_state_dict(d['u'])
        self.net_w.load_state_dict(d['w'])
        self.history = d.get('history', [])
        return self


def _outgoing_sqrt_t(v):
    """torch version of outgoing_sqrt (Im >= 0), for complex or real tensors."""
    v = v.to(torch.complex128)
    s = torch.sqrt(v)
    return torch.where(s.imag < 0, -s, s)
