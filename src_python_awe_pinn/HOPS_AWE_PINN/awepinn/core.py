"""HOPS/AWE + PINN (physics-informed summation) for EVERY refl_map.py scenario.

The best combined method of the previous studies (HOPS_PINN_Hybrid): the HOPS/AWE coefficient fields
u_{n,m}, w_{n,m} of each frequency window are the hidden layer of a least-squares interface PINN, and at
each (eps, omega) the output weights come from one least-squares solve of the PINN loss of eqs. (6)
(Helmholtz in physical coordinates through the TFE map, interface conditions, exact DtN TBCs).
Here it is generalised to all windows of a refl_map.py run (several bands, named dispersive materials
with max_delta windows, 'joint' windows, TE/TM, any profile) and made adaptive:

  1. refl_map.run(..., keep_fields=True)        -> the AWE results, exactly as refl_map.py computes them
  2. per window: indicator = PINN loss of the AWE (full-order Taylor) sum at every grid point (~5 ms)
  3. where indicator > tol: least-squares physics-informed weights over the (POD-compressed) basis
     with orders n, m <= basis_order (the higher orders add nothing once the weights are free)
  4. the window result gets ru, rl, ee (and RR) replaced by the hybrid values; the AWE values are kept as
     ru_awe, rl_awe, ee_awe, together with the indicator map and the per-window statistics

The result has the refl_map.py format, so refl_map.plot / refl_movie work on it unchanged.
"""
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_ROOT = os.path.dirname(HERE)
for _p in (os.environ.get('HOPS_PINN_HYBRID', os.path.join(_ROOT, 'HOPS_PINN_Hybrid')),
           os.environ.get('HOPS_PYTHON', os.path.join(_ROOT, 'HOPS_Python'))):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from hybrid.core import PISum, spectral_dx, TWO_PI, Grating2D    # noqa: E402
import refl_map as rm                                              # noqa: E402
from hops import cheb                                              # noqa: E402

DEFAULT_TOL = 1e-14          # indicator threshold (PINN loss); calibrated in error_analysis.py
DEFAULT_ORDER = 12           # basis orders n, m <= 12 used by the least-squares weights
COMPRESS = 1e-13             # POD tolerance of the basis


class WindowBand:
    """band-like object (what PISum needs) built from ONE refl_map window result r and the run info"""

    def __init__(self, r, info, basis_order=DEFAULT_ORDER):
        self.r = r
        self.Nx, self.Nz, self.N, self.M = info['Nx'], info['Nz'], info['N'], info['M']
        self.a, self.b = info['a'], info['b']
        self.n_u, self.n_w = complex(r['n_u']), complex(r['n_w'])
        self.mode, self.alpha_bar, self.omega_bar = info['mode'], info['alpha'], r['omega_bar']
        self.x = TWO_PI * np.arange(self.Nx) / self.Nx
        self.f = np.asarray(info['f'], float)
        self.fx = np.real(spectral_dx(self.f.astype(complex)))
        self.fxx = np.real(spectral_dx(self.f.astype(complex), 2))
        D, t = cheb(self.Nz)
        self.zeta = 0.5 * (t + 1)
        self.Dzeta = 2 * D
        self.Eps, self.delta, self.omega = r['Eps'], r['delta'], r['omega']
        # full-order Taylor weights use ALL coefficients; the least-squares basis only n, m <= basis_order
        self.nm = [(n, m) for m in range(self.M + 1) for n in range(self.N + 1)]
        nb = min(basis_order, self.N), min(basis_order, self.M)
        self.sub = np.array([k for k, (n, m) in enumerate(self.nm) if n <= nb[0] and m <= nb[1]])
        K = len(self.nm)
        self.basis, self.raw = {}, {}
        for layer, key in (('u', 'u_n_m'), ('w', 'w_n_m')):
            B = r[key].reshape(self.Nx, self.Nz + 1, K)
            bad = ~np.all(np.isfinite(B), axis=(0, 1))
            if bad.any():
                B = np.where(bad[None, None, :], 0, B)
            self.raw[layer] = B                                # all orders: for the (full) Taylor sum
            self.basis[layer] = self._derivs(B[..., self.sub])  # derivatives only of the least-squares basis

    def _derivs(self, B):
        Bx = spectral_dx(B)
        Dz = self.Dzeta
        return dict(f=B, x=Bx, z=np.einsum('lk,jkq->jlq', Dz, B), xx=spectral_dx(B, 2),
                    zz=np.einsum('lk,jkq->jlq', Dz @ Dz, B), xz=np.einsum('lk,jkq->jlq', Dz, Bx))

    def taylor_weights(self, eps, delta):
        return np.array([eps ** n * delta ** m for (n, m) in self.nm], dtype=complex)

    def grating(self, eps, omega):
        return Grating2D(n_u=self.n_u, n_w=self.n_w, omega=omega, eps=eps,
                         alpha=self.alpha_bar * omega / self.omega_bar, a=self.a, b=self.b, mode=self.mode)


class HybridSum(PISum):
    """PISum whose indicator uses the FULL AWE series (all n <= N, m <= M) while the least-squares
    weights use the truncated, POD-compressed basis."""

    def taylor_fields(self, eps, omega):
        tw = self.B.taylor_weights(eps, omega / self.B.omega_bar - 1)
        # sum first, differentiate the single summed field (linear: identical to summing the derivatives)
        return {layer: self.B._derivs((self.B.raw[layer] @ tw)[..., None]) for layer in ('u', 'w')}

    def fields(self, eps, omega, c=None):
        """nodal fields u, w (Nx, Nz+1) of the hybrid solution at (eps, omega) (for pictures / movies)"""
        if c is None:
            c = self.solve(eps, omega)['c']
        nu = self.nfeat['u']
        return self.feat['u']['f'] @ c[:nu], self.feat['w']['f'] @ c[nu:]

    def point(self, eps, omega, tol=DEFAULT_TOL):
        """adaptive evaluation at one point: Taylor where the indicator is below tol, else least squares"""
        s = self.indicator(eps, omega)
        if s['loss_awe_taylor'] <= tol:
            s.update(R=s['R_taylor'], T=s['T_taylor'], D=s['D_taylor'], loss=s['loss_awe_taylor'], solved=False)
            return s
        A, r, G = self.assemble(eps, omega)
        import scipy.linalg as sla
        cn = np.linalg.norm(A, axis=0)
        cn[cn == 0] = 1
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            c = sla.lstsq(A / cn, r, cond=self.rcond, lapack_driver='gelsd')[0] / cn
        nu = self.nfeat['u']
        s.update(self._energy(G, c[:nu], c[nu:]), loss=float(np.sum(np.abs(A @ c - r) ** 2)), solved=True, c=c)
        return s


def hybrid_window(args):
    """process one window (module level: runs in worker processes)"""
    r, info, tol, basis_order, keep = args
    t0 = time.time()
    B = WindowBand(r, info, basis_order)
    S = HybridSum(B, compress=COMPRESS)
    ne, nd = len(r['Eps']), len(r['delta'])
    R, T, D, ind, solved = (np.zeros((ne, nd)) for _ in range(5))
    for i, e in enumerate(r['Eps']):
        for j, dl in enumerate(r['delta']):
            s = S.point(float(e), float(r['omega_bar'] * (1 + dl)), tol)
            R[i, j], T[i, j], D[i, j], ind[i, j], solved[i, j] = s['R'], s['T'], s['D'], s['loss_awe_taylor'], s['solved']
    out = dict(R=R, T=T, D=D, indicator=ind, solved=solved.astype(bool), time=time.time() - t0,
               n_unknowns=int(S.nfeat['u'] + S.nfeat['w']))
    if keep:
        out['S'] = None
    return out


def run(scenario, qq=None, N_Eps=31, N_delta=31, tol=DEFAULT_TOL, basis_order=DEFAULT_ORDER, workers=1,
        keep_fields=False, verbose=True, **over):
    """refl_map.run for the hybrid; returns (results, info) in refl_map format (+ AWE copies, indicator)."""
    t0 = time.time()
    res, info = rm.run(scenario, qq=qq, N_Eps=N_Eps, N_delta=N_delta, verbose=False, workers=workers,
                       keep_fields=True, **over)
    t_awe = time.time() - t0
    args = [(r, info, tol, basis_order, False) for r in res]
    if workers > 1 and len(res) > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            hyb = list(ex.map(hybrid_window, args))
    else:
        hyb = [hybrid_window(a) for a in args]
    for r, h in zip(res, hyb):
        r['ru_awe'], r['rl_awe'], r['ee_awe'], r['RR_awe'] = r['ru'], r['rl'], r['ee'], r['RR']
        r['ru'], r['rl'], r['ee'] = h['R'] + 0j, h['T'] + 0j, h['D'] + 0j
        r['RR'] = r['ru'] / r['ru_flat'] if info['relative'] else r['ru']
        r['indicator'], r['solved'], r['t_hybrid'], r['n_unknowns'] = h['indicator'], h['solved'], h['time'], h['n_unknowns']
        if not keep_fields:
            for k in ('u_n_m', 'w_n_m'):
                r.pop(k, None)
    info['t_awe'] = t_awe
    info['t_hybrid'] = sum(h['time'] for h in hyb) / max(1, min(workers, len(res)))
    info['tol'], info['basis_order'] = tol, basis_order
    if verbose:
        ns = sum(int(r['solved'].sum()) for r in res)
        nt = sum(r['solved'].size for r in res)
        print(f'{scenario}: {len(res)} windows, AWE {t_awe:.1f} s, hybrid {info["t_hybrid"]:.0f} s '
              f'({ns}/{nt} points needed the least-squares weights)', flush=True)
    return res, info
