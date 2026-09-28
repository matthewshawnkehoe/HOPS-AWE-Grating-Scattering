"""Summation of the double Taylor series  sum_{n,m} c_{n,m} eps^n delta^m.

Ports of: taylorsum.m, taylorsum2.m, taylorsum_2_coeff.m, padesum.m, padesum_safe.m,
padesum2.m, padesum2_safe.m, padeapprox.m, fcn_sum.m, fcn_sum_fast.m, vol_fcn_sum.m.

SumType convention (as in MATLAB):  1 = Taylor, 2 = Pade, 3 = Pade (safe).
Extra (new):                         4 = robust SVD Pade (padeapprox, Gonnet-Guettel-Trefethen).

Polar re-summation (paper Sec. 6.3): eps = rho cos(theta), delta = rho sin(theta),
    ctilde_p = sum_{q=0}^{p} c_{p-q,q} cos^{p-q}(theta) sin^q(theta),  p = 0..min(N,M),
then a one-variable Taylor/Pade sum in rho.

Notes on faithfulness
---------------------
* padesum/padesum_safe solve the Toeplitz system with a dense LU solve exactly as
  MATLAB's backslash does (numpy.linalg.solve == LAPACK gesv).  A *Levinson*
  Toeplitz solver (scipy.linalg.solve_toeplitz) must NOT be used: it needs every
  leading principal minor to be nonsingular and blows up for tiny eps/delta,
  which is the "ill-conditioned Pade" failure seen in earlier ports.
  If the matrix is *exactly* singular MATLAB returns Inf/NaN with a warning; here
  we fall back to the minimum-norm SVD (lstsq) solution instead (set
  ``SINGULAR_FALLBACK = 'nan'`` to reproduce MATLAB's NaNs).
* taylorsum_2_coeff (SumType 1 in fcn_sum/fcn_sum_fast) calls
  taylorsum(coeff, rho, floor(min(N,M)/2)), i.e. the MATLAB code sums ctilde only
  up to rho^floor(min(N,M)/2).  This is reproduced by default; pass
  ``taylor_full_order=True`` to sum all min(N,M)+1 terms.
"""
import warnings
import numpy as np
from scipy.linalg import toeplitz

SINGULAR_FALLBACK = 'lstsq'     # or 'nan' (MATLAB-like)


# ----------------------------------------------------------------------------
# scalar ports (one-to-one with the .m files)
# ----------------------------------------------------------------------------
def polyval_asc(a, x):
    """sum_k a[..., k] x^k  via Horner (== MATLAB polyval(a(end:-1:1), x))."""
    a = np.asarray(a)
    if a.shape[-1] == 0:
        return np.zeros(np.broadcast_shapes(a.shape[:-1], np.shape(x)), dtype=complex)
    s = a[..., -1] * np.ones_like(x, dtype=complex)
    for k in range(a.shape[-1] - 2, -1, -1):
        s = s * x + a[..., k]
    return s


def taylorsum(c, Eps, N):
    return polyval_asc(np.asarray(c)[:N + 1], Eps)


def taylorsum2(c, Eps, delta, N, M):
    c = np.asarray(c)
    vander = np.outer(Eps ** np.arange(N + 1), delta ** np.arange(M + 1))
    return np.sum(c[:N + 1, :M + 1] * vander)


def _polar(Eps, delta):
    return np.sqrt(Eps ** 2 + delta ** 2), np.arctan2(delta, Eps)


def _ctilde(c, cos_t, sin_t, K):
    """c[..., N+1, M+1] -> ctilde[..., K+1]; cos_t/sin_t broadcast over the batch."""
    cos_t = np.asarray(cos_t, dtype=float)
    sin_t = np.asarray(sin_t, dtype=float)
    kk = np.arange(K + 1)
    c1 = cos_t[..., None] ** kk
    c2 = sin_t[..., None] ** kk
    batch = np.broadcast_shapes(c.shape[:-2], cos_t.shape)
    # Fast path: c is (Nx, N+1, M+1) and theta has a trailing singleton axis for x
    # -> each ctilde_p is one matrix product  W_p (..., p+1) @ C_p (p+1, Nx).
    if c.ndim == 3 and cos_t.ndim >= 1 and cos_t.shape[-1] == 1 and cos_t.size > 1:
        coeff = np.empty(batch + (K + 1,), dtype=complex)
        cc = c1[..., 0, :]
        ss = c2[..., 0, :]
        for p in range(K + 1):
            q = np.arange(p + 1)
            Wp = cc[..., p - q] * ss[..., q]                    # (..., p+1)
            Cp = c[:, p - q, q]                                  # (Nx, p+1)
            coeff[..., p] = Wp @ Cp.T
        return coeff
    coeff = np.zeros(batch + (K + 1,), dtype=complex)
    for p in range(K + 1):
        for q in range(p + 1):
            coeff[..., p] = coeff[..., p] + c[..., p - q, q] * c1[..., p - q] * c2[..., q]
    return coeff


def taylorsum_2_coeff(c, Eps, delta, N, M, full_order=False):
    rho, theta = _polar(Eps, delta)
    K = min(N, M)
    coeff = _ctilde(np.asarray(c), np.cos(theta), np.sin(theta), K)
    return taylorsum(coeff, rho, K if full_order else K // 2)


def padesum(c, Eps, M):
    """[L/M] = [M/M] Pade of a single-variable series c_0..c_{2M}; returns (psum, a, b)."""
    c = np.asarray(c).ravel()
    if M == 0:
        a = c[:1].astype(complex)
        b = np.ones(1)
    else:
        H = toeplitz(c[M:2 * M], c[M:0:-1])
        rhs = c[M + 1:2 * M + 1]
        bb = -_solve_one(H, rhs)
        b = np.concatenate([[1.0], bb])
        a = np.convolve(b, c[:M + 1])[:M + 1]
    return polyval_asc(a, Eps) / polyval_asc(b, Eps), a, b


def padesum_safe(c, Eps, M):
    """Pade sum that first discards trailing coefficients with |c_j| <= 1e-14."""
    c = np.asarray(c).ravel()
    idx = np.nonzero(np.abs(c[:2 * M + 1]) > 1e-14)[0]
    N_true = idx[-1] + 1 if idx.size else -1       # 1-based, as in MATLAB
    M_safe = int(np.floor(N_true / 2.0))
    if M_safe < 0:                                  # all coefficients ~ 0 -> psum = 0
        return 0.0 + 0j, np.zeros(0), np.ones(1)
    return padesum(c, Eps, M_safe)


def padesum2(c, Eps, delta, N, M):
    rho, theta = _polar(Eps, delta)
    K = min(N, M)
    coeff = _ctilde(np.asarray(c), np.cos(theta), np.sin(theta), K)
    return padesum(coeff, rho, K // 2)[0]


def padesum2_safe(c, Eps, delta, N, M):
    rho, theta = _polar(Eps, delta)
    K = min(N, M)
    coeff = _ctilde(np.asarray(c), np.cos(theta), np.sin(theta), K)
    return padesum_safe(coeff, rho, K // 2)[0]


def _solve_one(H, rhs):
    try:
        return np.linalg.solve(H, rhs)
    except np.linalg.LinAlgError:
        if SINGULAR_FALLBACK == 'nan':
            return np.full(rhs.shape, np.nan + 0j)
        return np.linalg.lstsq(H, rhs, rcond=None)[0]


# ----------------------------------------------------------------------------
# robust Pade (padeapprox.m, Chebfun / Gonnet-Guettel-Trefethen 2013)
# ----------------------------------------------------------------------------
def padeapprox(f, m, n, tol=1e-14, r=1.0, Nfft=2048):
    """Robust Pade approximation via SVD.

    Returns (r_fun, a, b, mu, nu, poles, residues) as in padeapprox.m.
    f is either a coefficient vector or a callable (coefficients via FFT on |z|=r).
    """
    if callable(f):
        z = r * np.exp(2j * np.pi * np.arange(Nfft) / Nfft)
        fc = np.fft.fft(f(z)) / Nfft
        tc = 1e-15 * np.linalg.norm(fc)
        fc[np.abs(fc) < tc] = 0
        if np.linalg.norm(fc.imag, np.inf) < tc:
            fc = fc.real
        f = fc / r ** np.arange(Nfft)
    f = np.asarray(f).ravel()
    c = np.concatenate([f, np.zeros(max(m + n + 1 - f.size, 0))])[:m + n + 1]
    ts = tol * np.linalg.norm(c)
    if np.linalg.norm(c[:m + 1], np.inf) <= tol * np.linalg.norm(c, np.inf):
        a = np.zeros(1); b = np.ones(1); mu = -np.inf; nu = 0
    else:
        row = np.concatenate([[c[0]], np.zeros(n)])
        col = c
        while True:
            if n == 0:
                a = c[:m + 1]; b = np.ones(1)
                break
            Z = toeplitz(col[:m + n + 1], row[:n + 1])
            C = Z[m + 1:m + n + 1, :]
            rho = int(np.sum(np.linalg.svd(C, compute_uv=False) > ts))
            if rho == n:
                break
            m = m - (n - rho)
            n = rho
        if n > 0:
            _, _, Vh = np.linalg.svd(C, full_matrices=True)
            b = Vh.conj().T[:, n]
            Dm = np.diag(np.abs(b) + np.sqrt(np.finfo(float).eps))
            Q, _ = np.linalg.qr((C @ Dm).conj().T, mode='complete')
            b = Dm @ Q[:, n]
            b = b / np.linalg.norm(b)
            a = Z[:m + 1, :n + 1] @ b
            lam = np.nonzero(np.abs(b) > tol)[0][0]
            b = b[lam:]
            a = a[lam:]
            b = b[:np.nonzero(np.abs(b) > tol)[0][-1] + 1]
        nz = np.nonzero(np.abs(a) > ts)[0]
        a = a[:nz[-1] + 1] if nz.size else a[:0]
        a = a / b[0]
        b = b / b[0]
        mu = a.size - 1
        nu = b.size - 1

    def r_fun(zz, a=a, b=b):
        return polyval_asc(a, zz) / polyval_asc(b, zz)

    poles = np.roots(b[::-1]) if b.size > 1 else np.zeros(0)
    t = max(tol, 1e-7)
    residues = t * (r_fun(poles + t) - r_fun(poles - t)) / 2
    return r_fun, a, b, mu, nu, poles, residues


def padesum_robust(c, Eps, M, tol=1e-14):
    """[M/M] Pade sum computed with padeapprox (SVD based, handles degeneracy)."""
    r_fun = padeapprox(np.asarray(c).ravel()[:2 * M + 1], M, M, tol)[0]
    return r_fun(Eps)


# ----------------------------------------------------------------------------
# vectorised engines (same arithmetic as the scalar versions, batched)
# ----------------------------------------------------------------------------
def _pade_batch(coeff, rho, Mp):
    """Batched padesum: coeff[..., >=2Mp+1], rho broadcastable to coeff.shape[:-1]."""
    rho = np.broadcast_to(rho, coeff.shape[:-1])
    if Mp == 0:
        return coeff[..., 0] + 0j
    B = coeff.shape[:-1]
    cf = coeff.reshape(-1, coeff.shape[-1])
    ii = np.arange(Mp)[:, None] - np.arange(Mp)[None, :]
    H = cf[:, Mp + ii]                             # H[i,j] = c[Mp+i-j]
    rhs = cf[:, Mp + 1:2 * Mp + 1]
    try:
        bb = -np.linalg.solve(H, rhs[:, :, None])[:, :, 0]
    except np.linalg.LinAlgError:
        bb = np.empty_like(rhs)
        nsing = 0
        for k in range(cf.shape[0]):
            try:
                bb[k] = -np.linalg.solve(H[k], rhs[k])
            except np.linalg.LinAlgError:
                nsing += 1
                bb[k] = -_solve_one(H[k], rhs[k])
        warnings.warn(f'padesum: {nsing} exactly singular Toeplitz system(s); '
                      f'used {SINGULAR_FALLBACK} fallback', RuntimeWarning)
    b = np.concatenate([np.ones((cf.shape[0], 1)), bb], axis=1)
    a = np.zeros((cf.shape[0], Mp + 1), dtype=complex)
    for k in range(Mp + 1):
        for i in range(k + 1):
            a[:, k] += b[:, i] * cf[:, k - i]
    x = rho.reshape(-1)
    return (polyval_asc(a, x) / polyval_asc(b, x)).reshape(B)


def _pade_safe_batch(coeff, rho, Mp):
    rho = np.broadcast_to(rho, coeff.shape[:-1])
    big = np.abs(coeff[..., :2 * Mp + 1]) > 1e-14
    last = np.where(big.any(-1), (2 * Mp) - np.argmax(big[..., ::-1], axis=-1) + 1, -1)
    M_safe = np.floor(last / 2.0).astype(int)
    out = np.zeros(coeff.shape[:-1], dtype=complex)
    for ms in np.unique(M_safe):
        if ms < 0:
            continue
        sel = M_safe == ms
        out[sel] = _pade_batch(coeff[sel], rho[sel], int(ms))
    return out


def sum_series(SumType, c, Eps, delta, N, M, taylor_full_order=False):
    """Evaluate the (truncated) double series for a batch of coefficient sets.

    c : array [..., >=N+1, >=M+1]  (index [n, m] = eps^n delta^m)
    Eps, delta : scalars or arrays broadcastable against c.shape[:-2]
    Returns an array of shape broadcast(c.shape[:-2], Eps, delta).
    """
    c = np.asarray(c)[..., :N + 1, :M + 1]
    Eps = np.asarray(Eps, dtype=float)
    delta = np.asarray(delta, dtype=float)
    rho, theta = _polar(Eps, delta)
    K = min(N, M)
    Mp = K // 2
    coeff = _ctilde(c, np.cos(theta), np.sin(theta), K)
    rho = np.broadcast_to(rho, coeff.shape[:-1])
    if SumType == 1:
        return polyval_asc(coeff[..., :(K if taylor_full_order else Mp) + 1], rho)
    if SumType == 2:
        return _pade_batch(coeff, rho, Mp)
    if SumType == 3:
        return _pade_safe_batch(coeff, rho, Mp)
    if SumType == 4:
        out = np.empty(coeff.shape[:-1], dtype=complex)
        for idx in np.ndindex(out.shape):
            out[idx] = padesum_robust(coeff[idx], rho[idx], Mp)
        return out
    raise ValueError('SumType must be 1 (Taylor), 2 (Pade), 3 (Pade safe) or 4 (robust Pade)')


def fcn_sum(SumType, f_n_m, Eps, delta, Nx, N, M, taylor_full_order=False):
    """f_n_m has the interface layout (Nx, M+1, N+1); returns f(x) of length Nx.
    Eps/delta may also be arrays (shape S) -> result shape S + (Nx,)."""
    c = np.transpose(np.asarray(f_n_m)[:Nx, :M + 1, :N + 1], (0, 2, 1))   # (Nx, N+1, M+1)
    Eps = np.asarray(Eps, dtype=float)[..., None]
    delta = np.asarray(delta, dtype=float)[..., None]
    return sum_series(SumType, c, Eps, delta, N, M, taylor_full_order)


def fcn_sum_fast(SumType, f_n_m, Eps, delta, Nx, N, M, taylor_full_order=False):
    """f_n_m has the permuted layout (N+1, M+1, Nx) used in energy_defect.m.
    SumType other than 1/2 falls to Pade-safe (as in MATLAB), except 4 = robust."""
    c = np.moveaxis(np.asarray(f_n_m), 2, 0)                             # (Nx, N+1, M+1)
    st = SumType if SumType in (1, 2, 4) else 3
    Eps = np.asarray(Eps, dtype=float)[..., None]
    delta = np.asarray(delta, dtype=float)[..., None]
    return sum_series(st, c, Eps, delta, N, M, taylor_full_order)


def vol_fcn_sum(SumType, u_n_m, Eps, delta, Nx, Nz, N, M):
    """Volume fields (Nx, Nz+1, M+1, N+1) -> u(x,z) of shape (Nx, Nz+1).
    SumType 1 uses the full double Taylor sum (taylorsum2), as in MATLAB."""
    c = np.transpose(np.asarray(u_n_m)[:, :, :M + 1, :N + 1], (0, 1, 3, 2))  # (Nx,Nz+1,N+1,M+1)
    if SumType == 1:
        vander = np.outer(Eps ** np.arange(N + 1), delta ** np.arange(M + 1))
        return np.sum(c * vander, axis=(-2, -1))
    return sum_series(SumType, c, Eps, delta, N, M)
