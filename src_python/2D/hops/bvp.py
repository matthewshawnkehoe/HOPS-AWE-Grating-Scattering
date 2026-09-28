"""Batched Chebyshev-collocation two-point BVP solvers:
solvebvp_colloc_fast.m and solvebvp_colloc_fast_lf.m.

Solves, for every Fourier mode p (batched, replaces MATLAB ``pagemldivide``):
    alpha*U'' + beta*U' + gamma_p*U = b_p   on the Chebyshev grid,
with Robin conditions on the first row (top, t=+1) and the last row (bottom, t=-1).
"""
import numpy as np


# Speed option: the collocation matrices depend only on the layer (gamma_p, D, the
# boundary operators) and NOT on the perturbation order (n,m), so they can be
# inverted once and re-used for every order.  This gives a ~2x speed-up; results
# agree with the plain LU solve (MATLAB pagemldivide) to ~1e-13.  Set to False to
# use a fresh LU solve on every call, exactly as MATLAB does.
USE_CACHED_INVERSE = True
_CACHE = {}


def bvp_matrix(alpha, beta, gamma, Nx, identy, D, D2, D_start, D_end, top_n, top_d, bot_n, bot_d):
    A = (alpha * D2 + beta * D)[None, :, :] + np.asarray(gamma).reshape(Nx, 1, 1) * identy[None, :, :]
    A = A.astype(complex)
    A[:, -1, :] = bot_n * D_end
    A[:, 0, :] = top_n * D_start
    A[:, -1, -1] += bot_d
    A[:, 0, 0] += top_d
    return A


def _inverse(A):
    key = (A.shape, A.tobytes())
    Ainv = _CACHE.get(key)
    if Ainv is None:
        if len(_CACHE) > 64:
            _CACHE.clear()
        Ainv = np.linalg.inv(A)
        _CACHE[key] = Ainv
    return Ainv


def _solve(b, alpha, beta, gamma, Nx, identy, D, D2, D_start, D_end,
           top_n, top_d, top_r, bot_n, bot_d, bot_r):
    A = bvp_matrix(alpha, beta, gamma, Nx, identy, D, D2, D_start, D_end, top_n, top_d, bot_n, bot_d)
    b = np.array(b, dtype=complex, copy=True)          # (Nz+1, Nx)
    b[-1, :] = bot_r
    b[0, :] = top_r
    if USE_CACHED_INVERSE:
        return np.einsum('pij,jp->pi', _inverse(A), b)   # (Nx, Nz+1)
    # b.T[..., None] -> (Nx, Nz+1, 1); solve page by page (LU, like MATLAB '\\')
    return np.linalg.solve(A, b.T[:, :, None])[:, :, 0]


def solvebvp_colloc_fast(b, alpha, beta, gamma, d_min, n_min, r_min, d_max, n_max, r_max,
                         Nx, identy, D, D2, D_start, D_end):
    """Upper layer: Dirichlet (d_min,n_min) at bottom row, DNO-Robin (d_max vector) on top.
    Returns Uhat of shape (Nx, Nz+1)."""
    return _solve(b, alpha, beta, gamma, Nx, identy, D, D2, D_start, D_end,
                  top_n=n_max, top_d=np.asarray(d_max), top_r=r_max,
                  bot_n=n_min, bot_d=d_min, bot_r=r_min)


def solvebvp_colloc_fast_lf(b, alpha, beta, gamma, d_min, n_min, r_min, d_max, n_max, r_max,
                            Nx, identy, D, D2, D_start, D_end):
    """Lower layer: DNO-Robin (d_min vector) at bottom row, Dirichlet on top.
    Returns What of shape (Nx, Nz+1)."""
    return _solve(b, alpha, beta, gamma, Nx, identy, D, D2, D_start, D_end,
                  top_n=n_max, top_d=d_max, top_r=r_max,
                  bot_n=n_min, bot_d=np.asarray(d_min), bot_r=r_min)
