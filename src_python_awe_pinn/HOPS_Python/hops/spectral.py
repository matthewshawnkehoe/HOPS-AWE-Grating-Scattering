"""Spectral building blocks: cheb.m, dx.m, dz.m, setup_2d.m.

Array conventions (used throughout the package)
-----------------------------------------------
* Arrays keep the MATLAB layouts, e.g. ``u_n_m`` is ``(Nx, Nz+1, M+1, N+1)``
  and interface data are ``(Nx, M+1, N+1)``; MATLAB index ``k`` is Python ``k-1``.
* MATLAB ``fft`` on a matrix acts along the FIRST dimension, so every FFT here
  uses ``axis=0``. (numpy's default ``axis=-1`` is a classic porting bug.)
* MATLAB ``.'`` (non-conjugate transpose) is numpy ``.T``; MATLAB ``'`` on a
  complex array would be ``.conj().T`` -- the HOPS code only uses ``.'`` on
  complex data.

About the Chebyshev "sign issue"
--------------------------------
``cheb`` returns exactly Trefethen's matrix: ``D[i, j] = d l_j / d t (t_i)`` on
the nodes ``t_i = cos(pi*i/Nz)``, which run from +1 down to -1. It is the
derivative with respect to ``t`` itself (not with respect to the node index),
so *no* extra minus sign is needed.  The physical map used by HOPS is
``z = (a/2)(t - 1) + a`` (upper layer, z in [0, a]) or
``z = (b/2)(t - 1)`` (lower layer, z in [-b, 0]); in both cases
``dz/dt = L/2 > 0``, hence ``d/dz = (2/L) d/dt`` -- exactly what dz.m does.
Adding a minus sign breaks the solver; this is verified in
``tests/test_against_matlab.py::test_dz_exact_derivative``.
"""
import numpy as np


def cheb(N):
    """CHEB  compute D = differentiation matrix, x = Chebyshev grid (Trefethen)."""
    if N == 0:
        return np.zeros((1, 1)), np.ones(1)
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * (-1.0) ** np.arange(N + 1)
    X = np.tile(x[:, None], (1, N + 1))
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))   # off-diagonal entries
    D = D - np.diag(D.sum(axis=1))                    # diagonal entries
    return D, x


def _col(p, ndim):
    """Reshape a length-Nx vector so it broadcasts along axis 0 of an ndim array."""
    return np.reshape(p, (-1,) + (1,) * (ndim - 1))


def dx(u, p):
    """u_x = ifft((1i*p).*fft(u))   (spectral x-derivative along axis 0)."""
    u = np.asarray(u)
    return np.fft.ifft(1j * _col(p, u.ndim) * np.fft.fft(u, axis=0), axis=0)


def dz(u, Dz, b):
    """u_z = ((2.0/b)*Dz*u.').'   for u of shape (Nx, Nz+1)."""
    return (2.0 / b) * (u @ Dz.T)


def csqrt(x):
    """MATLAB-style sqrt: sqrt of a negative real returns +i*sqrt(|x|)."""
    return np.lib.scimath.sqrt(x)


def setup_2d(Nx, L, alpha, beta):
    """[xx,kk,alphap,betap,eep,eem] = setup_2d(Nx,L,alpha,beta).

    betap is chosen with Im{betap} >= 0 (outgoing / bounded branch).
    """
    xx = (L / Nx) * np.arange(Nx)
    kk = (2.0 * np.pi / L) * np.concatenate([np.arange(0, Nx // 2), np.arange(-Nx // 2, 0)])
    alphap = alpha + kk
    kappa = csqrt(alpha ** 2 + beta ** 2)
    value = kappa ** 2 - alphap ** 2
    rho = np.abs(value)
    theta = np.angle(value)
    theta = np.where(theta < 0, theta + 2 * np.pi, theta)
    betap = np.sqrt(rho) * np.exp(1j * theta / 2.0)
    eep = np.exp(1j * alpha * xx)
    eem = np.exp(-1j * alpha * xx)
    return xx, kk, alphap, betap, eep, eem
