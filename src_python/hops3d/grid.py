"""3D (doubly periodic) spectral building blocks: the analogues of setup_2d.m, dx.m, dz.m.

Geometry
--------
Two layers separated by the doubly periodic interface  z = g(x, y) = eps * f(x, y),
periods d_x (in x) and d_y (in y), default 2*pi each.  The unknowns carry the Bloch phase
removed:  u(x,y,z) = e^{i(alpha x + beta y)} u~(x,y,z), u~ periodic.

Array conventions (3D analogues of the 2D package ``hops``)
-----------------------------------------------------------
* interface data / DNOs / traces:  (Nx, Ny, M+1, N+1)            [:, :, m, n]
* volume fields:                   (Nx, Ny, Nz+1, M+1, N+1)      [:, :, ell, m, n]
* Fourier modes (p, q) = (kx[j], ky[k]) are on axes 0 and 1 (numpy/MATLAB FFT order).
* Chebyshev index ell = 0 is the top of a layer (t = +1), ell = Nz the bottom (t = -1);
  d/dz = (2/L) d/dt, exactly as in the 2D code (no extra sign, see hops/spectral.py).
"""
import numpy as np
import scipy.fft as sfft

from hops.spectral import cheb, csqrt          # noqa: F401  (re-exported)

fft2 = lambda u: sfft.fft2(u, axes=(0, 1))
ifft2 = lambda u: sfft.ifft2(u, axes=(0, 1))


def wavenumbers(N, L):
    """(2 pi / L) * [0, 1, ..., N/2-1, -N/2, ..., -1]   (for N = 1: [0])."""
    if N == 1:
        return np.zeros(1)
    return (2.0 * np.pi / L) * np.concatenate([np.arange(0, N // 2), np.arange(-(N // 2), 0)])


def outgoing_sqrt(value):
    """sqrt with Im >= 0 (outgoing/decaying branch), the rule of setup_2d.m."""
    value = np.asarray(value, dtype=complex)
    theta = np.angle(value)
    theta = np.where(theta < 0, theta + 2 * np.pi, theta)
    return np.sqrt(np.abs(value)) * np.exp(1j * theta / 2.0)


def setup_3d(Nx, Ny, Lx, Ly, alpha, beta, gamma):
    """3D analogue of setup_2d.

    Returns xx (Nx,), yy (Ny,), kx (Nx,), ky (Ny,), alphap (Nx,), betap (Ny,), gammap (Nx, Ny)
    with  alphap = alpha + kx,  betap = beta + ky,  k^2 = alpha^2 + beta^2 + gamma^2,
    gammap = sqrt(k^2 - alphap^2 - betap^2),  Im gammap >= 0.
    """
    xx = (Lx / Nx) * np.arange(Nx)
    yy = (Ly / Ny) * np.arange(Ny)
    kx = wavenumbers(Nx, Lx)
    ky = wavenumbers(Ny, Ly)
    alphap = alpha + kx
    betap = beta + ky
    k2 = alpha ** 2 + beta ** 2 + gamma ** 2
    gammap = outgoing_sqrt(k2 - alphap[:, None] ** 2 - betap[None, :] ** 2)
    return xx, yy, kx, ky, alphap, betap, gammap


def dx3(u, kx):
    """d/dx along axis 0 (spectral)."""
    sh = (-1,) + (1,) * (np.ndim(u) - 1)
    return sfft.ifft(1j * kx.reshape(sh) * sfft.fft(u, axis=0), axis=0)


def dy3(u, ky):
    """d/dy along axis 1 (spectral)."""
    sh = (1, -1) + (1,) * (np.ndim(u) - 2)
    return sfft.ifft(1j * ky.reshape(sh) * sfft.fft(u, axis=1), axis=1)


def dz3(u, Dz, L):
    """d/dz along the Chebyshev axis 2 of u (Nx, Ny, Nz+1, ...):  (2/L) Dz u."""
    u = np.asarray(u)
    ut = np.moveaxis(u, 2, -2) if u.ndim > 3 else u[..., None]
    out = (2.0 / L) * np.matmul(Dz, ut)
    return np.moveaxis(out, -2, 2) if u.ndim > 3 else out[..., 0]


def grad_profile(f, kx, ky):
    """Spectral (f_x, f_y) of a real periodic profile f (Nx, Ny)."""
    fh = sfft.fft2(f)
    fx = np.real(sfft.ifft2(1j * kx[:, None] * fh))
    fy = np.real(sfft.ifft2(1j * ky[None, :] * fh))
    return fx, fy
