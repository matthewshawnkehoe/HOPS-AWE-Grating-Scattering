"""plot_errors.m and small helpers for MATLAB-like contour plots."""
import numpy as np
import matplotlib.pyplot as plt

LATEX_SAFE = {'text.usetex': False}


def safe_log10(x):
    x = np.abs(np.asarray(x, dtype=complex))
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.log10(np.where(x > 0, x, np.nan))


def plot_errors(fig_num, sum_type, N, M, err_G, err_U, err_ubar, err_J, err_W, err_wbar,
                names=('G', 'U', r'\bar{u}', 'J', 'W', r'\bar{w}')):
    """2x3 panel of log10 errors versus (n, m), as plot_errors.m."""
    fig, axs = plt.subplots(2, 3, num=fig_num, figsize=(14, 7), clear=True)
    nn, mm = np.arange(N + 1), np.arange(M + 1)
    for ax, err, nm in zip(axs.ravel(), (err_G, err_U, err_ubar, err_J, err_W, err_wbar), names):
        Z = safe_log10(err).T
        cs = ax.contourf(nn, mm, Z, levels=10)
        fig.colorbar(cs, ax=ax)
        ax.set_xlabel('$n$'); ax.set_ylabel('$m$')
        ax.set_title(f'Error in ${nm}$ ({sum_type})')
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------
# MATLAB look-alike contour plots
# ----------------------------------------------------------------------------
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import MaxNLocator

_PARULA = [(0.2422, 0.1504, 0.6603), (0.2780, 0.3556, 0.9777), (0.1540, 0.5902, 0.9218),
           (0.0265, 0.7160, 0.7865), (0.2178, 0.7964, 0.5982), (0.5709, 0.8032, 0.3453),
           (0.8805, 0.7598, 0.2154), (0.9934, 0.8228, 0.1780), (0.9769, 0.9839, 0.0805)]
parula = LinearSegmentedColormap.from_list('parula', _PARULA, N=256)


def matlab_levels(Z, nbins=10, step=None):
    """Automatic 'nice' contour levels over the finite range of Z (like MATLAB contourf).
    step: fixed level spacing (e.g. step=1 for log10 plots: one level per decade, as in the paper)."""
    Z = np.asarray(Z, dtype=float)
    zmin, zmax = np.nanmin(Z[np.isfinite(Z)]), np.nanmax(Z[np.isfinite(Z)])
    if zmax - zmin < 1e-14 * max(1.0, abs(zmax)):
        return np.array([zmin - 1e-12, zmax + 1e-12])
    if step:
        lev = np.arange(np.ceil(zmin / step) * step, zmax, step)
        lev = lev[(lev > zmin) & (lev < zmax)]
        return np.concatenate([[zmin], lev, [zmax]])
    lev = MaxNLocator(nbins=nbins, steps=[1, 2, 5, 10]).tick_values(zmin, zmax)
    lev = lev[(lev > zmin) & (lev < zmax)]
    return np.concatenate([[zmin], lev, [zmax]])


def matlab_contourf(ax, x, y, Z, cmap, norm, nbins=10, step=None):
    """One MATLAB contourf call: own automatic levels, shared colour limits (CLim)."""
    Z = np.where(np.isfinite(Z), Z, np.nan)
    lev = matlab_levels(Z, nbins, step)
    cs = ax.contourf(x, y, Z, levels=lev, cmap=cmap, norm=norm)
    ax.contour(x, y, Z, levels=lev[1:-1], colors='k', linewidths=0.5, linestyles='solid')   # MATLAB: solid also for negative levels
    return cs


def shared_colorbar(fig, ax, cmap, norm):
    import matplotlib.cm as cm
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    return fig.colorbar(sm, ax=ax)
