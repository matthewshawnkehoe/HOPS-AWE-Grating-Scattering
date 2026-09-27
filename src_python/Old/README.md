# HOPS/AWE Grating Scattering — Python port

A line-by-line Python port of the MATLAB code in `src/` of
[HOPS-AWE-Grating-Scattering](https://github.com/matthewshawnkehoe/HOPS-AWE-Grating-Scattering)
(Kehoe & Nicholls, *J. Sci. Comput.* 100:9, 2024). Every MATLAB function has a Python function
with the same name and argument list, and the array layouts are the same.

## Quick start (PyCharm or a terminal)

```bash
pip install -r requirements.txt
python refl_map.py                        # = refl_map.m  (silver, paper Fig. 10a)   ~6-10 min
python refl_map.py --scenario gold        # paper Fig. 10b
python refl_map.py --scenario dielectric  # = plots/test_scenarios.m, GitHub README, paper Fig. 9
python mms_error.py                       # = mms_error.m (N=M=16, Pade)               ~1-2 min
python mms_error.py --preset fig2 --epsmax 1e-2   # paper Fig. 2(a); try 1e-4, 1e-6, 1e-8
python test_single_eps_delta.py           # = test_single_eps_delta.m (all 7 plot_errors figures)
python test_mms_error.py                  # full validation suite (PASS/FAIL table), ~5 min
python paper_figures.py                   # paper Figs. 2 and 3 (MMS relative error, Taylor)
```
In PyCharm, open the folder as a project and run any of the scripts. The settings at the top of each
script are the same as the variables at the top of the matching `.m` file. Figures are shown and
also saved to `figures/`.

## File map (MATLAB → Python)

| MATLAB | Python |
|---|---|
| `cheb.m`, `dx.m`, `dz.m`, `setup_2d.m` | `hops/spectral.py` |
| `gamma_exp.m`, `E_exp.m`, `E_exp_lf.m`, `A_exp.m`, `T_dno.m` | `hops/expansions.py` |
| `setup_xi_u_nu_u_n_m.m`, `setup_xi_w_nu_w_n_m.m`, `setup_zeta_psi_n_m.m`, `fourier_repr_lipschitz.m`, `fourier_repr_rough.m` | `hops/setup_data.py` |
| `solvebvp_colloc_fast.m`, `solvebvp_colloc_fast_lf.m` (`pagemldivide`) | `hops/bvp.py` |
| `field_tfe_helmholtz_m_and_n(_lf).m`, `dno_tfe_helmholtz_m_and_n(_lf).m` | `hops/fields.py` |
| `AInverse.m`, `two_layer_solve.m`, `two_layer_solve_fast.m` | `hops/two_layer.py` |
| `taylorsum*.m`, `padesum*.m`, `padeapprox.m`, `fcn_sum.m`, `fcn_sum_fast.m`, `vol_fcn_sum.m` | `hops/summation.py` |
| `energy_defect.m` | `hops/energy.py` |
| `plot_errors.m` | `hops/plotting.py` |
| `refl_map.m`, `mms_error.m`, `test_single_eps_delta.m` | `refl_map.py`, `mms_error.py`, `test_single_eps_delta.py` |
| (`test_mms_error.m` — named in the GitHub README but not in the repo; `test_single_eps_delta.m` does that job) | `test_mms_error.py` (validation suite) |

Conventions: MATLAB index `k` is Python `k-1`; fields are `(Nx, Nz+1, M+1, N+1)` and interface data
`(Nx, M+1, N+1)`; **every FFT runs along axis 0**, because MATLAB's `fft` of a matrix works column by column
and numpy's default is the last axis; MATLAB `.'` is numpy `.T`; and `sqrt` of a negative real is
`hops.csqrt` (`numpy.lib.scimath.sqrt`), because `numpy.sqrt` would return NaN.

## How the conversion was checked

The original MATLAB sources were run under GNU Octave 8.4 (with a small `pagemldivide` stand-in), and the
results are saved in `reference_data/`. The scripts that make them are in `reference_data/octave/`, so you
can make the same files with real MATLAB. `test_mms_error.py` / `pytest tests` compare against them:

* Building blocks (cheb, dx, dz, setup_2d, the expansions, T_dno, AInverse, all Taylor/Padé sums, padeapprox):
  they match to 1e-13 to 1e-16.
* Fields u, w, DNOs G, J, and the two-layer U, W, ubar, wbar in the `test_single_eps_delta` setup:
  they match to 1e-12 to 1e-13. The whole (n, m, summation type) error table matches too.
* Reduced `mms_error.m` and `refl_map.m` runs (silver and dielectric): R matches to 1e-9 with Taylor sums.
  With Padé sums the median difference is 1e-8 or less and the worst is 1e-2. D (dielectric) also matches.
* Maths checks: dz is exact on polynomials, fast and slow solvers agree, the manufactured solution is recovered
  to about 1e-12, R + T = 1 for a dielectric, and Padé sums stay stable when the coefficients are degenerate.

### Why the highest-order coefficients differ in the last digits
The HOPS/AWE recursion strengthens round-off from one order to the next. For silver, the |U_{n,m}| grow to
about 1e16 by n = m = 15. Coefficients with n, m ≤ 5 match MATLAB to 1e-8 to 1e-13. The highest orders differ
by O(1) *relative* amounts. Two Python versions of the same code (LU solve and cached inverse) differ from each
other by the same amount as they differ from MATLAB. So this is a property of the floating-point recursion,
not a porting error. The quantities you actually use (R, D, and the summed U, W, ubar) agree.

## About the problems seen in earlier ports

1. **The Chebyshev "sign issue" is not a real problem — do not negate D.** `cheb` returns Trefethen's
   D = d/dt on the nodes t_j = cos(πj/Nz). These nodes run from +1 down to −1, but D is still the derivative
   with respect to t. The layer maps z = (a/2)(t−1)+a and z = (b/2)(t−1) both have dz/dt = L/2 > 0, so
   d/dz = (2/L)·D, which is exactly what `dz.m` computes. `tests::test_dz_exact_derivative` checks this on
   polynomials, and shows that −D gives O(1) errors. The earlier symptoms most likely came from an FFT on the
   wrong axis, from `'` being used where `.'` was meant, or from `numpy.sqrt` returning NaN.
2. **The "Taylor error 1.6e6 for U" does not appear here.** The data follow MATLAB exactly:
   ψ = −ν_u − τ²ν_w, where ν_w already carries the sign for −∂_N w (`setup_xi_w_nu_w_n_m`) and the
   lower DNO J is +∂_N w. U is recovered to about 1e-14 (`test_single_eps_delta.py`, and B3 in the suite).
3. **Padé with ill-conditioned matrices.** MATLAB's `H\c` is a dense LU solve. It warns "matrix singular to
   machine precision" (Octave shows rcond ~ 1e-20 in refl_map) but still returns an answer.
   `numpy.linalg.solve` is the same LAPACK routine and behaves the same way. A Levinson Toeplitz solver
   (`scipy.linalg.solve_toeplitz`) must **not** be used: it fails when a leading minor is singular. If a
   matrix is *exactly* singular, the port uses the minimum-norm SVD solution instead of MATLAB's NaNs
   (set `hops.summation.SINGULAR_FALLBACK = 'nan'` to reproduce MATLAB). A fully robust SVD Padé
   (`padeapprox`, Gonnet–Güttel–Trefethen) is available as **SumType 4**.

## Matching the paper's figures

* Figs. 2/3 (`paper_figures.py`): log10 relative errors come out as [-4, -14.8], [-8, -14.9],
  [-11.8, -14.9], [-14.4, -14.9] for (N=M, eps) = (4,1e-2), (8,1e-4), (12,1e-6), (16,1e-8). The paper
  shows [-4,-14], [-8,-14], [-11.5,-14], [-13.75,-14.15]. The errors do **not** go to ~1e-16 everywhere,
  because the error away from ω = 1.5 comes from truncating the δ-series, not from round-off.
  These numbers only match the paper with MATLAB's half-order Taylor sum (below). With every order summed,
  Fig. 3(a) would reach 1e-6 instead of 1e-4.
* The plots are drawn the way MATLAB does it. Each frequency band is its own `contourf` call with its own
  automatic levels, and all bands share one colour scale ('hot' for R/D, parula for errors). A few isolated
  Padé spikes above R = 1 next to the silver resonance are drawn in the top colour.

## Behaviours kept exactly as in MATLAB (worth knowing)

* `taylorsum_2_coeff.m` sums the polar series only up to ρ^⌊min(N,M)/2⌋. This is kept by default
  (`taylor_full_order=False`). Pass `taylor_full_order=True` to `fcn_sum` / `energy_defect` to use every order.
* In `energy_defect.m`, `alpha_p.^2 < k^2` with complex k compares **real parts** (the MATLAB rule). For metals,
  this means no lower-layer mode counts as propagating, so rl = 0. Octave compares by modulus instead, so D for
  metals differs between MATLAB and Octave. Only R is plotted for metals in the paper.
* With `mms_error.m`'s settings (Nx = 32, f = cos(4x)/4, ε up to 0.2, r = 4), the error in U levels off
  at ~1e-4 for every summation method. This is Fourier truncation/aliasing of e^{iγεf} (mode 4 + 4k reaches
  Nyquist at k = 3), not a summation error. The paper uses Nx = 256 for large ε (`--preset fig6`).

## Speed
The recursions run in about 1–2 minutes per frequency band q (N = M = 15–16, Nx = Nz = 32). The collocation
matrices do not change with the order (n, m), so their inverses are cached
(`hops.bvp.USE_CACHED_INVERSE`; set it to `False` to use a fresh LU solve on every call, as MATLAB does).
All (ε, δ) sums are vectorised.
