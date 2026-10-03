# HOPS/AWE Grating Scattering — Python port

A line-by-line Python port of the MATLAB code in `src/` of
[HOPS-AWE-Grating-Scattering](https://github.com/matthewshawnkehoe/HOPS-AWE-Grating-Scattering)
(Kehoe & Nicholls, *J. Sci. Comput.* 100:9, 2024). Every MATLAB function has a Python function
with the same name and argument list, and the array layouts are the same.

**3D (doubly periodic / crossed gratings):** see [README_3D.md](README_3D.md) — package `hops3d/`,
scripts `refl_map_3D.py`, `mms_error_3D.py`, `test_single_eps_delta_3D.py`, `test_mms_error_3D.py`,
`paper_figures_3D.py`, `material_survey_3D.py`, tests in `tests_3d/`, figures in `figures_3d/`.

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
python paper_figures.py                   # paper Figs. 2-5 (MMS relative error, Taylor)   ~1 min
python paper_figures.py --fig 6 7 8 --workers 6   # Figs. 6-8 at paper resolution (hours; see below)
python paper_figures.py --fig 6 7 8 --quick       # low-resolution preview of Figs. 6-8, ~5 min
```
In PyCharm, open the folder as a project and run any of the scripts. The settings at the top of each
script are the same as the variables at the top of the matching `.m` file. Figures are shown and
also saved under `figures/`, one subdirectory per driver:

| directory | written by |
|---|---|
| `figures/refl_map/` | `refl_map.py` (`refl_map_<tag>_R.png`, `_D.png`, `.npz`) |
| `figures/mms_error/` | `mms_error.py` |
| `figures/paper/` | `paper_figures.py` (Figs. 2-10, 14; cache in `figures/paper/paper_cache/`) |
| `figures/test_single_eps_delta/` | `test_single_eps_delta.py` (plot_errors Figs. 1-6, 11) |
| `figures/material_survey/` | `material_survey.py` |
| `figures/movies/` | `refl_movie.py` |

`figures_3d/` has the same layout for the 3D scripts.

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

## Materials, thesis examples and new scenarios (refl_map.py)
`refl_map.py` now takes any n^u / n^w: a number (`--nw 0.05+2.275i`, `--nw 20i`), a material from the curated
refractiveindex.info library (`--nw Ag`, `--nu water`), or any page of the full database (`--rii-db`). It also
has TE/TM, all thesis Chapter 6 examples (`--scenario thesis17` … `thesis33c`), new material presets, band
q = 0, dispersion per frequency window, and "joint" windows that split at the lower layer's Wood anomalies.
See **MATERIALS.md** for the database review, the survey of 54 materials and the command-line reference.

## Movies (refl_movie.py)
`python refl_movie.py --preset <name>` computes the HOPS/AWE coefficients once, together with the volume
fields u_{n,m} and w_{n,m}, and then animates a sweep in lambda (or in eps). Each frame shows:
- the R map and the D map, with a cursor at the current (lambda, eps);
- the spectrum at the current eps;
- the total field Re H_y (TM) in physical space over two periods: incident + reflected above the grating,
  transmitted below it. Beyond the artificial boundaries the field continues through the exact Rayleigh
  expansions.

Presets:
- `silver_spp`: Ag, P = 0.5 um. Rayleigh anomaly, then the SPP dip at 0.525 um.
- `silver_spp_eps`: eps sweep at the SPP wavelength. Critical coupling at eps ~ 0.1 (R -> 0, near-total
  absorption, |H| ~ 19x the incident field).
- `gold_spr`: Au under water.
- `dielectric`: paper Fig. 9.
- `sic_sphp`: SiC phonon polariton.
- `tir`: glass over air at 53 deg.

Options:
- `--scenario/--lam-range/--eps/--sweep eps --lam`: movie of any other case.
- `--quantity abs`: plot |H| instead of Re H.
- `--format gif`: GIF output. MP4 needs ffmpeg (`pip install imageio-ffmpeg` is enough); otherwise a GIF
  is written.

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

### Figs. 6-8 (large deformations, eps up to 2)
`python paper_figures.py --fig 6 7 8 --workers 6` (default `--res reduced`; `--res paper` uses exactly
(39)/(40)).  All panels now use the coupled solver (`hops/coupled.py`): the two-layer solve takes
~30 s at Nx = 256, Nz = 128 and ~1 min at Nx = 1024 (the lean classic solver took 30 min to 2.5 h), so the
paper resolution (`--res paper`) is now practical; the Pade summation of U on the 100 x 100 grid is
the remaining cost.  Accuracy vs. the lean classic solver (Fig. 6 at paper resolution): the coupled
solver's error is lower in 74-98 % of the cells (mean log10 error -12.47 vs -11.87 for eps > 1 in 6a),
i.e. the same method with slightly less round-off accumulation.  Panels are cached in `figures/paper/paper_cache/`; `--replot` redraws without recomputing.

**How much resolution is needed** (convergence study, relative error of U, N = M = 20, eps <= 2):

| panel | paper (Nx x Nz) | reduced | agreement with paper resolution | time/panel (1 core) |
|---|---|---|---|---|
| 6a f_s, n_w=1.1 | 256 x 128 | **128 x 64** | same picture; ~0.7 digit *lower* error (less round-off) | 31 -> 5 min |
| 6b f_s, n_w=10.1 | 256 x 128 | **128 x 128** | median difference 0.0 digits, 75 % of cells within 0.5 | 31 -> 15 min |
| 7a f_r, n_w=1.1 | 1024 x 128 | **512 x 64** | same as 1024 x 64 (86 % within 0.5 digit, same 95th pct.) | ~2.5 h -> 20 min |
| 7b f_r, n_w=10.1 | 1024 x 128 | (paper) | 512 x 128 not converged | ~2.5 h |
| 8a/b f_L | 1024 x 128 | (paper) | 512 x 128 is ~4 digits worse; Lipschitz needs both | ~1.5-2 h |

Not converged: Nx = 64 for f_s, Nx = 256 for the P = 120 profiles, Nz = 32 or 64 for n_w = 10.1,
Nz = 64 for the Lipschitz profile (its slowly decaying Fourier modes are strongly evanescent and form
thin boundary layers at the interface).  N = M = 20 is needed for eps up to 2.
`--quick` (Nx = 256, Nz = 64, N = M = 12) is only a layout preview.

**Differences from the published figures**: 6a, 7a, 8a match the paper (8a at paper resolution).
For n_w = 10.1 (6b; 7b), the error pattern (fan around omega = 3/2) matches, but at the edges
of the frequency window the errors are 1e-4 to 1e-6 instead of the paper's <= 1e-8.  The classic
MATLAB-equivalent solver gives the same (-4.1), so this is a property of the code for these
parameters (the lower layer, k_w = 15.2, has Rayleigh singularities at |delta| ~ 0.07, far inside
the delta window chosen for the upper layer), not of the Python port.

## Energy defect D: Python vs. MATLAB vs. the paper (Fig. 9)

**The numbers are identical.** The complete Fig. 9 computation was run with the *original, unmodified*
MATLAB code under Octave: `plots/test_scenarios.m` settings, all 6 bands, a 100 x 100 grid per band,
N = M = 16, Taylor summation. That took about 2.5 h of CPU time. The result is saved as
`reference_data/ref_refl_dielectric_full.mat` (script: `reference_data/octave/ref_refl_full.m`).

* Against Python, max |D_Python − D_MATLAB| = 4.7e-13 and max |R_Python − R_MATLAB| = 2e-14 over all
  60,000 (ε, λ) points.
* The ubar coefficients agree to 5e-12 (relative).
* Test A8 checks this. `python paper_figures.py --fig 9` also draws both maps with the same renderer,
  plus their difference: `figures/paper/paper_fig9_matlab_vs_python.png`.

**Why the pictures looked different.** The difference was entirely in how the plots were drawn:

1. **Dashed contour lines.** matplotlib draws the contour lines of *negative* levels dashed; MATLAB draws
   them solid. Since log10|D| is always negative, every line in our D maps was dashed. That is the
   speckled look in the old plots.
2. **Too few levels.** For log10 D, the paper's figure uses one contour level per decade. Our automatic
   levels took 2-decade steps, which halved the number of colour bands.

Both are fixed in `hops.plotting.matlab_contourf`. Lines are now solid everywhere. Maps of D that span at
least 4 decades (lossless layers) get one level per decade (`step=1`). This applies to `refl_map.py`,
`refl_map_3D.py` and `paper_figures.py`, and the D maps now match paper Fig. 9(b) band by band.

For an absorbing layer (silver, gold, W, ...), D is the absorptance 1 − R. It spans less than one decade,
so it keeps MATLAB-style automatic levels. A first version of this fix used decade levels for every D map,
which collapsed the metal absorption maps to one or two colours. That affected only the rendering, never
the numbers.

**Nothing in the computation changed for any material.** `hops/energy.py` and `hops/summation.py` are
untouched, and `--taylor-full-order` and `--windows joint` are opt-in (off by default = MATLAB).

* **Checks against the original MATLAB/Octave code:**
  * Silver (test A5) and gold (new test A9, paper Fig. 10b, reduced grid): R agrees to about 1e-13 in
    band q = 1.
  * Silver, band q = 3 (Taylor and Padé): up to 7e-4 / 1e-2 (median 1e-7). This is the documented
    round-off amplification in silver's highest-order coefficients (see "Why the highest-order
    coefficients differ ..."), not a porting error.
  * Gold, band q = 3: 1e-11 (Taylor) and 2e-9 (Padé).
* **D for metals:** MATLAB treats no lower-layer mode as propagating (its `<` compares real parts), so
  rl = 0 and D = 1 − R. The port reproduces this exactly (test A9).

**Which options help, and for which materials** (100 x 100 grid, 60 x 60 here):

| scenario | summation | effect of `--windows joint` / `--taylor-full-order` |
|---|---|---|
| silver, gold, thesis20a (W), Ag_disp, TiN (absorbing) | Padé | **none: bit-identical R, D.** `joint` only cuts at anomalies of a *lossless* lower layer (none for a metal), and `--taylor-full-order` only changes Taylor sums. |
| dielectric n = 1.1 (Fig. 9) | Taylor | both: median log10 D −8.2 → −14.0; `full-order` alone diverges (+9.7) |
| dielectric n = 1.1 | Padé | `joint`: median −9.5 → −13.0, max −1.8 → −3.6 |
| TiO2 (n ≈ 2.5), P = 1 µm | Padé | `joint` (the scenario default): median −12.3 → −12.6, 99th percentile −4.1 → −4.5 |
| high_contrast (Si, n = 3.48) | Padé | `joint` (the scenario default): 99th percentile −4.1 → −5.7 |
| thesis21 (ZnGeP2, n = 3.19, alpha = 0.01) | Padé | both: median −5.5 → −8.1, but one Padé spike next to an anomaly raises the max to ~1 |

In short: for metals and absorbing materials the options change nothing, and the defaults *are* the
MATLAB algorithm. For lossless dielectrics, `--windows joint` is the option that matters (plus
`--taylor-full-order` when Taylor summation is used). It is a genuine accuracy gain away from the
anomalies, but it is not MATLAB-faithful, so it stays opt-in.

**Why D is only 1e-5 to 1e-8 near the top of the bands** (in the paper as well). This is the method, not
the port:

* **Half-order Taylor sum.** `taylorsum_2_coeff.m` sums the polar series only up to ρ^⌊min(N,M)/2⌋, i.e.
  8 of the 17 orders. Even at ε = 0 this leaves D ≈ (γ̄ a δ)^9 / 9! ≈ 1e-8 at the edges of a band.
* **Branch points inside the bands.** Every band [q, q+1] with n^w = 1.1 contains a lower-layer
  Rayleigh/Wood anomaly (1.1 ω = p). The δ-series has a branch point there, so it cannot converge
  across it.

The table shows the effect over the whole Fig. 9 map (log10 of |D|):

| options | max | median | max for ε ≤ 0.1 |
|---|---|---|---|
| MATLAB default (half order, `--windows paper`) | −4.2 | −8.2 | −5.0 |
| half order, `--windows joint` | −4.4 | −8.9 | −5.0 |
| `--taylor-full-order`, `--windows paper` | **+9.7 (diverges)** | −10.4 | +6.6 |
| `--taylor-full-order --windows joint` | −4.5 | **−14.0** | −5.1 |

Summing all orders across a branch point diverges. So the MATLAB half-order truncation acts as a crude
regulariser. The clean fix is to split each band at the lower-layer anomalies (`--windows joint`) and sum
all orders (`--taylor-full-order`, new option): the median defect drops from 1e-8 to 1e-14. The remaining
maxima (~1e-5) sit right next to the anomalies. The same holds in 3D: `refl_map_3D.py` uses the same half-order default and has the same
`--windows joint` and `--taylor-full-order` options.

## Oblique incidence (alpha != 0): four fixes to the MATLAB formulation
For alpha != 0 the MATLAB code (src.zip and GitHub) is not consistent, which is why the energy defect of
paper Fig. 14 (alpha = 0.01) came out at 1e-5 to 1e-6 instead of the 1e-8 to 1e-15 of Fig. 9. The port now
corrects it (`hops/config.py`, `ALPHA_FIX = True`; set it to `False` to reproduce MATLAB exactly). None of
these fixes changes anything when alpha = 0.
1. `T_dno` was called with p instead of alpha_p = alpha + p, so the first frequency coefficient of gamma_p
   was off by alpha^2/gamma_p. This alone gives D ~ 1e-5 even for a flat interface.
2. The TFE right-hand side dropped the term 2 i alpha A_xz d_z u. The Bloch term 2 i alpha d_x also
   reaches d_z' through the change of variables. This affects both layers.
3. The interface condition needs R -= i alpha(delta) g_x (U - tau^2 W), because the DNOs use the
   phase-removed d_x. `two_layer_solve.m` has this with shifted indices; `two_layer_solve_fast.m` has it
   commented out.
4. `setup_zeta_psi_n_m` kept only alpha_bar of alpha(delta) = alpha_bar (1 + delta) in the psi term.

Checks:
- Manufactured solutions with alpha = 0.1: the fields and DNOs of both layers are now exact to 1e-12
  (before: 1e-5 to 1e-2).
- Two-layer MMS: the errors are the same as for alpha = 0.
- Fig. 14 setting: the median energy defect goes from 10^-5.3 to 10^-8.0, the same as alpha = 0.

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

| `refl_map.py --scenario gold` (6 bands, N=M=15, Nx=Nz=32, 100x100 grid) | time |
|---|---|
| classic `two_layer_solve_fast` (1:1 port of the MATLAB), serial | ~5-10 min |
| operator solver, serial, 1 core | ~27 s |
| operator solver, 6 bands in parallel | ~15 s on 2 cores, less on more |
| **coupled solver** (default everywhere in 2D), 6 bands in parallel | **two_layer_solve 0.5 s per band** (operator: 4.3 s, classic: ~55 s); the run is now dominated by energy_defect (~2.5 s per band) |

**Coupled solver** (`hops/coupled.py`, `SOLVER = 'coupled'`, also used by `solver='auto'`). Same algebra
as two_layer_solve_fast, organised differently: instead of one TFE recursion per data order U_{r,s}
(~(NM)^2/4 order-solves per layer), ONE recursion per layer runs on the full data, interleaved with the
interface solve. At order (n,m) the particular field (zero Dirichlet data) depends only on lower orders
and its DNO *is* sum_{(r,s)<(n,m)} G_{n-r,m-s}[U_{r,s}]; AInverse then gives U_{n,m}, and the exact
homogeneous part e^{i gamma z} U^_{n,m} is added. This works because the recursion is linear and
shift-invariant in (n,m). (N+1)(M+1) order-solves per layer, no Nx x Nx matrices, so it also scales to
the paper's Figs. 6-8: Nx = 256, N = M = 20 in 6.6 s; Nx = 1024, Nz = 128 in 60 s (the lean classic
solver takes ~2 h there). Checked: agrees with two_layer_solve_fast to 1e-10 (B7, alpha = 0 and 0.1),
silver R map (Pade) to 6e-13, and all MATLAB-reference tests pass with it.

What changed (all switchable, all validated against the MATLAB reference data):

1. **Operator formulation** (`hops/operators.py`, `SOLVER = 'operator'`). two_layer_solve_fast
   re-runs the field recursion for every order (q,s), ~(NM)^2/4 = 18,496 order-solves per layer.
   Because U -> G_{p,r}[U] is linear and the same for every (q,s), the port now runs the recursion
   **once** with all Nx Fourier modes e^{i p x} as a batch of right-hand sides (256 order-solves, each a batched
   matrix-matrix product) and stores the Nx x Nx matrices G_{p,r}, J_{p,r} and the trace maps
   U -> ubar, W -> wbar. The two-layer recursion is then just small matrix-vector products. Same maths,
   same result to round-off (~1e-12 on U; R agrees to 1e-13 except a handful of points at the
   ill-conditioned Pade resonance, where MATLAB/Octave also disagree at that level).
   About 19x faster per band. `SOLVER = 'fast'` gives the classic routine; `'auto'` picks the cheaper
   one (the operator form costs about 4*Nx/(N*M) times the flops, so for Nx=256 with N=M=20
   (paper Fig. 6) the classic routine is cheaper).
2. Inside the recursion: derivative calls are fused (dx(A u_x) + dx(B u_z) = dx(A u_x + B u_z));
   the 2i*alpha*d_x terms are skipped when alpha = 0; and only 3 n-levels of the volume field are kept.
3. **Bands in parallel** (`WORKERS = os.cpu_count()`, `--workers N`): the six frequency bands q are
   independent, so they run in separate processes (one BLAS thread each).
4. `energy_defect`: all (eps, delta) sums in one batched call, and the polar coefficients
   ctilde_p are formed with one matrix product per p.
5. BVP matrices do not depend on the order (n,m), so they are inverted once and cached.

Further ideas, if you need more:
* The upper-layer operators (G, ubar map) do not depend on n_w. When you sweep several materials
  (silver, gold, dielectric) at the same q and f, compute them once and reuse them.
* For big grids (N_eps = N_delta = 1000, as in paper Fig. 14), energy_defect is the dominant cost
  (~linear in N_eps*N_delta); run it in parallel over delta chunks or on a GPU (CuPy is a drop-in
  replacement for the numpy calls used).
* For Nx >= 256 the operator form becomes memory-heavy (3*(M+1)*Nx^2*(Nz+1) complex numbers); it
  could be run in chunks of unit vectors, or restricted to the Fourier modes actually excited by f.
