# HOPS/AWE in 3D — doubly periodic (crossed) gratings, two layers

`hops3d/` extends the 2D joint boundary/frequency perturbation method of Kehoe & Nicholls
(*J. Sci. Comput.* 100:9, 2024; package `hops/`) to three dimensions. The interface is
z = g(x, y) = ε f(x, y), with period 2π in both x and y. The frequency is ω = ω̄(1 + δ), and every
quantity is a joint Taylor series

    u(x, y, z; ε, δ) = Σ u_{n,m}(x, y, z) ε^n δ^m,        R(ε, δ) = Σ R_{n,m} ε^n δ^m,

summed with the same polar Taylor/Padé re-summation as in 2D (`hops.summation` is reused).
It is the direct analogue of the 2D code: two layers, scalar Helmholtz, plane-wave incidence
exp(i(αx + βy − γ^u z)), and τ² = (n^u/n^w)² ("TM") or 1 ("TE") in the transmission condition.
It is **not** the three-layer code in `src_thee_layers.zip`, and it is not the vector-Maxwell
formulation of `maxwellfie.anm.pdf`. For crossed gratings, the two polarisations couple in full
Maxwell, so the "TE/TM" names here mean the scalar model only.

The 2D files are untouched. Every 2D directory and script has a 3D twin:

| 2D | 3D |
|---|---|
| `hops/` (package) | `hops3d/` (package; reuses `hops.summation`, `hops.materials`, `hops.plotting`, `cheb`) |
| `tests/` | `tests_3d/` |
| `figures/` | `figures_3d/` (same subdirectories: `refl_map/`, `mms_error/`, `paper/`, `test_single_eps_delta/`, `material_survey/`, `movies/`) |
| `refl_map.py` | `refl_map_3D.py` |
| `mms_error.py` | `mms_error_3D.py` |
| `test_single_eps_delta.py` | `test_single_eps_delta_3D.py` |
| `test_mms_error.py` | `test_mms_error_3D.py` |
| `paper_figures.py` | `paper_figures_3D.py` |
| `material_survey.py` | `material_survey_3D.py` |
| `refl_movie.py` | `refl_movie_3D.py` |

## Quick start

```bash
python test_mms_error_3D.py                      # 12 validation checks, ~10 s
python test_single_eps_delta_3D.py               # all 7 plot_errors figures, 3D manufactured solution
python refl_map_3D.py                            # crossed silver grating, omega in [1, 7], ~50 s
python refl_map_3D.py --scenario dielectric      # crossed analogue of paper Fig. 9, ~6 s
python refl_map_3D.py --scenario silver_1d       # y-invariant: reproduces the 2D silver map, ~4 s
python refl_map_3D.py --nw Au --period 0.8 --profile cosx+cosy --theta 10 --phi 45
python mms_error_3D.py --preset fig2             # MMS error maps (presets matlab, fig2, fig3, fig6, oblique)
python paper_figures_3D.py                       # 3D analogues of Figs. 2-5, ~30 s
python paper_figures_3D.py --fig 6 7 8 --workers 2   # Figs. 6-8 (Nx = Ny = 64), ~4 min per figure
python material_survey_3D.py --keys Ag Au Al TiN TiO2 Si
python refl_movie_3D.py --preset au_crossed_spp         # movie (see below), ~4 min
```

In PyCharm, right-click `test_mms_error_3D.py` or `test_single_eps_delta_3D.py` and choose
*Run 'pytest in …'*, or run them as scripts.

## The method in 3D

**Transformed Field Expansions.** The upper layer uses z' = a(z − g)/(a − g) and the lower layer uses
z' = b(z − g)/(b + g). After multiplying by (a − g)²/a², the equation for ũ = e^{−i(αx+βy)} u is

    Δ'ũ + 2i(α∂x' + β∂y')ũ + γ̄²ũ = F(ε, δ).

The coefficients are the 2D ones with f_x → (f_x, f_y) and |∇f|² in place of f_x² (listed in
`hops3d/layer.py`). There is no x–y cross term. The terms 2i(αA_xz + βA_yz)∂z' are included; these are
the 3D form of fix 4 of the 2D `ALPHA_FIX`, and they are always on in 3D.

**Frequency expansion of the DNO symbol.** γ_pq(δ)² = γ̄_pq² + 2δ(k̄² − ᾱα_p − β̄β_q) + δ²γ̄², so
the 2D recursion (paper eqs. 12–16) carries over with αα_p → αα_p + ββ_q (`T_dno_3d`).

**DNOs.** With N = (−g_x, −g_y, 1), G = −∂_N u and J = +∂_N w (the 2D sign conventions). The phase-removed
DNOs miss i(α g_x + β g_y)(U − τ²W). This is added in the interface equation (`phase_correction_3d`),
as in 2D.

**Incidence.** ζ = −e^{−iγεf} and ψ = (iγ + iα g_x + iβ g_y) e^{−iγεf}, with the full δ-convolution of
α(δ), β(δ) and γ(δ).

**Frequency windows** (`hops3d/windows.py`). In 2D the bands [q, q+1] lie between consecutive Rayleigh
frequencies. In 3D the Rayleigh frequencies are |(α+p, β+q)|/n (for normal incidence: 1, √2, 2, √5, √8,
3, …), and each window between two of them gets its own expansion with |δ| < σ(w_b − w_a)/(w_b + w_a).
This is the same as eq. (30) for [q, q+1].

* `--windows joint` also cuts at the lower-layer anomalies (for lossless n^w), as in 2D.
* `--lattice profile` cuts only at the orders the profile can excite. For example, cos(4x)cos(4y) excites
  only (4a, 4b), so its first cut is at 5.66.
* `--theta/--phi` fix the incidence angles, so α and β scale with ω. The Rayleigh frequencies then come
  from a quadratic.
* `Ny = 1` (or `Nx = 1`) marks an invariant direction and gives back exactly the 2D bands.

**Energy.** R = Σ_prop (γ^u_pq/γ^u_00)|û_pq|² and T = τ² Σ_prop (γ^w_pq/γ^u_00)|ŵ_pq|², summed over the 2D
lattice of propagating orders. D = 1 − R − T.

* `--sum-domain physical` sums ubar(x, y) pointwise and then applies an FFT, exactly like `energy_defect.m`.
* The default `fourier` sums the series of the few propagating amplitudes û_pq directly. It is identical
  for Taylor and agrees for Padé (see the check below). It is also about Nx·Ny/#orders times cheaper.

## Solvers (`hops3d/two_layer.py`)

* **`coupled`** (default, new) runs one TFE recursion per layer on the full data
  U = Σ U_{r,s} ε^r δ^s, interleaved with the interface solve. At order (n, m), the particular field
  (zero Dirichlet data) depends only on lower orders, and its DNO is exactly
  Σ_{(r,s)<(n,m)} G_{n−r,m−s}[U_{r,s}]. AInverse then gives U_{n,m}, and the exact homogeneous part
  e^{iγz}Û_{n,m} is added. Because the recursion is linear and shift-invariant, this is algebraically the
  same as `two_layer_solve_fast.m`. It needs (N+1)(M+1) order-solves per layer instead of ~(N+1)²(M+1)²/4:
  **0.4 s instead of 17 s** for N = M = 10 on a 16×16×17 grid, with identical errors.
  (The same idea would speed up the 2D code as well.)
* **`lean`** follows the `two_layer_solve_fast.m` ordering, accumulating G[U_{r,s}] on the fly.
* **`operator`** uses the Fourier-mode-basis operators, as in 2D. It only pays off for small grids.

## How it was checked (`python test_mms_error_3D.py`, all 12 pass)

| check | result |
|---|---|
| E1 y-invariant gratings (Ny = 1 and 4, α = 0 and 0.1) vs the validated 2D solver | U, W, ubar, wbar agree to 1e-12 (relative) |
| E2 y-invariant R, T vs 2D `energy_defect` | 1e-12 |
| E3 coupled = lean = operator | 1e-11 |
| E4 manufactured solution, egg-crate profile (cos x + cos y + cos x cos y)/3, α = 0.13, β = 0.21, both layers | U, W, G, J, ubar, wbar ~1e-12 to 1e-13 |
| E5 `test_single_eps_delta_3D` (RunNumber 1) — every field, DNO, trace, all three summations | 1e-13 to 1e-15 |
| E6 flat interface (oblique, absorbing) vs Fresnel; R + T = 1 | 1e-14 |
| E7 crossed dielectric, f = cos4x cos4y/4, Nx = Ny = 32, Padé | max \|D\| ≈ 4e-13 |
| E8 x ↔ y symmetry | 1e-14 |
| E9 T_dno_3d vs direct γ_pq(δ) | 1e-12 |
| E10 windows (3D lattice, fixed-angle quadratic, Ny = 1 → 2D bands) | exact |
| E11 Fourier vs physical summation | identical (Taylor); 1e-8 (Padé, clean window) |
| E12 conical incidence (β = 0.3) on a 1D grating — beyond the 2D code | \|D\| < 1e-9 |

A full 3D `refl_map` of the y-invariant silver grating (`--scenario silver_1d --sum-domain physical`) agrees
with the 2D `refl_map.py` silver map to a median of 5e-15. The only larger differences are at the handful of
Padé-sensitive points next to the plasmon resonance, and the maps look the same.

## Paper figures in 3D (`paper_figures_3D.py`)

These use the manufactured mode (p, q) = (4, 0), f_s = cos(4x)cos(4y)/4, and window 1 (ω ∈ [1, √2]).

* **Figs. 2–5** (Nx = Ny = Nz = 32, Taylor) reproduce the 2D behaviour. For example, Fig. 3d has a
  relative error of 1e-14.4 to 1e-14.8.
* **Figs. 6–8** (Padé, a = b = 4) use Nx = Ny = 64, N = M = 16 and eps_max = 1 by default. In 3D the paper's
  Nx = 256/1024 would mean 65k–1M unknowns per Chebyshev level.
* The large-ε error is set by the spatial resolution, just as in 2D, where Nx = 64 was not enough and
  128 was. For f_s with n^w = 1.1, the relative error at ε = 0.5 / 1 / 2 is:

  | Nx = Ny | ε = 0.5 | ε = 1 | ε = 2 |
  |---|---|---|---|
  | 32 | 6e-4 | 5e-3 | 2e-2 |
  | 48 | 2e-6 | 5e-5 | 8e-4 |
  | 64 | 3e-9 | 3e-7 | 2e-5 |

* n^w = 10.1 also needs a finer z-grid, as in 2D, where Nz = 128 was used.
* `--eps-max 2 --nx 96` approaches the paper setting (about 4 GB of RAM per panel).

## Movies (`refl_movie_3D.py`)

Like `refl_movie.py`, the coefficients are computed once per frequency window, and every frame is only a
Taylor/Padé summation. To make the field pictures, the coupled solver keeps two things: the interface
data U, W, and one vertical slice (y = y0) of the volume fields u_{n,m} and w_{n,m}.

Each frame has six panels:

* **R map and D/A map.** Each has a cursor at the current (λ, ε).
* **Spectrum.** R, T and A at the current ε.
* **x–z field cut.** Re u_tot through two periods: incident + reflected above the grating and transmitted
  below it. Beyond z = a and z = −b, the field is continued with the exact Rayleigh expansion.
* **Surface top view.** |u_tot| on the crossed surface z = εf(x, y) over 2 × 2 periods. This is the 2D
  hot-spot / standing-wave pattern that a 1D grating cannot produce.
* **Diffraction orders in reciprocal space.** Each order (α+p, β+q) is drawn as a dot, together with the
  light circles n^u ω (solid) and n^w ω (dashed). As λ decreases, the orders enter the circles at the
  Rayleigh anomalies. Black discs show the reflected efficiencies R_pq and blue rings the transmitted
  efficiencies T_pq, with marker area proportional to √efficiency.

| preset | what it shows |
|---|---|
| `au_crossed_spp` | Gold cos x cos y, P = 0.8 µm. The (±1,±1) orders become evanescent at 0.566 µm and couple to surface plasmons at about 0.60 µm: R dips to about 0.35 and the field on the crests reaches about 7× the incident field. The (±1, 0) orders are propagating but never lit, because cos x cos y cannot excite them. |
| `au_crossed_spp_eps` | The same at λ = 0.60 µm, with ε growing from 0 to 0.2. |
| `silver_crossed` | Silver cos4x cos4y (nondimensional), with windows cut only at the excited orders. The (4,4)-order plasmon appears at λ ≈ 1.23 (2π/5.14). |
| `egg_silver` | Silver egg-crate profile: a sharp plasmon line at λ ≈ 4.9, and several orders active at once. |
| `dielectric_crossed` | n = 1.1 with joint windows. R ≈ 0.002, so R and T are shown on twin axes. Wood anomalies appear where orders enter the two light circles; the energy defect stays around 1e-10. |
| `conical` | A y-invariant cos x grating under conical incidence (β = 0.3). The orders lie on the line β + q = 0.3, off the p-axis. This case is beyond the 2D code. |
| `oblique_gold` | Gold at θ = 20°, φ = 30°. The order pattern is asymmetric, and the surface pattern is tilted. |

The gold spectra look like staircases because the index is held constant within each window (see
`--max-delta`), exactly as in 2D. The jumps shrink as the windows get smaller.

## Energy defect in 3D (same story as 2D)

The 2D Fig. 9 investigation (README.md, "Energy defect D: Python vs. MATLAB vs. the paper") applies here
unchanged:

* D maps are now drawn with one contour level per decade and solid lines, like MATLAB.
* The size of D is set by two things: the MATLAB-style half-order Taylor sum, and Rayleigh anomalies of
  the lower layer inside a window.

`refl_map_3D.py --scenario dielectric` (cos x cos y, n^w = 1.1, 50 x 50 grid per window), log10 of |D|:

| options | max | median |
|---|---|---|
| default (half order, `--windows paper`) | +0.2 | −8.5 |
| half order, `--windows joint` | −4.4 | −9.3 |
| `--taylor-full-order`, `--windows paper` | +11 (diverges) | −10.7 |
| `--taylor-full-order --windows joint` | −4.6 | **−13.6** |

In 3D the paper-style windows contain more lower-layer anomalies (e.g. |(1,1)|/1.1 = 1.286 inside
[1, √2]), so `--windows joint` matters even more than in 2D.
