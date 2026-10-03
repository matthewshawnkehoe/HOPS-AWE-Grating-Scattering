# 3D HOPS/AWE + PINN: every 3D scenario, with an error analysis

This directory re-runs the **3D HOPS/AWE code** (`Hops_3D_AWE` = the 3D part of `HOPS_Python`) with the combined
HOPS/AWE + PINN method that did best in 2D. The 3D code covers crossed and doubly periodic gratings and comprises
`hops3d/`, `refl_map_3D.py`, `refl_movie_3D.py`, `material_survey_3D.py` and `paper_figures_3D.py`.

The re-run consists of:

- every reflectivity-map / energy-defect scenario;
- the movies;
- an error analysis against converged references;
- a survey of the material library.

It is the 3D twin of `../HOPS_AWE_PINN`.

## 1. The method and why it was chosen

The 2D studies compared a PINN alone, PINN correction of the Padé sum, neural operators and the **physics-informed
summation (PI-sum)**. PI-sum was the clear winner, so it is the method used here.

PI-sum works as follows:

1. HOPS/AWE computes the coefficient fields u_{n,m}, w_{n,m} of the joint (ε, δ) expansion once per frequency
   window.
2. Those fields are used as the hidden layer of a least-squares interface PINN.
3. At each (ε, ω) the output weights come from **one linear least-squares solve of the PINN loss** at that point:
   - the Helmholtz equation in both layers;
   - the two interface conditions with the exact incident data;
   - the exact 2D-FFT transparent (DtN) conditions at z = a and z = −b.
4. The method is adaptive. The PINN loss of the full-order AWE Taylor sum is a cheap indicator. Where it is below
   tol = 1e-14, the Taylor sum is kept; elsewhere the least-squares weights are used.

Code: `awepinn3d/core.py`. Its `run(scenario, ...)` is a drop-in for `refl_map_3D.run`: the AWE numbers are exactly
those of `refl_map_3D.py`, and the hybrid values are added next to them.

### What had to change for 3D

**1. The residual must be HOPS's own discrete equation.** The 2D hybrid wrote the Helmholtz rows in physical
coordinates via the chain rule. That form is analytically identical to the TFE equation. On an under-resolved
product grid, however, it differs from hops3d's *discrete* equation, which uses conservative flux form with products
formed on the grid (`layer.py::_rhs_hat`). An example is cos 4x cos 4y on 32 × 32, where order-n fields contain the
modes 4n.

The two forms differ by aliasing errors of about 1e-8 at ε = 0.05, so the least-squares step was fitting a slightly
different problem. The interior rows are now **exactly the discrete equation hops3d solves**, summed over all
orders. With that form, the AWE Taylor sum has a residual at round-off level.

| scenario (summation error, order 12) | chain-rule rows (2D form) | discrete rows (3D default) |
|---|---|---|
| silver, cos4x cos4y, 32 × 32 | 1.9e-2 | **9.8e-6** |
| silver_lattice | 2.8e-2 | **1.3e-5** |
| silver_1d (= the 2D silver map) | 1.3e-2 | **6.9e-5** |
| gold | 1.1e-4 | **1.6e-7** |
| sum_gold | 9.9e-6 | **2.5e-9** |
| dielectrics, materials with cos x cos y | the same | the same |

The chain-rule value for silver_1d, 1.3e-2, is exactly the 2D hybrid's error on the 2D silver map, which the 2D
report attributed to "discretisation". **That was aliasing in the residual. The same fix applied to the 2D hybrid
should remove most of its metal errors on Nx = 32.**

**2. Exact compression of the PDE rows.** After multiplying by the layer thickness, the discrete equation is
polynomial in ε (degree 2) and in δ (degree 2):

  A_pde(ε, δ) = Σ_{k,j ≤ 2} ε^k δ^j E_kj

The column-stacked [E_00 … E_22] is QR-factorised once per window. As a result, the Nx·Ny·(Nz−1) PDE rows of every
map point collapse *exactly* to the small triangular factor, and only the 4·Nx·Ny boundary rows are assembled per
point. The cost is 0.05–0.4 s per least-squares point, on 16 × 16 × 16 and 32 × 32 × 16 grids respectively.

**3. A reduced basis does not work.** Building a reduced basis from snapshot solutions was tried. The optimal
weights do not lie on a low-dimensional manifold, and the result had 1e-6 errors, so it was dropped.

## 2. How to run

The code needs `../HOPS_Python` as a sibling directory, or the environment variable `HOPS_PYTHON`. It does not need
PyTorch.

```bash
pip install -r requirements.txt
python -m pytest tests -q                       # 4 tests (~3 min)
python run_scenarios.py                         # all 19 refl_map_3D scenarios + silver_1d_nx64 (~6 h, 2 processes)
python run_scenarios.py --only dielectric gold  # subset;  --force recomputes;  --replot
python error_analysis.py                        # comparison with converged pointwise HOPS (~3 h)
python error_analysis.py --report-only
python materials_hybrid_3D.py --gallery         # material library survey (~2 h)
python hybrid_movie_3D.py --only conical        # one movie per call is safest (each holds its windows in memory)
```

**One reflectivity map (like `refl_map_3D.py`):** `refl_map_hybrid_3D.py` runs one scenario with both methods
and plots it. In PyCharm, edit `SCENARIO` at the top and press Run, or:

```bash
python refl_map_hybrid_3D.py                                # dielectric, 15 x 15 per window
python refl_map_hybrid_3D.py --scenario gold --q 1 2        # first two bands only
python refl_map_hybrid_3D.py --nw Au --period 0.8 --profile cosx+cosy --q 1
```

Output goes to `figures/refl_map_hybrid/`.

```python
from awepinn3d import core
res, info = core.run('gold', qq=(1,), N_Eps=15, N_delta=15)   # refl_map_3D format, hybrid ru/rl/ee
```

## 3. Scenarios reproduced (`results/scenarios.md`, `figures/<scenario>/`)

All 19 `refl_map_3D.py` scenarios were run, plus `silver_1d_nx64`, which is silver_1d on a resolved x-grid. Each
scenario has:

- `refl_map_3D_<s>_awe_{R,D}.png` and `_hybrid_{R,D}.png`: the `refl_map_3D.plot` rendering;
- `compare.png`: both maps, |R_hyb − R_AWE|, both energy-defect maps, and the indicator with the corrected points.

The maps use 21² points per window for up to 8 windows, 15² for up to 25 windows, and 9² for SiC (97 windows),
against 100² in `refl_map_3D.py`.

**Energy defect, max log10|D|, for the lossless scenarios**

| scenario | AWE | hybrid |
|---|---|---|
| dielectric (paper-style windows, MATLAB half-order Taylor) | +0.2 | −10.2 |
| dielectric_joint | −4.4 | −9.4 |
| dielectric_alpha | +0.9 | −6.2 |
| dielectric_1d | −4.2 | −10.5 |
| conical_1d | −4.8 | −9.2 |
| TiO2_crossed | −2.9 | −8.2 |

`dielectric_1d` (y-invariant) reproduces the 2D hybrid's dielectric numbers exactly: 67 % of points corrected, max
|ΔR| 6.3e-5, and an energy defect going from 10^-4.2 to 10^-10.5. This cross-validates the 3D code against the 2D
code.

**Metals and materials.** max |R_hyb − R_AWE| is 3.1e-2 for silver, 1.1e-2 for gold, 8.3e-2 for silver_lattice,
2.7e-2 for water_over_gold_crossed and 1.7e-3 for Si_crossed. Section 4 shows which method is right at those
points.

## 4. Error analysis (`results/error_analysis.md`, `figures/error_analysis/`)

The analysis uses **2,160 samples from 20 scenarios**: in every window, δ = −0.8 δ_max and +0.45 δ_max, times
ε = 0.15, 0.4, 0.65 and 1.0 ε_max. Two references are used:

- **Total error:** a converged reference, namely HOPS solved pointwise in frequency (δ = 0, no frequency
  expansion), N = 24, Padé in ε, on 2Nx × 2Ny × (Nz + 16).
- **Summation error:** the same on the scenario's own grid. This isolates the only step the PINN changes.

| | 3D HOPS/AWE (refl_map_3D) | 3D HOPS/AWE + PINN |
|---|---|---|
| max summation error in R | 5.9e-3 | **6.9e-5** |
| median summation error in R | 7.8e-11 | **2.8e-13** |
| max total error in R | 5.8e-3 | 1.2e-3 (metals on 32 × 32: grid-limited) |
| median total error in R | 1.1e-10 | **4.7e-13** |
| indicator vs true error | – | Spearman ρ = 0.91 |

**Per family (summation error, max)**

| family | AWE | hybrid |
|---|---|---|
| crossed dielectrics (dielectric, _joint, _alpha, _1d, conical_1d) | 1.4e-6 … 3.4e-5 | 9e-13 … 2.3e-10, better at 81–96 % of samples |
| crossed silver / gold, cos4x cos4y, 32 × 32 | 1.6e-4 … 5.9e-3 | 1.6e-7 … 1.3e-5 |
| silver_1d → silver_1d_nx64 | 3.7e-3 → 2.2e-3 | 6.9e-5 → **2.1e-6** |
| named materials (Au, Ag, Al, water/Au, TiO2, oblique gold) | 1.4e-6 … 2.3e-4 | 4e-9 … 1.2e-7 |
| egg_silver | 1.2e-4 | 5.9e-8 |
| Si_crossed | 1.4e-5 | **4.5e-5 (worse)** |
| SiC reststrahlen | 5.7e-10 | 4.8e-8 (both negligible) |

Further findings:

- **Metals on 32 × 32 are limited by the grid, not the summation.** The summation error drops 100–1000×, but both
  methods differ from the refined reference by about 1e-3, so the total error is set by Nx. Silver_1d at Nx = 64
  shows that the hybrid then tracks the refinement (2.1e-6) while AWE does not (2.2e-3).
- **Where the hybrid is worse**, after excluding round-off (both errors below 1e-12 at 27 % of samples), it is
  worse by more than 2× at 11 % of samples. This is almost always where Padé is already at the 1e-10 level and the
  least-squares weights stop at a floor of about 1e-8 (Au, Ag, oblique gold). The worst such case is 5.8e-7
  (silver).
- **Si_crossed is the one real loss.** The error does not depend on ε, so it lies in the δ direction. The substrate
  is weakly absorbing and has a high index, with near-grazing orders of the substrate inside the windows, and the
  basis does not converge with order (5e-5 at every order).

  Where AWE and the hybrid disagree most on the Si map (the λ ≈ 0.71 µm resonance), the pointwise reference shows
  AWE wrong by 1.7e-3 and the hybrid right to 6e-11.
- **Tolerance sweep.** Changing tol from 1e-20 to 1e-8 changes the fraction of points solved from 76 % to 30 %. The
  maximum error does not change, and the median error rises slowly from 6.7e-14 to 4.2e-12.
- **Basis order.** Errors fall from order 6 to 12: silver 7.5e-3 → 9.8e-6, gold 4e-4 → 1.6e-7, dielectric_1d
  8.7e-8 → 1.2e-14, Al 5.9e-7 → 4.1e-12.

## 5. Movies (`figures/movies/`)

| file | content |
|---|---|
| `hybrid_au_crossed_spp.mp4`, `hybrid_au_crossed_spp_eps.mp4`, `hybrid_silver_crossed.mp4`, `hybrid_egg_silver.mp4`, `hybrid_dielectric_crossed.mp4`, `hybrid_conical.mp4`, `hybrid_oblique_gold.mp4` | the seven `refl_movie_3D.py` presets, re-made with the hybrid (see below) |
| `compare_dielectric_eps.mp4` | **new.** Crossed dielectric, ε swept from 0 to 0.4 (twice the map range). The AWE energy defect grows to 10^-5.5 and its R error to about 1e-7; the hybrid stays at 10^-14 to 10^-12.3 in |D| and below 1e-13 in R. |
| `compare_dielectric_lambda.mp4` | **new.** Joint windows, λ swept through the Rayleigh anomalies of both layers at ε = 0.2. AWE |D| reaches 10^-4.5 at the anomalies; the hybrid stays near 1e-14 and stays on the truth curve. |
| `compare_TiO2_lambda.mp4` | **new.** Lossless high-index TiO2 crossed grating (period 1 µm), λ from 1.0 to 0.51 µm at ε = 0.15. AWE |D| spikes to 10^-3.5 at λ ≈ 0.7 and 0.9 µm; the hybrid stays at 10^-8 or below. The staircase in R is the per-window refractive index (dispersion), as in `refl_movie_3D.py`. |

In the preset movies, the maps, the spectra, the diffraction-order efficiencies and the field pictures (x–z cut and
crossed-surface top view) all come from the physics-informed weights. `refl_movie_3D.make_movie` is reused unchanged,
with only its module and field renderer swapped.

The comparison movies are new. Each frame shows:

- R along the sweep: AWE, the hybrid and pointwise HOPS (truth, same grid);
- |D| for both methods;
- |R − R_true| for both methods, together with the indicator;
- the hybrid |u| on the crossed surface;
- log10|u_hyb − u_AWE| on the surface;
- the hybrid x–z field.

## 6. What the materials survey shows (`results/materials.md`, `figures/materials/`)

The survey covers 54 materials from `hops/materials.py`, all except the superstrates. The setup matches
`material_survey_3D.py`:

- vacuum above a crossed cos x cos y grating;
- Padé summation, N = M = 10, on a 16 × 16 × 16 grid;
- joint windows, ω ∈ [1, 3], with the material's own period;
- 11 × 11 points per window.

For every material, pointwise HOPS on the edge columns of every window gives the **true** summation error of both
methods.

- **The hybrid is never worse than AWE for any material.**
  - The median of the per-material maximum true error falls from 6.4e-3 (AWE) to 2.3e-9 (hybrid).
  - For 49 of the 54 materials the hybrid is below 1e-6 everywhere sampled; AWE is below 1e-6 for none of them.
  - The worst AWE cases are Mo (0.59), W (0.29), Ag_IR and Au_IR (0.16), and Pt (0.11).
- **Absorption decides how much AWE needs the correction**, as in 2D.
  - Spearman ρ with the true AWE error is 0.92 for Im n, −0.78 for Re ε, −0.69 for being lossless, and −0.31 for
    Re n.
  - The surface-plasmon distance |Re ε + 1| is not predictive (ρ = −0.16).
- **Rule of thumb** (a depth-2 decision tree, 96 % training accuracy): if Re ε ≤ −2.9 (a good metal), the 3D AWE is
  off by more than 0.01 in R somewhere on the map. Otherwise it is not.
- **Group medians of the true maximum error:**
  - lossless dielectrics: 3.3e-4 → 1.6e-9;
  - absorbing with Re ε > 0: 2.1e-3 → 2.2e-8;
  - metallic: 4.1e-2 → 2.3e-9.
- **Lossless dielectrics** (TiO2, Si3N4, diamond, GaN, AlN, ZnS, LiNbO3, sapphire, SiO2, …) get many joint
  windows (13–21). Their energy defect goes from 10^-2.3 to 10^-3.7 with AWE to 10^-6.7 to 10^-8.8 with the
  hybrid.
- **Limits.** The PEC-like IR metals (Ag_IR, Au_IR, n ≈ 36i) improve only 0.16 → 7e-3: the basis is already at
  full order M = 10. As in 2D, they need a higher expansion order. Si (n ≈ 4.6 + 0.7i) improves 4.2e-3 → 3.5e-5.
- **Non-physical R.** AWE gives R > R_flat or R < 0 for Ag (0.8 % of the map) and SiC (0.3 %); the hybrid never
  does.

The figures are `materials_overview.png` (where the correction is needed; true error AWE vs hybrid per material)
and `materials_gallery.png` (hybrid R/R_flat and |R_hyb − R_AWE| maps for Ag, Au, Al, TiN, Si, TiO2, SiC and
water).

## 7. Limitations

- **Grids.** The maps use 9²–21² points per window instead of 100² (about 0.05–0.4 s per corrected point). The
  32 × 32 metal scenarios take about 40 min each.
- **Metals with cos4x cos4y on 32 × 32** are limited by the grid: both methods sit about 1e-3 from the refined
  reference. A resolved 3D grid (Nx = Ny = 64) needs about 1.5 GB per window for the basis; it is shown here on the
  y-invariant silver_1d_nx64.
- **Si_crossed** (see section 4) and other cases where Padé is already below 1e-9: the least-squares floor of about
  1e-8 can be slightly worse.
- `paper_figures_3D.py`, `mms_error_3D.py` and `test_single_eps_delta_3D.py` validate the HOPS *solver* with
  manufactured solutions, which this method does not change. They are unchanged and still pass in `HOPS_Python`.
