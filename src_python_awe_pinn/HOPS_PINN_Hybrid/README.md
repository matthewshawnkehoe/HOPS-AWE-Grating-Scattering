# HOPS/AWE + least-squares PINN: a hybrid for 2D grating scattering

This directory combines the 2D HOPS/AWE method of Kehoe & Nicholls (*J. Sci. Comput.* 2024) with the
least-squares interface PINN of `../PINN_HOPS`, the PINN variant that worked best. It then compares the
hybrid with HOPS/AWE on reflectivity maps, energy defects and the inverse problem of Kaplan & Nicholls
(*Appl. Numer. Math.* 143, 2019).

The code is separate from `HOPS_Python` and `PINN_HOPS` and imports both (siblings by default; override
with `HOPS_PYTHON=...`, `PINN_HOPS=...`). The AWE numbers come from `refl_map.py` itself, so the baseline is
exactly the paper's code.

```bash
pytest -q tests                  # 5 tests, ~15 s
python compare_maps.py           # 12 map scenarios: AWE vs hybrid vs pointwise truth (~25 min)
python compare_orders.py         # accuracy vs expansion order N = M
python summarize.py              # results/summary_all.md, adaptive.md, summary.png
python inverse_paper.py          # Kaplan & Nicholls grating identification (GN / LM)
```

## The idea: physics-informed summation of the HOPS/AWE series

For each frequency band, HOPS/AWE computes once the coefficient fields of the joint expansion

  u(x, z'; ε, δ) = Σ_{n≤N} Σ_{m≤M} u_{n,m}(x, z') ε^n δ^m,  and the same for w,

on the TFE grid, with ω = ω̄(1 + δ). It then sums that series (Taylor or Padé) at every point of the map.
The summation is where AWE loses accuracy:

* the series is truncated (and `refl_map.m` sums only to half order in the Taylor case);
* the δ-series has branch points at the Rayleigh/Wood anomalies, including those of the lower layer
  inside the band;
* the ε-series has a finite radius of convergence, and Padé can place spurious poles.

**The hybrid keeps HOPS's coefficient fields and uses them as the hidden layer of a least-squares PINN**
(`hybrid/core.py`, class `PISum`):

  u_θ = Σ_j c^u_j φ^u_j,  w_θ = Σ_j c^w_j φ^w_j,  {φ_j} = {u_{n,m}} (plus, optionally, random Fourier × sin features)

At each (ε, ω) the output weights c come from one complex least-squares solve of the PINN loss of
eqs. (6a)–(6f) *at that (ε, ω)*. The loss uses:

* the Helmholtz equation in physical coordinates, via the chain rule through the ε-dependent TFE map;
* the interface conditions with the exact incident data;
* the exact FFT transparent (DtN) conditions;
* the same row weighting as the PINNs in `../PINN_HOPS`.

The Taylor weights c_{n,m} = ε^n δ^m are one admissible choice, so the hybrid's PDE residual can never be
larger than the Taylor sum's. The weights are free, so the result is not tied to the joint (ε, δ)
perturbation form, only to the span of the HOPS fields. This is the same move as in the attached
Kaplan–Nicholls paper, which replaces a fixed formula by a least-squares fit of the HOPS residual. Here the
residual is linear in the unknowns, so no Gauss–Newton iteration is needed.

Two refinements:

* **POD compression** (`compress=1e-13`) removes the numerically dependent directions of the 17 × 17 basis
  (578 → ~360 unknowns). R does not change (test H5), and the solve is 2× faster.
* **A-posteriori error indicator.** The PINN loss of the AWE Taylor sum costs about 5 ms per point and
  needs no reference solution. It ranks the true AWE error with Spearman correlation 0.94–1.00 in every
  scenario (bottom-right panel of each `results/<scenario>/maps.png`). The **adaptive hybrid** solves the
  least-squares problem only where the indicator exceeds a tolerance.

## Results: reflectivity maps and energy defects

The truth is HOPS solved separately at every grid point (δ = 0, so there is no frequency expansion), with
N = 24, Padé summation, Nx = 64 (or 128 for cos 4x) and Nz = 48. The maps use the refl_map band grids
(21 × 21, or 15 × 15 for metals). Full table: `results/summary_all.md`; figure: `results/summary.png`.

| scenario | AWE (refl_map.py) max err R | **hybrid** max err R | median err R AWE / hybrid | max err D AWE / hybrid |
|---|---|---|---|---|
| dielectric n = 1.1, band q = 1 (paper Fig. 9) | 1.8e-5 | **6.5e-15** | 2.8e-9 / 2.7e-15 | 4.3e-5 / 7.6e-14 |
| dielectric, q = 2 | 4.9e-7 | **7.7e-15** | 5.3e-10 / 4.3e-15 | 5.9e-6 / 4.6e-14 |
| dielectric, q = 3 | 2.6e-6 | **4.3e-15** | 2.4e-9 / 4.7e-16 | 5.5e-6 / 9.2e-14 |
| dielectric, TE (thesis Fig. 29) | 1.5e-4 | **7.4e-15** | 1.5e-8 / 3.0e-15 | 2.6e-4 / 1.5e-13 |
| n_w = 1.8: two lower-layer Wood anomalies inside the band | 3.0e-4 | **2.5e-13** | 1.3e-7 / 2.9e-14 | 3.4e-3 / 8.0e-14 |
| dielectric, ε up to 0.6 (3× the paper's range) | 1.5e-4 | **2.8e-13** | 2.0e-7 / 3.2e-15 | 2.3e-3 / 7.7e-12 |
| silver, cos 4x, q = 1 (paper Fig. 10a), Nx = 32 | 1.2e-5 | 6.0e-6 | 5.1e-8 / 6.7e-9 | 1.2e-5 / 6.0e-6 |
| gold, cos 4x, q = 1 (paper Fig. 10b), Nx = 32 | 1.4e-5 | 4.2e-6 | 1.1e-7 / 2.3e-8 | 1.4e-5 / 4.2e-6 |
| gold, q = 3, Nx = 32 | 7.0e-3 | 4.0e-5 | 9.7e-7 / 4.4e-7 | 7.0e-3 / 4.0e-5 |
| silver, q = 1, **Nx = 64** | 2.6e-5 | **2.9e-8** | 1.8e-10 / 2.0e-11 | 2.6e-5 / 2.9e-8 |
| gold, q = 1, **Nx = 64** | 6.2e-7 | **8.0e-11** | 3.6e-11 / 3.0e-12 | 6.2e-7 / 8.0e-11 |
| gold, q = 3, **Nx = 64** | 7.6e-3 | **3.7e-7** | 3.7e-8 / 1.4e-9 | 7.6e-3 / 3.7e-7 |

* **Dielectrics: HOPS/AWE accuracy up to round-off everywhere.** The AWE error sits at the band edges and
  along the lower-layer anomalies. The hybrid removes both: in the `high_index_q1` case the band contains
  two Wood anomalies (n_w ω = 2, 3), where the δ-series cannot converge, and the hybrid still reaches 2.5e-13
  from the same coefficient fields. The energy defect falls from up to 3e-3 (AWE) to 1e-13 to 1e-11 across
  the map (lossless cases). Adding random features barely helps (2.3e-14): the HOPS basis already spans the
  solution.
* **Metals.** With the paper's Nx = 32 the hybrid gains only 2–3×, except 175× for gold q = 3. What
  remains is the x-discretisation error of the Nx = 32 grid itself (cos 4x on a metal, which the PINN study
  had already shown), and no summation can remove that. With Nx = 64 the hybrid gains 10^3 to 2 × 10^4: gold
  q = 1 reaches 8e-11 against 6e-7 for AWE, and gold q = 3 reaches 3.7e-7 against 7.6e-3.
* **Beyond the paper's ε range** (ε ≤ 0.6), the hybrid is still at 3e-13, where AWE reaches 1.5e-4
  (D up to 2e-3).

### Lower expansion orders suffice (`compare_orders.py`, `results/orders.md`)

| scenario | N = M | AWE max err R | hybrid max err R | band solve / hybrid map |
|---|---|---|---|---|
| dielectric q = 1 | 8 | 2.3e-5 | 7.5e-13 | 0.3 s / 26 s |
| dielectric q = 1 | 10 | 2.6e-5 | 6.7e-15 | 0.6 s / 40 s |
| n_w = 1.8 (anomalies) | 10 | 6.5e-4 | 1.3e-12 | 0.6 s / 40 s |
| gold q = 1, Nx = 64 | 10 | 3.6e-5 | 4.6e-9 | 0.9 s / 42 s |
| gold q = 1, Nx = 64 | 16 | 5.2e-8 | 1.9e-10 | 2.3 s / 121 s |

AWE's own error barely improves with the order: it is limited by the band edges and anomalies. The
hybrid converges quickly in the order, and N = M = 10 already gives round-off for dielectrics.

### Cost, and the adaptive hybrid (`results/adaptive.md`)

| | per point | 21 × 21 map (dielectric) |
|---|---|---|
| HOPS/AWE (one band solve + summation) | – | 1.4 s |
| AWE error indicator (PINN loss of the AWE sum) | ~5 ms | +2.5 s |
| hybrid, N = M = 16 / N = M = 10 | 0.17 s / 0.09 s | 77 s / 40 s |
| HOPS solved separately at every point (δ = 0, same Nx = 32, N = 16) | ~0.09 s | ~40 s |

The honest comparison: **per corrected point, the hybrid costs about as much as solving HOPS again at
that frequency with the same grid.** So its advantage is not raw speed. It gets that accuracy from the
single band solve, and the indicator says *which* points need it. The adaptive hybrid solves only
where the indicator exceeds a tolerance:

| scenario | tol 1e-18: solved / max err | tol 1e-10 | tol 1e-6 |
|---|---|---|---|
| dielectric q = 1 | 62 % / 3.2e-14 | 39 % / 4.7e-10 | 22 % / 1.8e-7 |
| n_w = 1.8 (anomalies) | 77 % / 2.5e-13 | 63 % / 2.4e-10 | 54 % / 7.5e-9 |
| gold q = 1, Nx = 64 | 80 % / 8.0e-11 | 47 % / 6.8e-10 | 27 % / 1.6e-7 |

(The 1e-18 run is the one measured by `compare_maps.py`; the other tolerances are evaluated offline
from the saved maps by `summarize.py`.)

A practical recipe is to run HOPS/AWE as usual, compute the indicator map (a few seconds), and apply the
physics-informed sum only where the indicator is above the accuracy you need. This gives HOPS/AWE with a
built-in error estimate and a local fix.

## The inverse problem of the attached paper (`inverse_paper.py`, `results/inverse/inverse.md`)

This is Kaplan & Nicholls's grating identification: recover g(x_j) at 32 points from u_a(x) = u(x, a) by
Gauss–Newton or Levenberg–Marquardt, with a finite-difference Jacobian, N = 10, and the profiles (6.1),
(6.3), (6.4) with their physical parameters. The study has three parts.

**A. Does the hybrid improve the fixed-frequency forward map?** No. At one frequency the basis is only the
N + 1 fields u_n of the ε-series. Minimising the PINN residual over that small span lowers the residual
100× but does not lower the error in u_a. For (6.1) at ε = 0.2 and N = 10, the errors are 1.0e-5 for
Taylor, 4.0e-5 for the hybrid and 2.0e-3 for Padé; the other cases behave the same way. On the paper's
Nx = 32 grid, with k_w = 5.5, the residual is dominated by a discretisation floor of about 1e-13. There
the least-squares weights respond to that floor rather than to the physics, which also makes the
finite-difference Jacobian unreliable. So the hybrid is **not** used for the inversion. **Its gain comes
from the frequency direction**, where the AWE basis is rich (N × M fields) and the δ-truncation is the
dominant error.

**B. The paper's experiment**, with data from the same N = 10 model (as in the paper): in most cases GN
and LM reach absolute L∞ errors of 1e-9 to 1e-13, taking 6–19 iterations to get below 1e-7, against the
paper's 2–4. The exceptions:
* (6.3) at ε = 0.001 stops at 1–2e-7;
* (6.4) with GN stalls at 6.5e-7 at ε = 0.01 and diverges at ε = 0.05 and 0.1;
* (6.4) with LM reaches 1.5e-10 at ε = 0.1 but stalls at 6e-6 at ε = 0.05. I use the reduced formulation with a black-box forward map u_a(g) and a finite-difference
Jacobian. The paper's residual (5.1) keeps the DNO equation as a separate block, so the difference in
iteration counts probably comes from the formulation.

**C. Data from the converged model** (no inverse crime): the recovery error is set by the forward-model
error. The relative error is 4e-6 to 2e-4 for (6.1) and 7e-4 to 5e-2 for (6.3) and (6.4). GN diverges with the
Padé-summed model everywhere except at ε = 0.01 for (6.1) and (6.3), and with the Taylor model for the sandbar (6.4)
at ε ≥ 0.05. LM stays robust. An accurate forward model matters more than the choice of
optimiser.

## Take-aways

* **Yes, the PINN improves HOPS/AWE for reflectivity maps and energy defects.** It turns the fixed series
  summation into a physics-informed least-squares summation of the same HOPS/AWE fields:
  * dielectrics reach round-off (1e-15) instead of 1e-4 to 1e-7 at the band edges, even with Wood
    anomalies inside the band and at 3× the paper's ε;
  * the energy defect drops from up to 3e-3 to 1e-13 to 1e-11;
  * resolved metal maps (Nx = 64) gain 10^3 to 2 × 10^4.
* The PINN loss of the AWE sum is a cheap, reliable **error indicator** (Spearman 0.94–1.00) and makes
  **adaptive correction** possible.
* Cost: each corrected point costs about one HOPS point solve at the same resolution, so the hybrid is
  worthwhile when the indicator says a small fraction of the map needs it, or when only the AWE band data
  are available.
* **The hybrid does not help** for the fixed-frequency inverse problem (the ε-series basis is too small),
  nor for metals on the paper's under-resolved Nx = 32 grid (there the grid is the error).

## Files

| file | content |
|---|---|
| `hybrid/core.py` | `AWEBand` (a refl_map.py band + its coefficient fields as a basis on the TFE grid), `PISum` (physics-informed sum, indicator, adaptive map, POD compression, ridge option), `RandomFeatures` (LSQ-PINN enrichment) |
| `hybrid/series.py` | ε-series for one sampled profile (`ShapeSeries`); forward maps `forward_taylor/pade/hybrid` for the inverse problem |
| `compare_maps.py` | the 12 map scenarios → `results/<scenario>/{maps.png, summary.json, maps.npz, truth.npz}` |
| `compare_orders.py` | accuracy vs N = M → `results/orders.{md,png,json}` |
| `summarize.py` | `results/summary_all.md`, `adaptive.md`, `summary.png` |
| `inverse_paper.py` | Kaplan & Nicholls GN/LM → `results/inverse/` |
| `tests/test_hybrid.py` | H1 chain rule, H2 consistency with HOPS, H3 band-edge fix, H4 adaptive map, H5 compression |
