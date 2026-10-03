# HOPS/AWE + PINN: every 2D scenario, with an error analysis

This directory re-runs **every scenario of the 2D HOPS/AWE code** (`Hops_2D_AWE`: `refl_map.py`, `paper_figures.py`,
`material_survey.py`, `refl_movie.py`) with the best combined HOPS/AWE + PINN method from `HOPS_PINN_Hybrid`. It
compares the two methods against converged references and uses a 64-material survey to see when the standard AWE
summation fails.

## 1. The method: physics-informed summation (PI-sum)

HOPS/AWE expands the fields in surface height (ε^n) and frequency (δ^m). It then sums u = Σ u_{n,m} ε^n δ^m with
Taylor or Padé. The PI-sum treats the coefficient fields u_{n,m} and w_{n,m} as the **hidden layer of an interface
PINN**. The output weights come from one linear least-squares solve of the PINN loss at each (ε, ω). The loss
includes:

- the Helmholtz equation in each layer, written in physical coordinates through the TFE map;
- the two interface conditions;
- the exact FFT Dirichlet–Neumann (DtN) boundary conditions at the top and bottom.

This was the clear winner of the earlier hybrid study, ahead of a PINN on its own, PINN correction of the Padé sum and
neural operators. It needs no training, is deterministic, and reduces to the AWE Taylor sum when that sum is already
exact.

The approach was generalised here to all `refl_map.py` options: multiple frequency windows, `max_delta` and `joint`
windows for dispersive materials, TE/TM, and any profile. It is also made **adaptive**:

1. `refl_map.run(..., keep_fields=True)` computes the standard AWE result, unchanged.
2. **Indicator:** the PINN loss of the full-order AWE Taylor sum, at every map point (about 5 ms each).
3. If the indicator is above `tol = 1e-14`, the weights come from least squares over the POD-compressed basis with
   n, m ≤ 12. Otherwise the Taylor sum is kept.
4. R, T and the energy defect D = 1 − R − T are replaced in the window result. The AWE values are kept as
   `ru_awe`, `ee_awe`, … so `refl_map.plot` and `refl_movie` work on the result unchanged.

Code: `awepinn/core.py` (`run(scenario, ...)` is a drop-in for `refl_map.run`).

## 2. How to run

The code requires `HOPS_Python` and `HOPS_PINN_Hybrid` as sibling directories. Alternatively, set the environment
variables `HOPS_PYTHON` and `HOPS_PINN_HYBRID`.

```bash
python run_scenarios.py                       # all refl_map scenarios -> figures/<scenario>/, results/scenarios/
python run_scenarios.py --only gold thesis22  # a subset;  --procs N for parallel scenarios;  --replot
python error_analysis.py                      # rigorous comparison vs converged HOPS -> results/error_analysis.md
python hybrid_movie.py                        # movies -> figures/movies/*.mp4 (needs ffmpeg)
python materials_hybrid.py --procs 2          # 64-material survey + learning -> results/materials.md
python materials_hybrid.py --gallery          # gallery of 8 materials
python -m pytest tests -q                     # 4 tests
```

```python
from awepinn import core
res, info = core.run('gold', N_Eps=31, N_delta=31)   # same output as refl_map.run, with hybrid R, T, D
```

## 3. Scenarios reproduced

The tables below cover 41 of the 43 `refl_map.py` scenarios, plus silver and gold at Nx = 64 (`silver_nx64`, `gold_nx64`). Full table:
`results/scenarios.md`.

For every scenario there are two kinds of output:

- **Figures** in `figures/<scenario>/`:
  - `*_awe.png` and `*_hybrid.png`: the reflectivity map and energy-defect map in the `refl_map.plot` format;
  - `compare.png`: R and |D| spectra and maps, side by side.
- **Data** in `results/scenarios/<scenario>.npz` and `.json`.

**Energy defect (max log10|D| for lossless structures)**

| scenario(s) | AWE | hybrid |
|---|---|---|
| dielectric, thesis17, thesis18 | −4.2 | −10.5 / −10.1 |
| thesis29 / 30 (TE) | −3.6 | −9.8 / −9.2 |
| thesis21 | −1.0 | −3.6 |
| thesis22 | −0.0 | −7.3 |
| thesis23 | −2.0 | −4.4 |
| TiO2 | −2.5 | −8.4 |
| high_contrast | −2.5 | −7.9 |

**Other notable scenarios**

- **glass_over_air:** AWE gives a non-physical R/R_flat of up to 2.6. The hybrid result is physical.
- **water_over_gold:** the AWE Padé sum has a spike with ΔR = 1.1, which the hybrid removes.
- **thesis26 / 27** (Nx = 512 and 256): these were run on a 9×9 grid with basis order 8 to fit in memory. At that
  order the hybrid is slightly *worse* (|D| 10^-1.5 → 10^-0.9). The full order would be needed; see the limitations
  section.
- **thesis28a / c** (Nx = 1024): not run. The basis fields do not fit in memory in this environment.

## 4. Error analysis (`results/error_analysis.md`, `figures/error_analysis/`)

The comparison uses **4,970 samples** from 41 scenarios (39 + the Nx = 64 reruns of silver and gold). Each sample is checked against two references:

- **Total error:** a converged HOPS reference with no frequency expansion, N = 24 with Padé, 2·Nx, and Nz + 16.
- **Summation error:** HOPS on the scenario's *own* grid. This isolates the error made by the summation step, which
  is the only step the PINN changes.

**Pooled results**

| | AWE (refl_map) | HOPS/AWE + PINN |
|---|---|---|
| max error in R | 0.58 | 0.095 |
| median error in R | 2.1e-10 | **5.4e-14** |
| max summation error | 0.58 | 0.11 |
| median summation error | 6.2e-11 | 4.8e-14 |
| better by > 2× | – | 53 % of samples (worse at 13 %, almost all metals on Nx = 32) |

The hybrid is not worse than the AWE at the remaining samples; there both methods agree to within a factor of 2.

**By scenario family**

| family | max error in R: AWE | max error in R: hybrid | notes |
|---|---|---|---|
| dielectric / thesis17 / thesis18 | 6.3e-5 | 2.2e-11 | hybrid better at 96 % of samples; Wilcoxon p = 2e-24 |
| TE dielectric (thesis29 / 30) | 1.8e-4 | 1.0e-10 | |
| thesis21 / 22 / 23 (high index) | 1.9e-2 / 2.6e-2 / 5.3e-3 | 1.1e-4 / 4.1e-7 / 9.4e-7 | |
| named dispersive materials (Ag, Al, Na, TiN, ITO, AZO, VO2, Si, TiO2, SiC, sapphire) | 4e-6 … 3.5e-3 | 9e-12 … 8e-9 | |
| glass_over_air | 0.58 | 1.4e-9 | |
| PEC-like (thesis25, 32, 33a/c) | 3e-3 … 3.5e-2 | 3e-7 … 7e-5 | |
| thesis24 (n = 20i, TM) | 7.6e-2 | 9.5e-2 | needs basis order 15; see below |
| metals with cos4x profile, Nx = 32 (silver, gold, thesis19/20/31) | 2e-3 … 0.13 | 3e-4 … 1.3e-2 | **limited by the discretization**, not the summation |

The **full-order AWE Taylor sum** (the "AWE full Taylor" column of the report) diverges, up to 1e+69. The
`refl_map.py` result is usable only because of Padé. The PINN weights do not need Padé.

Further findings:

- **Metals on Nx = 32.** The cos4x profile under-resolves the surface plasmon. Even the same-grid HOPS reference
  differs from the converged one by about 1e-2. The hybrid reduces the maximum error by about 10×, but its median is
  no better (Wilcoxon p ≈ 0.8). This case is limited by the grid, not by the summation. The Nx = 64 reruns of silver
  and gold test this directly (see below).
- **Indicator calibration.** Spearman ρ between the indicator and the actual error is 0.92, with
  log10|err| ≈ 0.71·log10(indicator) − 0.58. The indicator can therefore be used as an a-posteriori error estimate
  for the AWE map at almost no cost.
- **Tolerance sweep.** Setting tol between 1e-20 and 1e-6 changes the fraction of points that are solved from 64 %
  to 34 %. The maximum error does not change. Points the indicator accepts are already accurate.
- **Basis order.** Errors fall steadily from order 6 to order 12 (for example, dielectric 1.2e-7 → 5e-13). PEC-like
  substrates need high order: thesis24 gives 0.61 (order 6), 0.095 (order 12) and 9e-5 (order 15, the full basis).
  When |n_w| ≥ 20, use `basis_order = N`.

**The resolution test confirms the metal diagnosis.** Silver and gold were rerun with Nx = 64; everything else was
unchanged.

| scenario | AWE max error | hybrid max error | hybrid better at | Wilcoxon p |
|---|---|---|---|---|
| silver, Nx = 32 | 0.11 | 1.3e-2 | 41 % | 0.8 |
| **silver, Nx = 64** | 0.19 | **7.1e-5** | 71 % | 8e-15 |
| gold, Nx = 32 | 2.2e-2 | 9.0e-4 | 45 % | 7e-5 |
| **gold, Nx = 64** | 2.7e-2 | **3.6e-7** | 71 % | 3e-15 |

On the finer grid, the hybrid error falls by 180× (silver) and 2,500× (gold). The AWE error does not fall at all,
because the AWE error comes from the summation, and the PINN step removes it. On a resolved grid, metals behave like
the dielectrics.

## 5. Movies (`figures/movies/`)

| file | content |
|---|---|
| `hybrid_dielectric.mp4`, `hybrid_silver_spp.mp4`, `hybrid_silver_spp_eps.mp4`, `hybrid_gold_spr.mp4`, `hybrid_sic_sphp.mp4`, `hybrid_tir.mp4` | the `refl_movie.py` presets, with fields and R/D from the hybrid method |
| `compare_wood_anomalies.mp4` | **new.** Sweep through the Wood anomalies (n_w = 1.8). Shows the AWE vs hybrid spectra, \|D\|, the indicator, the hybrid scattered field and log10\|u_hyb − u_AWE\| |
| `compare_paper_fig9_edges.mp4` | **new.** Same layout, at the band edges of paper Fig. 9 (ε = 0.2), where the AWE degrades |
| `compare_large_eps.mp4` | **new.** Sweep in ε for thesis23 (n = 8.1). The AWE breaks down as ε grows; the hybrid stays physical |

## 6. What the materials survey shows (`results/materials.md`, `figures/materials/`)

The survey covers 64 materials from `hops/materials.py`, each on its own band and period, with a cos x grating and a
full map. At every point the hybrid and the AWE were compared, along with the indicator.

- **Absorption (Im n) is the strongest predictor of AWE trouble.** Spearman ρ = 0.84 with the max
  |R_hyb − R_AWE|; ρ = 0.70 with the fraction of the map that is off. Re ε (ρ = −0.71) and being lossless
  (ρ = −0.70) predict the opposite. The SPP condition |Re ε + 1| is a weak predictor (ρ = −0.34).
- **A two-level decision tree** classifies 92 % of the materials correctly for "AWE is off by more than 0.01
  somewhere":
  - if Re ε ≤ −3.1 (a good metal), AWE is off;
  - otherwise, if Re n > 2.7 (high index), AWE is off;
  - otherwise AWE is acceptable.

  This is a practical rule for when to switch on the PINN step.
- **Worst materials:**
  - Mo: ΔR = 1.3, and the AWE gives 0.7 % non-physical R > R_flat points;
  - Ag_IR and Au_IR: 0.21 (n ≈ 20i, the PEC-like limit);
  - Pt and Ag: 0.18.

  Air is at 2e-11. Lossless dielectrics have an AWE energy defect of 10^-1.3 to 10^-2.2, which the hybrid reduces to
  10^-7.5 to 10^-8.
- **The hybrid never produced a non-physical R, in any material.** The AWE produced them for Mo, Ag, W, Ag_IR and
  Au_IR.
- The earlier cos4x / Nx = 32 survey is kept in `results/materials_cos4x_nx32*`. It shows that when the grid is too
  coarse, both methods fail together, and the indicator flags about 89 % of the map. That high flag rate is itself a
  useful warning to refine the grid.

## 7. Limitations

- **thesis28a / c** (Nx = 1024) were not run. They are too large for the memory available here, though they fit on a
  workstation with more than 32 GB of RAM.
- **thesis26 / 27** used basis order 8 on a 9×9 grid, which makes the hybrid slightly worse than AWE there. Rerun them
  with the full basis order: `python run_scenarios.py --only thesis26 thesis27`, after removing them from the
  `COARSE` tuple in `run_scenarios.py`.
- The hybrid cannot fix **discretization error**. That comes from a too-small Nx for metals with a cos4x profile.
- `mms_error.py` and `test_single_eps_delta.py` test the HOPS *solver*, which this method does not change. Those
  results are unchanged and still pass in `HOPS_Python`.
- **Cost.** The least-squares step costs about 0.1–0.5 s per point. A full 6-window 31² map takes about 7 min,
  compared with about 7 s for the AWE alone. The adaptive indicator avoids the least-squares step on 45 % of the
  error-analysis samples, and on more than 90 % of the points for the named dispersive-material maps.
