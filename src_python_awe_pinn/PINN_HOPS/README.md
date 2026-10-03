# PINN vs HOPS/AWE for 2D grating scattering

This is a set of Physics-Informed Neural Network (PINN) solvers for exactly the problem solved by the 2D
HOPS/AWE code of Kehoe & Nicholls, *J. Sci. Comput.* 100:9 (2024), eqs. (6a)–(6h). It includes the
tools to compare them with HOPS/AWE on the paper's reflectivity map R and energy defect D.

**Short answer** (section 3): no gradient-trained PINN variant does much better. That covers tanh, sin and
adaptive activations, I-PINN, DeepXDE, and least-squares or variable-projection output layers, which all
land at 1e-4 to 1e-2 in R. What does do better is a PINN built around the *linearity* of the problem: the
least-squares interface PINN (`pinn_hops/lsq_pinn.py`) solves its output layer exactly and uses features in
the TFE-flattened coordinates. It agrees with HOPS/AWE to 5e-15 (dielectric) and 4e-9 (gold) in 0.5–7 s.
HOPS/AWE is still 5–25× faster per point and far faster per reflectivity map.

The project is independent of `HOPS_Python`. HOPS is only imported, by `pinn_hops/hops_reference.py`,
to produce the reference results; set `HOPS_PYTHON=/path/to/HOPS_Python` if it is not a sibling
directory. Requirements: `numpy`, `scipy`, `matplotlib`, `torch` (CPU is enough), `pytest`, `deepxde` (only for
`deepxde_solver.py`), and HOPS_Python's own requirements.

```bash
pytest -q tests                          # 10 tests (~4 min; the training test dominates)
python compare_point.py                  # PINN vs HOPS at single (eps, omega) points, ~10 min per case
python compare_refl_map.py               # parametric PINN vs refl_map.py: R and D maps, band q = 1
python compare_refl_map.py --scenario gold
python compare_lsq.py                    # least-squares interface PINN vs HOPS (points, convergence, maps)
python study_activations.py              # tanh / sin / LAAF / I-PINN / LSQ output layer / LSGD
python deepxde_solver.py --case gold     # the same problem in DeepXDE
python inverse_demo.py                   # recover (eps, n_w) from efficiencies, HOPS vs LSQ-PINN
python summarize_variants.py             # table + figure of all variants
./run_all.sh                             # everything
```

## The problem (same as HOPS)

The unknowns are u (scattered field, above the interface z = g(x) = εf(x)) and w (transmitted field,
below it), with the Bloch phase e^{iαx} removed. They satisfy:

* **(6a, 6b) Helmholtz equations:** Δu + 2iα u_x + (γ^u)² u = 0 above, and Δw + 2iα w_x + (γ^w)² w = 0
  below.
* **(6c) Continuity:** u − w = ζ on z = g, with ζ = −e^{−iγ^u g}.
* **(6d) Flux condition:** [∂_N u − iα g_x u] − τ² [∂_N w − iα g_x w] = ψ on z = g, with N = (−g_x, 1) and
  ψ = (iγ^u + iα g_x) e^{−iγ^u g}.
* **(6e, 6f) Transparent boundary conditions:** ∂_z u = T^u[u] at z = a and ∂_z w = T^w[w] at z = −b,
  where T^u and T^w are the exact DtN Fourier multipliers ±iγ_p.
* **(6g, 6h) Periodicity:** u and w are 2π-periodic in x.

R, T and D = 1 − R − T use the formulas of paper Sect. 2.2 / `energy_defect.m` (`problem.py`).

## How the PINN is built, and why (`pinn_hops/pinn.py`)

| design choice | reason |
|---|---|
| **Two networks**, u_θ (upper) and w_θ (lower), each outputs (Re, Im) | the normal derivative jumps at z = g (τ² ≠ 1 in TM, and the incident field is not part of u, w), so one smooth network for both layers fights the interface conditions |
| **Exact periodicity**: x enters only through Fourier features cos(kx), sin(kx), k ≤ K | (6g, 6h) hold by construction; K is chosen from the profile (6 for cos x, 12 for cos 4x) |
| **Exact transparent BCs** (6e, 6f): on z = a and z = −b the network is evaluated on a uniform x-grid, then FFT → ±iγ_p → inverse FFT, inside torch | the loss contains the paper's nonlocal DtN operators themselves, not an approximate absorbing layer, so the PINN and HOPS solve the *same* truncated problem |
| **Flat-interface ansatz**: u = u_flat + ε N_u, w = w_flat + ε N_w, with u_flat = r e^{iγ^u z} and w_flat = t e^{−iγ^w z} (Fresnel) | u_flat, w_flat satisfy (6a, b, e, f) exactly and ε = 0 is reproduced exactly; the networks learn only the O(ε) correction, the analogue of HOPS expanding about the flat interface. In the dielectric test the error dropped from 29 % to 0.1 % in U, with 2× less training time than bare networks (`ansatz='plain'`) |
| **Taylor-mode (forward) derivatives** of the tanh MLP | u_xx, u_zz and the first derivatives come out of one forward pass, about 2–3× faster than nested reverse-mode autograd; checked against autograd to 1e-15 (test P3) |
| **Collocation uniform in the TFE variable z'** (z = g + s(a − g), s ~ U(0, 1)), resampled every 250 Adam steps | follows the moving interface, like the paper's change of variables |
| **Nondimensional residuals** (PDE / (1 + \|k\|²), fluxes / (1 + \|k\|)), weights λ = 10 on (6c–f) | balances terms of different physical size; the interface and boundary conditions carry the data |
| **float64; Adam (lr 2e-3 decaying to 1e-4), then L-BFGS (strong Wolfe)** on larger fixed point sets | the standard recipe to reach PINN accuracy limits; float32 stalls around 1e-4 |
| **Parametric PINN** u_θ(x, z; ε, ω) (`eps_range`, `omega_range`) | one network covers a whole band, the PINN counterpart of the HOPS/AWE joint (ε, δ) expansion; every coefficient (γ, TBC multipliers, the moving interface) is evaluated per sample |

**Verification that the loss *is* problem (6).** The HOPS solution, turned into a differentiable function
of physical (x, z) through its Fourier–Chebyshev interpolant (`HOPSFieldTorch`), has PINN loss 1e-22 to
1e-24 for the dielectric, silver, gold and oblique (α = 0.1) cases (test P2), and the exact Fresnel
solution has loss 1e-33 (test P1). So PINN and HOPS solve the same boundary-value problem, and every
remaining difference is PINN training error.

This check also works in reverse, as an independent check of HOPS resolution: for cos 4x, HOPS with
Nx = 32 has a PINN-loss of 1e-4 (aliasing), Nx = 64 gives 7e-10, and Nx = 128 gives 3e-15.

## Results

All numbers come from the scripts in this repository. They were run on 2 CPU cores; the files are in
`results/`.

### 1. Single points: `compare_point.py` (budget "standard": 3000 Adam + 3000 L-BFGS iterations)

The HOPS reference is `hops_point`: δ = 0, N = 16, Padé in ε, Nx = 128, Nz = 48. Its PINN-loss is at
round-off level, so it is a converged reference.

| case | R HOPS | R PINN | rel. err R | D HOPS | D PINN | U rel. err | field L2 err | PINN loss | HOPS loss | time PINN / HOPS |
|---|---|---|---|---|---|---|---|---|---|---|
| dielectric n = 1.1, cos x, ε = 0.1, ω = 1.5 (Fig. 9) | 0.00224486 | 0.00224753 | 1.2e-3 | −1.9e-15 | −2.5e-5 | 1.3e-3 | 1.2e-4 | 2.5e-7 | 2.2e-22 | 582 s / 1.8 s |
| dielectric, ε = 0.2 | 0.00218053 | 0.00218374 | 1.5e-3 | −6.0e-15 | 5.8e-6 | 3.9e-3 | 2.3e-4 | 8.8e-7 | 2.3e-22 | 560 s / 1.5 s |
| dielectric, oblique α = 0.1 | 0.00222656 | 0.00222699 | 1.9e-4 | 6.7e-16 | −3.5e-6 | 1.5e-3 | 1.1e-4 | 2.4e-7 | 2.2e-22 | 541 s / 1.6 s |
| silver 0.05+2.275i, cos 4x (Fig. 10a) | 0.963003 | 0.964079 | 1.1e-3 | 0.0370 (absorption) | 0.0359 | 6.4e-3 | 3.4e-3 | 6.8e-5 | 2.4e-15 | 500 s / 1.6 s |
| gold 1.48+1.883i, cos 4x (Fig. 10b) | 0.358613 | 0.359541 | 2.6e-3 | 0.641 (absorption) | 0.641 | 2.6e-3 | 1.1e-3 | 2.6e-5 | 2.4e-16 | 571 s / 1.4 s |

(HOPS time includes computing the full volume fields for the comparison; R alone takes ~0.1 s.)
The per-case figures `results/point/<case>.png` show the HOPS and PINN total fields, the pointwise
difference, and the PINN loss history with the HOPS solution's loss as a horizontal line.

An ablation run with the "quick" budget (1500 Adam + 500 L-BFGS) on the dielectric case:

| variant | rel. err R | D | U rel. err | train time |
|---|---|---|---|---|
| bare networks (no flat ansatz), reverse-mode derivatives, 2000 Adam + 500 L-BFGS | 2.5e-1 | 1.2e-2 | 2.9e-1 | 407 s |
| **flat-interface ansatz + Taylor-mode derivatives (default)** | 1.3e-3 | −1.2e-4 | 4.5e-2 | 211 s |
| + per-layer balancing of the losses by the Fresnel amplitudes (tried, removed) | 1.2e-2 | 1.5e-3 | 3.1e-1 | 166 s |

### 2. Reflectivity map and energy defect: `compare_refl_map.py` (parametric PINN, band q = 1)

One parametric PINN u_θ(x, z; ε, ω) is trained over the whole refl_map band, ω ∈ 1.5(1 ± 0.99/3) and
ε ∈ [0, 0.2]. It is compared with `refl_map.py` (HOPS/AWE, one joint expansion) on a 21 × 21 (ε, ω) grid.

| scenario | HOPS/AWE | parametric PINN |
|---|---|---|
| dielectric (Fig. 9, Taylor, N = M = 16), budget "standard" | 0.66 s for the band; median log10 \|D\| = −8.0, max −4.4 | 2034 s training + 0.7 s evaluation; median \|R_PINN − R_HOPS\| = 1.4e-5 (0.6 % of R), max 8.2e-5; median log10 \|D\| = −3.8, max −2.7 |
| gold (Fig. 10b, Padé, N = M = 15, cos 4x), budget "quick" | 0.83 s for the band | 833 s training; median \|R_PINN − R_HOPS\| = 0.012 (3.4 %), max 0.15; absorptance median agrees (10^−0.19) but the λ-dependence at large ε is wrong |

`results/refl_map/dielectric_q1.png`: R/R_flat varies by only 7 % over this band (R ≈ 0.0022). The
parametric PINN's R error (~1e-5) is comparable to that variation, so its R/R_flat map does **not**
reproduce the smooth HOPS map. Its energy defect (~1e-4) is 4–11 orders of magnitude above HOPS/AWE.
In contrast, the single-point PINNs above resolve R to 1e-3 relative: learning the whole two-parameter
family with one network costs roughly an order of magnitude in accuracy, at 3–4× the training time.

`results/refl_map/gold_q1.png`: for gold the parametric PINN (loss plateau ~5e-4, versus 2.6e-5 for the
single-point gold PINN) is right near the band centre, where a single-point PINN at ω = 1.5, ε = 0.1
agrees with HOPS to 0.3 %. Towards the band edges and at large ε it drifts by up to 0.15 in R, and at
ε = 0.2 it even gets the trend of R(λ) wrong. Across the band |k_w|² changes by a factor of 4 and the
skin layer by a factor of 2, which the network did not capture in this budget.

(The container this was run in restarts roughly every 15–25 minutes, which killed the gold run with the
"standard" budget; hence the "quick" budget here. `python compare_refl_map.py --scenario gold` runs the
standard budget, which reached a loss of ~3e-4 before being interrupted: a modest improvement.)


## 3. Better PINN variants: does any of them do better?

The first PINN reached about 1e-3 in R. To test whether that was a limit of the design, the same problem
(6a)–(6h) was solved with the variants proposed in the literature, and with one that exploits the
structure of the problem:

| variant | where | idea |
|---|---|---|
| activation `sin` | `pinn.py` (`activation='sin'`) | SIREN-like, oscillatory basis for a wave problem |
| adaptive activation (LAAF) | `pinn.py` (`activation='laaf'`) | layer-wise locally adaptive tanh(n a h), trainable a, n = 10 (Jagtap, Kawaguchi & Karniadakis 2020) |
| **I-PINN** | `pinn.py` (`activation=('tanh', 'sin')`) | Sarma et al. (CMAME 2024): one network per subdomain with a *different activation function* in each; the interface conditions couple them. Our two-network design already was XPINN/I-PINN-like; this adds the per-subdomain activation |
| PINN + least-squares output layer | `pinn_hops/hybrid.py` `lsq_refine` | after Adam/L-BFGS, solve the (linear) last layer exactly |
| LSGD / variable projection | `pinn_hops/hybrid.py` `lsgd_train` | alternate an exact LSQ solve of the output layer with Adam on the hidden layers (Cyr et al. 2020) |
| **DeepXDE** | `deepxde_solver.py` | the problem built from DeepXDE components: `PFNN` (one sub-network per output), `apply_feature_transform` (periodic Fourier features), `apply_output_transform` (flat ansatz), masked PDE residuals, the interface and FFT-DtN conditions as `PointSetOperatorBC`, `PDEPointResampler`, Adam + L-BFGS, float64; also DeepXDE's `"LAAF-10 tanh"` |
| **least-squares interface PINN** | `pinn_hops/lsq_pinn.py`, `compare_lsq.py` | see below |

### The least-squares interface PINN

Problem (6) is **linear** in (u, w). If the last layer of each subdomain network is linear,
u_θ = u_flat + Σ_j c^u_j φ^u_j(x, z), then every residual of (6a)–(6f) is affine in the output weights c.
The PINN loss (with exactly the weights of `pinn.py`) is then a linear least-squares problem, whose
*global* minimiser is one complex `lstsq` solve. No Adam, no L-BFGS, and no non-convex landscape.
This is the physics-informed extreme learning machine / random feature method idea (Dwivedi & Srinivasan 2020;
Chen, Chi, E & Yang 2022). Here the hidden layer is random and fixed, and everything else is kept from the
first PINN: two subdomain networks (I-PINN style, activation per subdomain optional), exact periodicity, the
exact FFT-DtN transparent conditions, and the flat-interface ansatz.

Two design choices decide the accuracy (both measured in `compare_lsq.py`):

1. **Separable Fourier features.** φ_{k,j}(x, z) = e^{ikx} σ(a_j z + b_j) for |k| ≤ K, with random a_j, b_j
   (a Fourier-feature / separable PINN). The generic one-hidden-layer network in (x, z)
   (`basis='rfm'`: σ(Σ_k A_jk cos kx + B_jk sin kx + a_j z + b_j)) drives the loss at the collocation points
   down to 1e-10, but its error stays at 1e-2: it fits the rows and oscillates between them
   (purple curves in `results/lsq/conv.png`).
2. **Features in the flattened coordinate of each layer (`coords='tfe'`).** ζ = (z − g)/(a − g) above and
   ζ = (z + b)/(g + b) below, the change of variables of the paper's TFE method, with the chain rule for the
   derivatives. With features in physical (x, z), the network must continue u *below* the crests of the
   interface (and w above its troughs). That is the Rayleigh-hypothesis difficulty, and it grows with the slope of
   g. For gold with cos 4x, the TFE coordinates are 10–100× more accurate at the same size
   (`results/lsq/coords.png`).

With these, the error falls exponentially with the network size (`results/lsq/conv.png`). The dielectric
case reaches round-off with 504 complex output weights (sin activation, 21 Fourier modes × 12 z-features per
layer); gold reaches 1e-9. `sin` converges fastest for the dielectric case; tanh, sin and Gaussian are
comparable for gold. The I-PINN pairing (tanh above, sin below) sits between its two parents.

### Results (quick budget for the gradient-trained variants: 1500 Adam + 500 L-BFGS; DeepXDE 1000 + 300)

All variants use the same equations, the same R/T/D formulas, and the same HOPS reference (N = 16, Padé,
Nx = 128) as sections 1–2. Source: `results/variants/summary.md` and `summary.png`.

| variant | dielectric: rel. err R | abs. err D | gold: rel. err R | abs. err D |
|---|---|---|---|---|
| PINN tanh (baseline, section 1 design) | 1.3e-3 | 1.2e-4 | 2.1e-2 | 7.6e-3 |
| PINN sin | 3.1e-3 | 1.3e-3 | 3.4e-3 | 1.2e-3 |
| PINN adaptive tanh (LAAF) | 5.7e-3 | 4.5e-4 | 1.8e-2 | 6.5e-3 |
| I-PINN (tanh above, sin below) | 6.7e-3 | 1.2e-3 | 2.2e-2 | 7.8e-3 |
| PINN tanh + exact LSQ output layer | 4.1e-4 | 1.2e-4 | 2.3e-2 | 8.2e-3 |
| PINN, LSGD (1500 Adam steps with LSQ output layer) | 2.1e-3 | 1.0e-3 | 6.0e-2 | 2.2e-2 |
| DeepXDE, PFNN tanh | 2.9e-4 | 1.1e-4 | 8.6e-2 | 3.1e-2 |
| DeepXDE, PFNN "LAAF-10 tanh" | 3.4e-3 | 1.3e-4 | – | – |
| **least-squares interface PINN** | **4.6e-15** | **3.3e-16** | **4.4e-9** | **1.6e-9** |

Least-squares PINN on all five single-point cases (`results/lsq/points.json`):

| case | rel. err R | abs. err D | max rel. err U (interface) | PINN loss | unknowns | time |
|---|---|---|---|---|---|---|
| dielectric, ε = 0.1 | 4.6e-15 | 3.3e-16 | 2.4e-11 | 6.6e-23 | 504 | 0.5 s |
| dielectric, ε = 0.2 | 1.0e-11 | 3.7e-13 | 2.1e-8 | 3.6e-17 | 504 | 0.6 s |
| dielectric, α = 0.1 | 2.1e-14 | 2.4e-15 | 2.5e-11 | 7.2e-23 | 504 | 0.5 s |
| silver, cos 4x | 8.0e-8 | 7.7e-8 | 1.8e-6 | 4.1e-8 | 1560 | 7 s |
| gold, cos 4x | 4.4e-9 | 1.6e-9 | 1.2e-6 | 1.2e-8 | 1560 | 7 s |

Cost per step on one thread (`results/variants/timing.json`):

| | dielectric | gold |
|---|---|---|
| HOPS point solve (δ = 0, N = 16) | 0.09 s | 0.29 s |
| least-squares PINN, one solve | 0.52 s | 7.2 s |
| gradient PINN, one Adam step (Taylor-mode derivatives, 2000 points) | 0.10 s | 0.13 s |
| DeepXDE, one Adam step (reverse-mode autodiff, 1000 points) | 0.96 s | 0.97 s |

The wall times in `summary.md` were measured with 3–4 jobs sharing 2 cores, so compare iteration budgets
and the single-thread table above, not those times.

### Reflectivity maps with the least-squares PINN (`compare_lsq.py --parts map`)

One least-squares solve per (ε, ω) grid point, on the refl_map band q = 1:

| scenario | median \|R_LSQ − R_HOPS/AWE\| | max | median log10 \|D\| (LSQ / HOPS/AWE) | LSQ time |
|---|---|---|---|---|
| dielectric (Fig. 9), 21 × 21 | 2.8e-9 | 1.8e-5 | −14.3 / −8.0 | 391 s (0.9 s per point) |
| gold (Fig. 10b), 11 × 11 | 1.0e-7 | 1.3e-5 | absorptance agrees, −0.19367 | 1118 s (9 s per point) |

The maximum differences sit at the **band edges**. Re-solving those points with HOPS at δ = 0 (no
frequency expansion, N = 24; `map_*_edges.json`) shows the least-squares PINN is the more accurate of the
two maps there. For the dielectric at (ε, ω) = (0.2, 1.005), HOPS point solve − LSQ = 9e-15, while
HOPS point solve − AWE map = 1.8e-5. That 1.8e-5 is the truncation error of the frequency (δ) series of
AWE, which `refl_map.m` sums to half order (see HOPS_Python's `--taylor-full-order`). For gold the
LSQ–point difference at the edges is 3e-9 to 1e-6, against 4e-8 to 1.4e-5 for the AWE map. The energy
defect of the LSQ-PINN map is 1e-14 everywhere, including near the Rayleigh anomalies where the AWE series
degrades (`results/lsq/map_dielectric_q1.png`).

### An inverse problem (`inverse_demo.py`)

The task: recover the grating amplitude ε and the substrate index n^w from 15 measured reflected
efficiencies (orders −1, 0, 1 at 5 frequencies). Each forward model is used inside
`scipy.optimize.least_squares`:

| forward model | noise | recovered ε (true 0.13) | recovered n^w (true 1.1) | forward solves | time |
|---|---|---|---|---|---|
| HOPS | 0 | 0.13 ± 2e-14 | 1.1 ± 1e-15 | 32 | 20 s |
| least-squares PINN | 0 | 0.13 ± 9e-14 | 1.1 ± 2e-15 | 28 | 128 s |
| HOPS | 1 % | 0.12991 | 1.09992 | 36 | 22 s |
| least-squares PINN | 1 % | 0.12991 | 1.09992 | 33 | 143 s |

Both give the same answer to 9 digits, so the error comes only from the noise. The classical "PINN for
inverse problems" (ε as a trainable parameter inside a gradient-trained PINN) would cost at least one
10-minute training and would be limited to the 1e-3 forward accuracy above.

### Neural operators (FNO, DeepONet)

These were not implemented. They learn the *map* from parameters (ε, ω, n^w, profile) to solutions from
training data, which here would come from HOPS itself. The result is a fast surrogate of HOPS whose accuracy
is bounded by the training data and typically 1e-2 to 1e-3. That is useful when millions of evaluations are
needed, but it is not a competitor for accuracy. DeepXDE's `DeepONet` could be trained on
`hops_refl_map` output if such a surrogate is wanted.

## Take-aways

* **Changing the activation function, using I-PINN, using DeepXDE, or changing the optimiser does not
  change the picture.** Every gradient-trained variant lands at 1e-4 to 1e-2 in R (tanh, sin, LAAF, I-PINN,
  DeepXDE's PFNN/LAAF, LSQ output layer, LSGD). Which one is best depends on the case: DeepXDE tanh and the
  LSQ refinement for the dielectric, sin for gold. None of them gains two orders of magnitude. The limit is
  the non-convex optimisation, not the network architecture. DeepXDE also costs 7–9× more per step here,
  because it uses reverse-mode autodiff for the Laplacian, while `pinn.py` uses Taylor-mode derivatives.
* **The variant that does better exploits the linearity of Maxwell/Helmholtz.** The least-squares interface
  PINN replaces training by one linear least-squares solve. It matches HOPS/AWE to round-off for the paper's
  dielectric case (R to 5e-15, D to 3e-16) and to 1e-9 for gold, in 0.5–7 s instead of 10 minutes. That is about
  10^11× (dielectric) and 10^6× (gold) more accurate than the best gradient-trained variant. On the reflectivity map it is more
  accurate than the AWE frequency expansion at the band edges.
* **HOPS/AWE is still faster.** One HOPS point solve costs 0.1–0.3 s, 5–25× less than the least-squares
  PINN. A whole band map costs HOPS/AWE 0.8 s against 391–1118 s for the least-squares PINN (one solve per
  grid point). The two use the same ingredients: Fourier modes in x, the TFE flattening, and exact DtN
  conditions. HOPS adds the perturbation series in (ε, δ), which gives the whole map almost for free.
  The least-squares PINN needs no series: no radius of convergence and no degradation at anomalies. That
  makes it a good independent check of HOPS/AWE, which is how it found the band-edge truncation error.
* **For metals** the gradient PINNs have the most trouble (skin layer; 3e-3 to 9e-2 in R), and so does the
  physical-coordinate least-squares PINN at large slopes (the Rayleigh effect). The TFE-coordinate version
  fixes the latter. At ε = 0.2 with cos 4x it still needs K ≈ 32–40 Fourier modes, as HOPS needs Nx = 128.

### First-study take-aways (gradient-trained PINN, sections 1–2)

* **Accuracy.** Once trained, the PINN reproduces R and the fields to about 1e-3 relative, and D to
  about 1e-5 (dielectric). HOPS/AWE reaches round-off (D ~ 1e-15) at the same point.
* **Cost.** One PINN solve takes 3–10 minutes of training on 2 CPU cores; HOPS takes under 1 s. For a
  whole reflectivity map, HOPS/AWE needs one solve per band plus cheap summation. The parametric PINN
  also needs just one (longer) training per band, but it is 3–4 orders of magnitude less accurate.
* **Where the PINN is structurally different.** It needs no series in ε or δ, so it has no radius of
  convergence. Its accuracy does not collapse at the lower-layer Rayleigh anomalies inside a band,
  where the HOPS/AWE Taylor series has a branch point (see the D maps). Its error is instead set by
  optimisation, spread fairly uniformly.
* **Metals** (thin skin-depth boundary layer, |k_w| ≈ 3.4 at ω = 1.5, K = 12 modes for cos 4x) are the
  hardest case for the PINN. HOPS handles them with no extra effort.

## References for the variants

* A. K. Sarma, S. Roy, C. Annavarapu, P. Roy, S. Jagannathan, *Interface PINNs (I-PINNs): A physics-informed
  neural networks framework for interface problems*, CMAME 429 (2024) 117135.
* A. D. Jagtap, K. Kawaguchi, G. E. Karniadakis, *Locally adaptive activation functions with slope recovery
  for deep and physics-informed neural networks*, Proc. R. Soc. A 476 (2020).
* L. Lu, X. Meng, Z. Mao, G. E. Karniadakis, *DeepXDE: a deep learning library for solving differential
  equations*, SIAM Review 63 (2021).
* J. Chen, X. Chi, W. E, Z. Yang, *Bridging traditional and machine learning-based algorithms for solving
  PDEs: the random feature method*, J. Mach. Learn. 1 (2022).
* V. Dwivedi, B. Srinivasan, *Physics informed extreme learning machine (PIELM)*, Neurocomputing 391 (2020).
* E. C. Cyr, M. A. Gulian, R. G. Patel, M. Perego, N. A. Trask, *Robust training and initialization of deep
  neural networks: an adaptive basis viewpoint* (LSGD), MSML 2020.
