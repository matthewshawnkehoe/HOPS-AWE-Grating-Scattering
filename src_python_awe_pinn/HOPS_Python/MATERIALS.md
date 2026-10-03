# Materials for the HOPS/AWE reflectivity map (refractiveindex.info)

## 1. What the database contains
The refractiveindex.info database (M. N. Polyanskiy, CC0 public domain; GitHub
`polyanskiy/refractiveindex.info-database`, release 2026-05-24) is a library of YAML files. It is organized in
**shelves → books (materials) → pages (one data set from one paper)**:

| shelf | content | books | pages |
|---|---|---|---|
| main | simple inorganic materials (elements, oxides, nitrides, semiconductors, halides) | 256 | 1158 |
| organic | liquids, solvents, polymers, dyes | 91 | 216 |
| glass | fused silica, BK7, soda-lime, ZBLAN, … | 7 | 18 |
| other | alloys, doped/mixed crystals (ITO, AZO, AlGaAs), perovskites, liquid crystals, resists, biological tissue | 112 | 313 |
| specs | manufacturer sheets (SCHOTT, OHARA, HOYA, CDGM, …) | 18 | 1802 |
| 3d, popular_glass | curated subsets | 14 | 69 |

Each page is stored as **tabulated n,k** (919 files), **tabulated n** (254), **tabulated k** (1739, usually next to a
formula for n), or one of **nine dispersion formulas** (Sellmeier 1/2, polynomial, "RefractiveIndex.INFO",
Cauchy, gases, Herzberger, Retro, Exotic; about 2,300 uses). `hops/materials.py` reads all of these.

## 2. Which materials give interesting reflectivity maps
The HOPS/AWE map depends on the material only through ε = n² of each layer. The database spans these regimes:

| regime | ε | what the map shows | good choices (curated key) |
|---|---|---|---|
| low-loss plasmonic metals | Re ε ≪ −1, small Im ε | sharp **surface-plasmon (SPP) lines** next to every Rayleigh anomaly; they deepen with ε and bend with it | `Ag`, `Na`, `K`, `Al` (UV), `Mg` (UV), `Au`/`Cu` (red–NIR), `AuAg50` |
| lossy (transition) metals | n ≈ k | smooth, broadband absorption; weak grating features (thesis Figs. 20 and 31) | `W`, `Fe`, `Co`, `Cr`, `Ni`, `Ti`, `Pt`, `Pd`, `Mo`, `steel`, `brass` |
| alternative plasmonics / ENZ | Re ε crosses 0 | the character changes inside the window: dielectric → epsilon-near-zero → metal | `TiN` (VIS), `ITO` (ENZ ≈ 1.3–1.6 µm), `AZO` (mid-IR) |
| phase-change | switches with temperature | the same grating, two very different maps | `VO2_cold` / `VO2_hot` |
| high-index lossless dielectrics | ε ≈ 4–16 | many **lower-layer Wood anomalies** inside each band (use `--windows joint`) | `Si_IR`, `Ge_IR`, `GaAs` (NIR), `TiO2`, `GaP`, `diamond`, `Si3N4`, `LiNbO3`, `GaN`, `ZnGeP2` |
| absorbing semiconductors (above the gap) | large Re ε, moderate Im ε | high R with smooth structure; absorption changes across the gap | `Si` (VIS/UV), `Ge`, `GaAs`, `InP`, `CdTe`, `MAPbI3`, `MoS2` |
| low-index dielectrics | ε ≈ 1.8–2.4 | weak R (2–10 %); Wood anomalies dominate | `SiO2`, `sapphire`, `HfO2`, `ZnO`, `Ta2O5`, `AlN`, `ZnS` |
| polar crystals (Reststrahlen band) | Re ε < 0 in the mid-IR | **surface phonon polaritons**: the mid-IR analogue of silver | `SiC` (Lorentz model; the database only has an amorphous film), `sapphire_IR`, `SiO2_IR` |
| superstrates n^u | 1.0–2.4 | sensing (SPR shift), total internal reflection with α ≠ 0 | `water`, `ethanol`, `PDMS`, `PMMA`, `fused_silica`, `BK7`, `CaF2`, `MgF2`, `ZnSe` |
| non-physical (thesis) | purely imaginary n, n^u > 1 | the thesis' "interesting patterns" | presets `thesis23`–`thesis25`, `thesis32`, `thesis33a/c` |

## 3. Survey of the curated library
From `python material_survey.py`: vacuum over the material, TM, f = cos(4x), ε ≤ 0.2, Padé, N = M = 12, 40×40 grid,
joint windows, each material at its suggested period P. Columns:
- **R range**: absolute reflectivity.
- **min R/R_flat**: depth of the grating resonances.
- **D**: energy defect, for lossless materials only; it measures accuracy.
- **Absorptance**: = D for opaque substrates.

The gallery is in `figures/material_survey/material_survey.png`.

| material | class | P (µm) | n (band 1) | R range (1–99 %) | min R/R_flat | D (lossless, ε≤0.05) | median absorptance |
|---|---|---|---|---|---|---|---|
| Ag | plasmonic metal | 1 | 0.0484+4.54j | 0.02–0.99 | 0.08 | – | 0.75 |
| Au | plasmonic metal | 1 | 0.138+3.76j | 0.01–0.96 | 0.03 | – | 0.77 |
| Cu | plasmonic metal | 1 | 0.218+3.82j | 0.01–0.95 | 0.04 | – | 0.68 |
| Al | plasmonic metal | 0.5 | 0.33+4j | 0.20–0.93 | 0.60 | – | 0.13 |
| Na | plasmonic metal | 1 | 0.0486+2.79j | 0.20–0.98 | 0.21 | – | 0.10 |
| K | plasmonic metal | 1.5 | 0.04+3.01j | 0.31–0.98 | 0.49 | – | 0.30 |
| Mg | plasmonic metal | 0.5 | 0.22+2.86j | 0.15–0.91 | 0.17 | – | 0.16 |
| AuAg50 | plasmonic metal | 1 | 0.308+4.18j | 0.01–0.94 | 0.04 | – | 0.70 |
| Au_IR | plasmonic metal | 10 | 5.19+45.3j | 0.94–0.99 | 0.95 | – | 0.01 |
| Ag_IR | plasmonic metal | 10 | 4.58+46.8j | 0.95–0.99 | 0.95 | – | 0.01 |
| TiN | alternative plasmonic | 1 | 1.43+2.94j | 0.01–0.61 | 0.05 | – | 0.73 |
| ITO | ENZ / TCO | 2 | 0.278+0.737j | 0.02–0.49 | 0.20 | – | 0.15 |
| AZO | ENZ / TCO | 10 | 2.24+5.24j | 0.02–0.76 | 0.04 | – | 0.60 |
| W | lossy metal | 1 | 0.914+7.13j | 0.14–0.93 | 0.18 | – | 0.32 |
| Fe | lossy metal | 1 | 2.91+3.12j | 0.01–0.53 | 0.03 | – | 0.73 |
| Co | lossy metal | 1 | 2.26+4.3j | 0.01–0.69 | 0.02 | – | 0.65 |
| Cr | lossy metal | 1 | 3.08+3.35j | 0.01–0.56 | 0.03 | – | 0.65 |
| Ni | lossy metal | 1 | 2+4.3j | 0.01–0.71 | 0.03 | – | 0.65 |
| Ti | lossy metal | 1 | 2.78+3.86j | 0.01–0.62 | 0.03 | – | 0.64 |
| Pt | lossy metal | 1 | 0.482+6.53j | 0.03–0.96 | 0.05 | – | 0.62 |
| Pd | lossy metal | 1 | 1.81+4.46j | 0.01–0.74 | 0.03 | – | 0.61 |
| Mo | lossy metal | 1 | 0.954+9.33j | 0.21–0.96 | 0.25 | – | 0.33 |
| brass | lossy metal | 1 | 0.444+3.83j | 0.01–0.89 | 0.04 | – | 0.73 |
| steel | lossy metal | 1 | 2.51+4.26j | 0.01–0.67 | 0.02 | – | 0.65 |
| VO2_cold | phase change | 5 | 3.16+0.221j | 0.03–0.27 | 0.12 | – | 0.74 |
| VO2_hot | phase change | 5 | 4.35+5.19j | 0.03–0.69 | 0.05 | – | 0.67 |
| Si | semiconductor | 1 | 3.81+0.0131j | 0.15–0.70 | 0.21 | – | 0.36 |
| Si_IR | semiconductor | 5 | 3.43+2.34e-11j | 0.05–0.33 | 0.16 | – | 0.00 |
| Ge | semiconductor | 1 | 5.19+0.528j | 0.11–0.67 | 0.16 | – | 0.44 |
| Ge_IR | semiconductor | 10 | 3.96+0j | 0.07–0.38 | 0.19 | – | 0.00 |
| GaAs | semiconductor | 1 | 3.68+0.152j | 0.03–0.55 | 0.07 | – | 0.52 |
| InP | semiconductor | 0.5 | 3.1+1.81j | 0.01–0.46 | 0.02 | – | 0.61 |
| GaP | semiconductor | 0.5 | 4.96+2.47j | 0.03–0.58 | 0.05 | – | 0.48 |
| InSb | semiconductor | 5 | 3.9+0.0985j | 0.05–0.43 | 0.15 | – | 0.46 |
| CdTe | semiconductor | 1 | 2.95+0.252j | 0.01–0.37 | 0.04 | – | 0.66 |
| MAPbI3 | semiconductor | 1 | 2.48+0.255j | 0.02–0.20 | 0.13 | – | 0.86 |
| MoS2 | semiconductor | 1 | 5.02+1.46j | 0.03–0.55 | 0.05 | – | 0.53 |
| SiO2 | dielectric | 1 | 1.46+0j | 0.01–0.04 | 0.30 | 1e-12 | – |
| sapphire | dielectric | 1 | 1.76+0j | 0.01–0.10 | 0.16 | 1e-11 | – |
| Si3N4 | dielectric | 1 | 2.04+0j | 0.02–0.14 | 0.18 | 1e-11 | – |
| TiO2 | dielectric | 1 | 2.57+0j | 0.03–0.23 | 0.13 | 1e-10 | – |
| ZnO | dielectric | 1 | 1.98+0j | 0.02–0.13 | 0.14 | 1e-11 | – |
| ZnGeP2 | dielectric | 5 | 3.13+0j | 0.04–0.29 | 0.14 | 1e-10 | – |
| diamond | dielectric | 1 | 2.41+0j | 0.01–0.30 | 0.07 | – | 0.00 |
| HfO2 | dielectric | 0.5 | 1.98+0j | 0.02–0.17 | 0.14 | 1e-11 | – |
| Ta2O5 | dielectric | 1 | 2.13+0j | 0.02–0.16 | 0.11 | – | 0.00 |
| LiNbO3 | dielectric | 1 | 2.28+0j | 0.03–0.18 | 0.18 | 1e-10 | – |
| GaN | dielectric | 1 | 2.38+0j | 0.03–0.21 | 0.16 | 1e-10 | – |
| AlN | dielectric | 0.5 | 2.24+0j | 0.02–0.19 | 0.10 | 1e-10 | – |
| ZnS | dielectric | 5 | 2.26+0j | 0.01–0.16 | 0.10 | 1e-11 | – |
| SiO2_IR | phonon polariton | 20 | 1.9+0.171j | 0.00–0.24 | 0.29 | – | 0.00 |
| sapphire_IR | phonon polariton | 25 | 1.06+4.92j | 0.00–0.85 | 0.25 | – | 0.24 |
| SiC_film | phonon polariton | 20 | 3.36+1.01j | 0.04–0.33 | 0.15 | – | 0.16 |
| SiC | phonon polariton | 20 | 5.84+0.116j | 0.02–0.50 | 0.12 | – | 0.00 |

## 4. How materials enter the computation (and its limits)
- **Units.** Period d = 2π and c0 = 1, so a grating of physical period P (µm) has λ = P/ω. The bands q cover
  ω ∈ [q, q+1]/n^u. The MATLAB code used q = 1…6; **q = 0** (ω < 1, sub-wavelength grating) is where
  first-order surface plasmons live at normal incidence, so the new plasmonic presets include it.
- **Dispersion.** The AWE expansion in δ holds n fixed inside each window. `--index-mode band` (default)
  evaluates n at every window centre. `--max-delta 0.05` cuts the windows finer (|δ| ≤ 0.05), which gives
  piecewise-constant dispersion re-evaluated about every 10 % in frequency. `--index-mode fixed --lambda-ref`
  gives the thesis' single "representative value".
- **Joint windows** (`--windows joint`). For a real n^w > n^u, the lower layer has its own Rayleigh points
  n^w ω = |α + p| inside the bands, and the δ-series cannot converge across them. Splitting each band at those
  points (paper eq. (29) applied to both layers) improves the energy defect by 1–3 orders of magnitude:
  - n^w = 1.1: median D goes from 1e-8.5 to 1e-10.8.
  - GaN: median D goes from 1e-4.2 to 1e-5.3, and to 1e-10 at small ε.
- **Validation of the physics.** The reflectivity dips sit where the flat-interface SPP condition
  n_spp P/λ = 1 predicts:
  - Ag, P = 0.5 µm: the dip is at 0.5254 µm versus 0.5238 µm predicted (`tests::test_surface_plasmon_position`).
  - Gold under water, P = 0.6 µm: the dip is at 0.8252 µm versus 0.8256 µm predicted.
- **Resolution warnings.**
  - IR noble metals (`Au_IR`, `Ag_IR`: |n| ≈ 45) have skin depths far below the Chebyshev spacing at
    Nz = 32, which gives Padé spikes with R > 1 at about 0.1 % of the points. Use `--Nz 64` or a smaller `--b`.
  - Lossless high-index substrates need `--windows joint`.

## 5. Command line
```
python refl_map.py --list-materials                    # curated library with n,k at 0.4/0.6/1/1.55/10 um
python refl_map.py --list-scenarios                    # paper, thesis (Figs. 17-33) and material presets
python refl_map.py --scenario Ag_disp                  # silver SPP lines (q = 0..2)
python refl_map.py --scenario water_over_gold          # SPR biosensor grating
python refl_map.py --nw Na --period 0.5 --q 0 1 2 --max-delta 0.05
python refl_map.py --nu water --nw Au --profile cosx --period 0.6 --q 0 1
python refl_map.py --nw Si_IR --period 5 --windows joint
python refl_map.py --nw VO2_hot --period 5 --max-delta 0.05     # then VO2_cold: a switchable grating
python refl_map.py --nw 3.8313+2.9043i --profile sin4x --M 15    # literal index (thesis Fig. 20a)
python refl_map.py --mode TE --nw Cu --period 1 --profile sin5x --a 0.6366 --b 0.6366
python refl_map.py --nw main/Au/nk/Olmon-sc.yml --rii-db <path>/database/data --period 10
python material_survey.py --keys Ag Au Al Na TiN ITO   # screen any subset of the library
```
Other options: `--mode TE|TM`, `--alpha`, `--profile` (cosx, sinx, cos3x, sin3x, cos4x, sin4x, sin5x, fs1, fs2,
rough, lipschitz, rough120, lipschitz120, or `expr:<numpy expression in x>`), `--eps-max`, `--a`, `--b`,
`--M`, `--N`, `--Nx`, `--Nz`, `--summation taylor|pade`, `--band-center`, `--sigma`, `--neps`, `--ndelta`,
`--absolute`, `--tag`.
