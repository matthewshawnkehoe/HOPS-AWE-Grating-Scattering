Errors against the converged reference (pointwise HOPS at delta = 0, N = 24, Pade, 2 Nx, 2 Ny, Nz + 16) -- total error -- and against the same pointwise HOPS on the scenario's own grid (same Nx, Ny, Nz) -- the SUMMATION error, which is what the PINN step addresses.

| scenario | n_w | grid | samples | total max err R: AWE (refl_map_3D) | AWE full Taylor | **hybrid** | summation max err R: AWE / **hybrid** | median total err R: AWE / hybrid | max err D: AWE / hybrid | hybrid better (> 2x) at | Wilcoxon p |
|---|---|---|---|---|---|---|---|---|---|---|---|
| silver | (0.05+2.275j) | 32x32x16 | 152 | 5.8e-03 | 1.0e+15 | **1.2e-03** | 5.9e-03 / **9.8e-06** | 4.7e-07 / 4.4e-07 | 5.8e-03 / 1.2e-03 | 22 % | 3e-04 |
| gold | (1.48+1.883j) | 32x32x16 | 152 | 2.5e-04 | 1.2e+01 | **1.5e-04** | 1.6e-04 / **1.6e-07** | 3.8e-07 / 3.9e-07 | 2.5e-04 / 1.5e-04 | 14 % | 3e-03 |
| silver_lattice | (0.05+2.275j) | 32x32x16 | 56 | 3.0e-03 | 5.2e+32 | **1.1e-03** | 2.0e-03 / **1.3e-05** | 1.3e-07 / 1.1e-07 | 3.0e-03 / 1.1e-03 | 27 % | 7e-02 |
| dielectric | 1.1 | 16x16x16 | 64 | 1.0e-05 | 1.2e+07 | **1.2e-12** | 1.0e-05 / **1.2e-12** | 5.7e-09 / 8.2e-16 | 6.8e-04 / 1.9e-12 | 100 % | 4e-12 |
| dielectric_joint | 1.1 | 16x16x16 | 120 | 1.8e-06 | 1.2e-07 | **2.3e-11** | 1.8e-06 / **2.3e-11** | 2.5e-10 / 1.2e-15 | 3.1e-06 / 3.1e-11 | 100 % | 2e-21 |
| dielectric_alpha | 1.1 | 16x16x16 | 64 | 1.4e-05 | 3.3e+10 | **2.3e-10** | 1.4e-05 / **2.3e-10** | 4.9e-09 / 2.0e-14 | 2.7e-03 / 4.4e-09 | 100 % | 4e-12 |
| egg_silver | (0.05+2.275j) | 16x16x16 | 64 | 1.2e-04 | 4.4e+25 | **5.9e-08** | 1.2e-04 / **5.9e-08** | 1.6e-08 / 1.5e-11 | 1.2e-04 / 5.9e-08 | 80 % | 9e-11 |
| sum_gold | (1.48+1.883j) | 32x32x16 | 64 | 1.0e-05 | 1.6e-02 | **6.8e-07** | 1.0e-05 / **3.0e-09** | 3.2e-09 / 3.1e-09 | 1.0e-05 / 6.8e-07 | 11 % | 2e-01 |
| silver_1d | (0.05+2.275j) | 32x1x32 | 48 | 3.7e-03 | 6.0e+39 | **7.2e-04** | 3.7e-03 / **6.9e-05** | 7.2e-07 / 3.3e-07 | 3.7e-03 / 7.2e-04 | 23 % | 1e-01 |
| dielectric_1d | 1.1 | 32x1x32 | 48 | 3.4e-05 | 4.1e+12 | **1.4e-12** | 3.4e-05 / **1.4e-12** | 1.3e-08 / 5.7e-16 | 3.4e-05 / 4.9e-13 | 100 % | 7e-15 |
| conical_1d | 1.1 | 32x1x32 | 56 | 1.4e-06 | 8.1e-08 | **9.1e-13** | 1.4e-06 / **9.1e-13** | 1.0e-09 / 1.2e-15 | 1.4e-06 / 3.7e-11 | 100 % | 8e-11 |
| Au_crossed | Au | 16x16x16 | 168 | 9.7e-06 | 3.6e+08 | **1.7e-08** | 9.7e-06 / **1.7e-08** | 1.6e-13 / 1.6e-13 | 9.7e-06 / 1.7e-08 | 26 % | 2e-01 |
| Ag_crossed | Ag | 16x16x16 | 168 | 7.5e-06 | 4.1e+38 | **1.2e-07** | 7.5e-06 / **1.2e-07** | 8.7e-14 / 1.5e-14 | 7.5e-06 / 1.2e-07 | 45 % | 9e-07 |
| Al_uv_crossed | Al | 16x16x16 | 168 | 1.0e-05 | 7.8e+05 | **7.6e-09** | 1.0e-05 / **7.6e-09** | 3.0e-13 / 1.5e-14 | 1.0e-05 / 7.6e-09 | 51 % | 3e-07 |
| water_over_gold_crossed | Au | 16x16x16 | 184 | 2.3e-04 | 3.7e+27 | **2.2e-09** | 2.3e-04 / **2.2e-09** | 4.3e-13 / 3.2e-14 | 2.3e-04 / 2.2e-09 | 49 % | 1e-06 |
| TiO2_crossed | TiO2 | 16x16x16 | 168 | 2.3e-05 | 3.3e-01 | **6.2e-09** | 2.3e-05 / **5.3e-09** | 3.6e-12 / 4.0e-13 | 1.5e-04 / 1.2e-08 | 38 % | 1e-01 |
| Si_crossed | Si | 16x16x16 | 112 | 1.6e-05 | 3.2e-05 | **5.2e-05** | 1.4e-05 / **4.5e-05** | 1.6e-09 / 1.1e-09 | 2.0e-05 / 5.2e-05 | 30 % | 4e-02 |
| SiC_reststrahlen_crossed | SiC | 16x16x16 | 104 | 5.7e-10 | 2.1e-01 | **4.8e-08** | 5.7e-10 / **4.8e-08** | 1.0e-15 / 8.9e-16 | 5.7e-10 / 4.8e-08 | 13 % | 5e-01 |
| oblique_gold | (1.48+1.883j) | 16x16x16 | 152 | 1.4e-06 | 2.5e-05 | **4.3e-09** | 1.4e-06 / **4.3e-09** | 1.1e-10 / 1.7e-11 | 1.4e-06 / 4.3e-09 | 49 % | 9e-05 |
| silver_1d_nx64 | (0.05+2.275j) | 64x1x32 | 48 | 2.2e-03 | 6.0e+39 | **2.1e-06** | 2.2e-03 / **2.1e-06** | 1.5e-08 / 8.1e-11 | 2.2e-03 / 2.1e-06 | 79 % | 9e-09 |

Pooled over 2160 samples (20 scenarios): max error AWE 5.8e-03, hybrid 1.2e-03; median AWE 1.1e-10, hybrid 4.7e-13; hybrid better by > 2x at 47 % of the samples, worse by > 2x at 15 %.
Summation error alone (against HOPS on the same grid): max AWE 5.9e-03, hybrid 6.9e-05; median AWE 7.8e-11, hybrid 2.8e-13; hybrid better by > 2x at 57 %, worse by > 2x at 18 %.
Excluding round-off (samples where the larger of the two errors is below 1e-12): hybrid better by > 2x at 40 %, worse by > 2x at 11 % of all samples (total error); errors below 1e-12 for both at 27 %.
Indicator calibration (full-order AWE Taylor sum vs its summation error): Spearman rho = 0.91; log10|R error| ~ 0.63 log10(indicator) -0.97.

| indicator tolerance | 1e-20 | 1e-16 | 1e-14 | 1e-12 | 1e-10 | 1e-08 | 1e-06 |
|---|---|---|---|---|---|---|---|
| fraction of samples solved | 76 % | 59 % | 51 % | 43 % | 35 % | 30 % | 22 % |
| max summation error | 6.9e-05 | 6.9e-05 | 6.9e-05 | 6.9e-05 | 6.9e-05 | 6.9e-05 | 8.0e-04 |
| median summation error | 6.7e-14 | 1.3e-13 | 2.8e-13 | 8.1e-13 | 2.2e-12 | 4.2e-12 | 7.6e-12 |

Basis order (max summation error of the least-squares weights, every sample solved):

| scenario | order 6 | order 8 | order 10 | order 12 |
|---|---|---|---|---|
| silver | 7.5e-03 | 8.6e-04 | 4.0e-05 | 9.8e-06 |
| gold | 4.0e-04 | 3.8e-05 | 9.3e-07 | 1.6e-07 |
| silver_lattice | 7.1e-03 | 7.0e-04 | 6.3e-05 | 1.3e-05 |
| dielectric | 2.3e-09 | 3.9e-11 | 7.1e-13 | 6.8e-14 |
| dielectric_joint | 1.3e-09 | 7.0e-12 | 8.7e-14 | 2.5e-14 |
| dielectric_alpha | 3.0e-07 | 7.3e-09 | 5.8e-10 | 2.3e-10 |
| egg_silver | 5.9e-07 | 1.3e-07 | 1.4e-07 | 5.9e-08 |
| sum_gold | 1.1e-05 | 1.8e-07 | 1.3e-08 | 2.5e-09 |
| silver_1d | 6.0e-03 | 2.2e-04 | 7.3e-05 | 6.9e-05 |
| dielectric_1d | 8.7e-08 | 1.0e-09 | 3.5e-12 | 1.2e-14 |
| conical_1d | 2.7e-09 | 9.2e-13 | 2.7e-15 | 4.4e-15 |
| Au_crossed | 7.8e-08 | 3.2e-09 | 5.2e-11 | 6.2e-11 |
| Ag_crossed | 9.3e-08 | 1.3e-07 | 1.2e-07 | 1.2e-07 |
| Al_uv_crossed | 5.9e-07 | 4.5e-09 | 2.8e-10 | 4.1e-12 |
| water_over_gold_crossed | 2.4e-07 | 1.5e-09 | 1.3e-10 | 4.9e-11 |
| TiO2_crossed | 8.2e-07 | 1.9e-08 | 8.6e-09 | 5.7e-09 |
| Si_crossed | 5.1e-05 | 3.4e-05 | 4.4e-05 | 4.5e-05 |
| SiC_reststrahlen_crossed | 1.1e-06 | 3.3e-08 | 1.0e-09 | 1.4e-11 |
| oblique_gold | 1.9e-06 | 2.0e-08 | 1.9e-09 | 1.2e-09 |
| silver_1d_nx64 | 1.1e-02 | 3.2e-04 | 1.4e-05 | 2.1e-06 |

Residual form (order 12, every sample solved): hops3d's discrete flux form (default) vs the chain-rule form of the 2D hybrid:

| scenario | max summation err R: discrete | chain rule | median indicator: discrete | chain rule |
|---|---|---|---|---|
| silver | 9.8e-06 | 1.9e-02 | 3.2e-06 | 3.6e-04 |
| gold | 1.6e-07 | 1.1e-04 | 3.5e-08 | 2.2e-05 |
| silver_lattice | 1.3e-05 | 2.8e-02 | 3.1e+03 | 3.1e+03 |
| dielectric | 6.8e-14 | 1.6e-12 | 6.0e-08 | 6.0e-08 |
| dielectric_joint | 2.5e-14 | 1.5e-12 | 3.4e-19 | 3.6e-18 |
| dielectric_alpha | 2.3e-10 | 2.3e-10 | 4.1e-07 | 4.1e-07 |
| egg_silver | 5.9e-08 | 5.8e-08 | 8.3e-04 | 8.3e-04 |
| sum_gold | 2.5e-09 | 9.9e-06 | 2.1e-12 | 4.1e-07 |
| silver_1d | 6.9e-05 | 1.3e-02 | 2.0e-02 | 3.6e-02 |
| dielectric_1d | 1.2e-14 | 1.2e-14 | 2.6e-04 | 2.6e-04 |
| conical_1d | 4.4e-15 | 4.6e-15 | 5.0e-19 | 5.0e-19 |
| Au_crossed | 6.2e-11 | 6.2e-10 | 2.6e-20 | 2.8e-19 |
| Ag_crossed | 1.2e-07 | 1.2e-07 | 2.9e-20 | 3.7e-20 |
| Al_uv_crossed | 4.1e-12 | 3.3e-10 | 8.6e-19 | 5.8e-17 |
| water_over_gold_crossed | 4.9e-11 | 4.8e-10 | 3.6e-19 | 7.2e-17 |
| TiO2_crossed | 5.7e-09 | 5.7e-09 | 1.0e-15 | 7.5e-15 |
| Si_crossed | 4.5e-05 | 4.5e-05 | 2.1e-08 | 2.1e-08 |
| SiC_reststrahlen_crossed | 1.4e-11 | 1.3e-09 | 3.4e-24 | 1.8e-21 |
| oblique_gold | 1.2e-09 | 1.2e-09 | 5.0e-16 | 9.5e-16 |
| silver_1d_nx64 | 2.1e-06 | 1.2e-05 | 2.7e-02 | 2.7e-02 |
