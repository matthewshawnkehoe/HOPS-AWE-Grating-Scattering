| scenario | profile | n_w | summation | grid | windows x map grid | corrected points | max \|R_hyb - R_AWE\| | max log10\|D\| AWE -> hybrid (lossless) | AWE / hybrid time |
|---|---|---|---|---|---|---|---|---|---|
| silver | cos4xcos4y | (0.05+2.275j) | Pade 12/12 | 32x32x16 | 19 x 15^2 | 82 % | 3.1e-02 | absorbing (A ~ 0.05) | 52 / 2294 s |
| gold | cos4xcos4y | (1.48+1.883j) | Pade 12/12 | 32x32x16 | 19 x 15^2 | 78 % | 1.1e-02 | absorbing (A ~ 0.67) | 55 / 2458 s |
| silver_lattice | cos4xcos4y | (0.05+2.275j) | Pade 12/12 | 32x32x16 | 7 x 21^2 | 79 % | 8.3e-02 | absorbing (A ~ 0.04) | 19 / 1425 s |
| dielectric | cosxcosy | 1.1 | Taylor 14/14 | 16x16x16 | 8 x 21^2 | 52 % | 6.9e-05 | 0.2 -> -10.2 | 21 / 1408 s |
| dielectric_joint | cosxcosy | 1.1 | Taylor 14/14 | 16x16x16 | 15 x 15^2 | 20 % | 8.2e-06 | -4.4 -> -9.4 | 41 / 789 s |
| dielectric_alpha | cosxcosy | 1.1 | Taylor 14/14 | 16x16x16 | 8 x 21^2 | 53 % | 1.8e-04 | 0.9 -> -6.2 | 25 / 1589 s |
| egg_silver | egg | (0.05+2.275j) | Pade 12/12 | 16x16x16 | 8 x 21^2 | 68 % | 3.6e-03 | absorbing (A ~ 0.03) | 16 / 1046 s |
| sum_gold | cos4x+cos4y | (1.48+1.883j) | Pade 12/12 | 32x32x16 | 8 x 21^2 | 54 % | 8.0e-03 | absorbing (A ~ 0.65) | 23 / 1345 s |
| silver_1d | cos4x | (0.05+2.275j) | Pade 15/15 | 32x1x32 | 6 x 21^2 | 76 % | 3.1e-01 | absorbing (A ~ 0.04) | 7 / 461 s |
| dielectric_1d | cosx | 1.1 | Taylor 16/16 | 32x1x32 | 6 x 21^2 | 67 % | 6.3e-05 | -4.2 -> -10.5 | 5 / 309 s |
| conical_1d | cosx | 1.1 | Taylor 14/14 | 32x1x32 | 7 x 21^2 | 24 % | 9.8e-06 | -4.8 -> -9.2 | 10 / 312 s |
| Au_crossed | cosxcosy | Au | Pade 12/12 | 16x16x16 | 21 x 15^2 | 16 % | 4.0e-03 | absorbing (A ~ 0.16) | 41 / 662 s |
| Ag_crossed | cosx+cosy | Ag | Pade 12/12 | 16x16x16 | 21 x 15^2 | 24 % | 3.5e-03 | absorbing (A ~ 0.07) | 39 / 766 s |
| Al_uv_crossed | cosxcosy | Al | Pade 12/12 | 16x16x16 | 21 x 15^2 | 29 % | 3.6e-03 | absorbing (A ~ 0.08) | 37 / 788 s |
| water_over_gold_crossed | cosxcosy | Au | Pade 12/12 | 16x16x16 | 23 x 15^2 | 26 % | 2.7e-02 | absorbing (A ~ 0.67) | 21 / 418 s |
| TiO2_crossed | cosxcosy | TiO2 | Pade 12/12 | 16x16x16 | 21 x 15^2 | 34 % | 4.1e-04 | -2.9 -> -8.2 | 19 / 633 s |
| Si_crossed | cosxcosy | Si | Pade 12/12 | 16x16x16 | 14 x 15^2 | 76 % | 1.7e-03 | absorbing (A ~ 0.07) | 33 / 1966 s |
| SiC_reststrahlen_crossed | cosxcosy | SiC | Pade 12/12 | 16x16x16 | 97 x 9^2 | 3 % | 2.0e-06 | absorbing (A ~ 0.00) | 84 / 639 s |
| oblique_gold | cosxcosy | (1.48+1.883j) | Pade 12/12 | 16x16x16 | 19 x 15^2 | 36 % | 5.2e-04 | absorbing (A ~ 0.64) | 41 / 1419 s |
| silver_1d_nx64 | cos4x | (0.05+2.275j) | Pade 15/15 | 64x1x32 | 6 x 21^2 | 76 % | 1.2e-01 | absorbing (A ~ 0.04) | 11 / 634 s |
