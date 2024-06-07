# HOPS/AWE Grating Scattering
The scattering of electromagnetic radiation by a layered periodic diffraction grating is a central model in engineering and the sciences.  The
numerical simulation of this experiment has been widely explored in the literature and we advocate for a novel interfacial method which is perturbative in
nature. This repository contains a High–Order Perturbation of Surfaces/Asymptotic Waveform Evaluation (HOPS/AWE) algorithm which is well–suited for PDEs
posed on piecewise homogeneous domains. The algorithm was built and designed by Matthew Kehoe (mskehoe001@gmail.com) and David Nicholls at the University of Illinois at Chicago.

## Installation
Download the Matlab code in the src directory. The three main files are

1. `refl_map.m`: Calculates the Reflectivity Map, R, which measures the response (reflected energy) of a periodically corrugated grating structure as a
function of illumination frequency, ω, and corrugation amplitude, h.
2. `mms_error.m`: Validates the accuracy of the algorithm through the Method of Manufactured Solutions (MMS).
3. `test_mms_error.m`: Rigorously validates all of the important Matlab code.

## Plotting 

Example 1: The Reflectivity Map and Energy Defect (D) for vacuum over a dielectric.

![alt text](https://axion004.files.wordpress.com/2022/10/refl_map_vacuum_dielectric.png)

Example 2: The Reflectivity Map for vacuum over silver and gold.

![alt text](https://axion004.files.wordpress.com/2022/10/refl_map_vacuum_metals.png)

More examples (with instructions on how to run the code) can be found in the plots directory.

## References
 
[1]  M. Kehoe and D. P. Nicholls, [*Joint Geometry/Frequency Analyticity of Fields Scattered by Periodic Layered Media*](https://epubs.siam.org/doi/10.1137/22M1477568), SIAM Journal on Mathematical Analysis, Volume 55, Issue 3, 1737-1765 (2023). https://doi.org/10.1137/22M1477568

[2] M. Kehoe and D. P. Nicholls, [*A stable high-order perturbation of surfaces/asymptotic waveform evaluation method for the numerical solution of grating scattering problems*](https://link.springer.com/article/10.1007/s10915-024-02566-6). Journal of Scientific Computing., 100(1) (2024). https://doi.org/10.1007/s10915-024-02566-6

[3] M. Kehoe, [*Joint Analyticity of the Transformed Field and Dirichlet–Neumann Operator in Periodic Media*](https://matthewshawnkehoe.github.io/files/kehoe_thesis.pdf), PhD Thesis, University of Illinois at Chicago (2022).


## Questions

Email Matthew Kehoe (mskehoe001@gmail.com) about any questions related to the HOPS/AWE algorithm or the above references.
