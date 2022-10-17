# HOPS/AWE Grating Scattering
The scattering of electromagnetic radiation by a layered periodic diffraction grating is a central model in engineering and the sciences.  The
numerical simulation of this experiment has been widely explored in the literature and we advocate for a novel interfacial method which is perturbative in
nature. This repository contains a High–Order Perturbation of Surfaces/Asymptotic Waveform Evaluation (HOPS/AWE) algorithm which is well–suited for PDEs
posed on piecewise homogeneous domains. The algorithm was built and designed by Matthew Kehoe (mkehoe@mtu.edu) and David Nicholls at the University of Illinois at Chicago.

## Installation
Download the Matlab code in the src directory. The three main files are

1. refl_map.m: Calculates the Reflectivity Map, R, which measures the response (reflected energy) of a periodically corrugated grating structure as a
function of illumination frequency, ω, and corrugation amplitude, h.
2. mms_error.m: Validates the accuracy of the algorithm through the Method of Manufactured Solutions (MMS).
3. test_mms_error.m: Rigorously validates all of the important Matlab code.

## Plotting 

Example 1: The Reflectivity Map and Energy Defect (D) for vacuum over a dielectric.

![alt text](https://axion004.files.wordpress.com/2022/10/refl_map_vacuum_dielectric.png)

Example 2: The Reflectivity Map for vacuum over silver and gold.

![alt text](https://axion004.files.wordpress.com/2022/10/refl_map_vacuum_metals.png)

More examples (with instructions on how to run the code) can be found in the plots directory.

## References
Interested readers can review the following papers

* General Overview of the Algorithm: [AWE Paper](http://homepages.math.uic.edu/~nicholls/papers/Submitted/HOPSAWEComput.pdf)
* Theory behind the Algorithm: [Theory Paper](http://homepages.math.uic.edu/~nicholls/papers/Submitted/HOPSAWEAnal.pdf)
* Overview of the Method: [Wave Scattering Periodic Media](https://axion004.files.wordpress.com/2022/10/wave_scattering_hops_awe-2-1.pdf)

## Questions

Email Matthew Kehoe (mkehoe@mtu.edu) about any questions related to the HOPS/AWE algorithm or the above references. Feel free to request a Zoom/Skype meeting if you are interested in the material.
