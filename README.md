NMR MOLECULAR EIGENMODES SOLVER
================================

This project contains Python codes for modeling molecular eigenmodes in the NMR relaxation of fluids.
The project contains two sets of codes: one in LJ- and NMR-reduced units, and another one in real SI units (International System of Units).
The codes numerically solve the Fokker-Planck equation for dipole pairs to compute multiexponential
time decays, amplitudes, autocorrelation functions, and NMR frequency dispersions.
The codes account for radial relative diffusion and potential of interaction between the dipoles.
Each code can be used independently or in combination for comprehensive evaluations.


Codes Overview
--------------

Code A | Radial Component of Diffusion Propagator
- Reads input from input.txt.
- If the code is executed in potential_def = "read", it also reads rdf.txt.
- Computes numerically the radial component of the diffusion propagator for two magnetic dipoles in real units.
- Solves the Fokker-Planck equation and the diffusion propagator.
- Calculates multiexponential time decays (τ's), amplitudes P(τ)'s, the dipole-dipole NMR
  autocorrelation function, and 1/T1(0) = 1/T2(0).
- Can be used alone or together with Codes B and C.

Code B | Multi-case Evaluation of τ's and P(τ)'s in different spherical shells
- Computes multiexponential time decays and amplitudes for different ranges of ri/rf (different spherical shells), optionally icluding interaction potentials.
- Reads and modifies input.txt to adjust ri/rf values.
- Designed to be used with Code A for evaluating multiple cases.

Code C | NMR Frequency Dispersion
- Reads input.txt.
- Computes dipole-dipole NMR frequency dispersion for 1/T1 and 1/T2 for single cases of dipole
  pairs with either like- or unlike-spins.
- Used together with Code A to evaluate NMR frequency-dependent properties.


Best Practices
--------------

To obtain reliable and converged results, users should follow three essential best practices:

(1) Radial discretization (dr): Ensure that the discretization step in the radial direction is sufficiently small so that results do not change appreciably when dr is further decreased. This guarantees that the radial grid is fine enough for accurate calculations.

(2) Eigenvalue range (vector x): Define an appropriate search range for the eigenvalues. The correctness of this choice can be verified by checking that the code captures the lowest eigenmode first (i.e., the highest characteristic time τ), whose eigenfunction R(r) has no nodes (i.e., R(r) does not change sign over the evaluated interval). In general, the l-th mode will exhibit l-1 nodes in its eigenfunction.

(3) Eigenvalue discretization (dx): Use a sufficiently small discretization step in the eigenvalue search so that results remain stable when dx is further decreased. This ensures that the eigenvalue grid is fine enough for accurate convergence.

(4) Identification of spherical regions: Each spherical shell assumes that all molecules within it share the same average diffusivity; if diffusivity changes significantly beyond a certain distance, a new region must be defined for those dipoles. Notice that the spherical definition with one dipole in the center also applies to dipole pairs within the same molecule (e.g., in a linear hydrocarbon).


Requirements
------------

- Python 3.x
- Required libraries: numpy, scipy, matplotlib
- Install dependencies via pip: pip install numpy scipy matplotlib


Usage
-----

Running the codes:
- python A_Numeric_radial_FP_realunits.py
- python B_Range_analysis.py
- python C_NMR_disp.py

Running the examples:
- bash run_example.sh


File Structure
--------------

- Code A      : Computes radial component of diffusion propagator and τ/P(τ) for a single case.
- Code B      : Computes τ/P(τ) for multiple ri/rf values.
- Code C      : Computes NMR frequency dispersion for dipole pairs.
- input.txt   : Input data file used by all codes.
- rdf.txt     : Input data file with r versus pair distribution function g(r) used in code A.
- README.txt  : This file.


Citation
--------

If you use these codes, please cite:

[1] T. J. Pinheiro dos Santos, B. Orcan-Ekmecki, W. G. Chapman, P. M. Singer, and D. N. Asthagiri,
"Extended molecular eigenmodes treatment of dipole-dipole NMR relaxation in real fluids," J. X., vol. X, p. X, 2025.

[2] T. J. Pinheiro dos Santos, B. Orcan-Ekmecki, W. G. Chapman, P. M. Singer, and D. N. Asthagiri,
"Theory and modeling of molecular modes in the NMR relaxation of fluids," J. Chem. Phys., vol. 064108, p. 160, 2024.
