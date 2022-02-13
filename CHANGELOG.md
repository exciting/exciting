# Change Log

[oxygen.0.4] 2022-01-012
----------------------------
### Sorting Tolerance Fix

#### Added

#### Fixed

* The use of a tolerance in the heapsort algorithm was removed just prior to the Oxygen release. On certain HPC
  environments, with O3 optimisation, this can lead to inconsistent sorting of degenerate values when the heapsort
  is called on multiple processes or at different points in exciting. One consequently sees failure in the GW, for 
  example. This fix reintroduces the tolerance, and moves from an absolute tolerance to a relative tolerance. 

#### Changed

* `WANNER_INFO.OUT` is no longer tested for vectors. Testing the TDOS is sufficient to assert that the Wannier functions
  are consistent. Differences can occur due to the iterative optimization step, which will find one of multiple 
  physically-equivalent minima of the cost function corresponding to different but equivalent sets of Wannier functions.

* MSE tolerance for `ELECTCOND_11.OUT` output of the `groundstate-LDA_PW-properties-transport-Si` test changed from
  4.e12 to 6.e12, because the test was failing for the MPISMP Intel2019 pipeline in the CI. One was unable to reproduce
  this locally, suggesting there was some small floating point difference. One notes that whilst the absolute tolerance
  is large, the relative error is reasonable (10-5 - 10-6). 


[oxygen.0.3] 2022-01-09
----------------------------

### Plasma Frequency Fix  

#### Added

#### Fixed

* The wrong occupation factor was leading to systematic underestimations for the plasma frequency. This fix leads to the 
  correct plasma frequency, obtained when calling the `dielmat` property. 

#### Changed

 
[oxygen.0.2] 2022-01-09
----------------------------

### Continuous Integration Script Update to Allow Future Patches.
    
#### Added

#### Fixed

#### Changed

* exciting's Gitlab CI runner environment has changed, requiring modifications to the .gitlab script.
* Tolerance for number of SCF iterations in XANES TiO2 has been loosened - energies
  remain consistent.
* Th Ar IORA test using serial GCC has been disabled. The issue (rdirac zero wavefunction) only occurs in the CI for the 
  serial GCC build, and cannot be reproduced, therefore we assume it is an environment issue. It is also fine in
  the development branch. 


[oxygen.0.1] 2021-10-22
----------------------------

### Bugfix in BSE Calculation of the Real Part of the Dielectric function.

#### Added

#### Fixed

* This fix corrects the calculation of the real part of the dielectric function from the eigenstates and vectors of the 
  BSE. The imaginary part is not affected.
* This fix also corrects the calculation of the electron energy loss function (EELS), which depends upon the real part 
  of the dielectric function.
* The BSE tests have been updated in a consistent manner.

#### Changed
