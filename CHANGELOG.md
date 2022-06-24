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
