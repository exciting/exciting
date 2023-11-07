project: The exciting FP-LAPW Code  
summary: An all-electron (L)APW+lo package for DFT and excited-state calculations.  
author: SOL group  
project_website: http://exciting-code.org  
project_download: http://exciting-code.org/download-exciting  
project_github: https://github.com/exciting/exciting  
predocmark: >  
docmark: !  
display: public  
         protected
         private
graph: false  
warn: false  
search: false  
proc_internals: false  
src_dir: ../src/src_gw  
         ../src/src_hybrids
         ../src/input
         ../src/testframework
         ../src/lapack_wrappers
         ../src/math 
         ../src/constants
         ../src/errors_warnings
         ../src/src_mpi
         ../src/src_xs/src_rttddft/
         ../src/simplified_input
         ../src/structure
         ../src/wavefunctions
         ../src/xstring
         ../src/writers
         ../src/potential
         ../src/matrix_elements
         ../src/src_xc
         ../src/basis
         ../src/file_io
         ../src/dfpt
         ../src/src_phonon
         ../src/matrix_fourier_interpolation
         ../src/xgrid/
         ../src/xhdf5/
         ../src/xhdf5/hdf5_wrappers/ 
	       ../src/svlo
         ../src/isdf
         ../src/src_xs/fastBSE/
         ../src/src_xs/fastBSE/BSH/
         ../src/src_xs/fastBSE/utils/

[//]: # "Note, ford commands can not be separated by whitelines."  
[//]: # "More information on ford's project file options can be found at:"  
[//]: # "https://github.com/Fortran-FOSS-Programmers/ford/wiki/Project-File-Options"  

Release: exciting fluorine.  
Copyright (C) 2002-2021 The exciting team.  

## Description

exciting is an all-electron full-potential computer package for
first-principles calculations, based on (linearized) augmented
planewave + local orbital [(L)APW+lo] methods. This family of
basis sets is known as the most precise numerical scheme to
solve the Kohn-Sham equations of density-functional theory (DFT),
reaching up to micro-Hartree precision [[1](#citing-exciting-in-your-work)].

## Documentation for developers

The documentation for developers can be found in the [exciting gitlab wiki](https://git.physik.hu-berlin.de/sol/exciting/-/wikis/home).

## Citing exciting in Your Work

If you use exciting in your work, please cite the following paper:

[1]: [exciting — a full-potential all-electron package implementing density-functional theory and many-body perturbation theory. Andris Gulans *et al* 2014 J. Phys.: Condens. Matter 26 363202.](https://doi.org/10.1088/0953-8984/26/36/363202)

Additionally, if you use specific features, please cite the corresponding papers where applicable:  

[2]: [Probing the LDA-1/2 method as a starting point for G0W0 calculations. Phys. Rev. B **94**, 235141 (2016).](https://doi.org/10.1103/PhysRevB.94.235141)

[3]: [Time-dependent density functional theory versus Bethe–Salpeter equation: an all-electron study. ***Phys. Chem. Chem. Phys.***, 2009, **11**, 4451-4457.](https://doi.org/10.1039/B903676H)

[4]: [All-electron full-potential implementation of real-time TDDFT in exciting. Electron. Struct. **3**, 037001 (2021).](https://doi.org/10.1088/2516-1075/ac0c26)

[5]: [Accurate all-electron G0W0 quasiparticle energies employing the full-potential augmented plane-wave method. Phys. Rev. B **94**, 035118 (2016).](https://doi.org/10.1103/PhysRevB.94.035118)

[6]: [Addressing electron-hole correlation in core excitations of solids: An all-electron many-body approach from first principles. Phys. Rev. B **95**, 155121 (2017).](https://doi.org/10.1103/PhysRevB.95.155121)
