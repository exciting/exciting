module rttddft_file_names
  implicit none
  
  private

  !> Default name of the file where the vector potential is printed out
  character(len=*), public, parameter :: filename_avec = 'VECTOR_POTENTIAL'
  !> Default name of the file where the external electric field strength is printed out
  character(len=*), public, parameter :: filename_evec = 'ELECTRIC_FIELD'  
  !> Default name of the file where the current density is printed out
  character(len=*), public, parameter :: filename_jind = 'CURRENT'
  !> Default name of the file where the number of excited electrons is printed out
  character(len=*), public, parameter :: filename_nexc = 'N_EXCITATIONS'
  !> Default name of the file with general information about the RT-TDDFT calculation
  character(len=*), public, parameter :: filename_info = 'RTTDDFT_INFO'
  !> Default name of the file where `pmat` is printed out
  character(len=*), public, parameter :: filename_pmat = 'PMATBASIS'
  !> Default name of the file where `pmat_mt` is printed out
  character(len=*), public, parameter :: filename_pmat_mt = 'PMATMTBASIS'
  !> Default name of the file where the projection coefficients are printed out
  character(len=*), public, parameter :: filename_projection_coefficients = 'PROJECTION_COEFFS_'
  !> Default name of the file where the eigenvalues are printed out
  character(len=*), public, parameter :: filename_eigenvalues = 'EIGVAL_'
  !> Default name of the file where the occupation factors are printed out
  character(len=*), public, parameter :: filename_occupations = 'OCCSV_TXT_'
  !> Default name of the file where the (initial) electron density is printed out
  character(len=*), public, parameter :: filename_density = 'RHO3D'
  !> Default name of the file where the changes in electron density are printed out
  character(len=*), public, parameter :: filename_density_changes = 'DELTARHO3D'
  !> Default name of the file to store the core and valence densities and the KS potential
  character(len=*), public, parameter :: filename_rho_vks = "RHO_VKS"
  !> Typical suffix to differentiate RT-TDDFDT files from ground state files
  character(len=*), public, parameter :: RTTDDFT_suffix = "_RTTDDFT"
  !> Suffix referring to groundstate
  character(len=*), public, parameter :: GND_sufix = "_GND"
  !> Suffix used when performing single-shot ground state calculation (before RT-TDDFDT)
  character(len=*), public, parameter :: RTTDDFT_GND_sufix = RTTDDFT_suffix // GND_sufix
  !> Default name of the file where there wavefunction coefficients are printed out
  character(len=*), public, parameter :: filename_wavefunction = 'EVECFV'
  !> Default name of the file where there wavefunction coefficients are printed out
  character(len=*), public, parameter :: filename_wavefunction_second_variation = 'EVECSV' 
  !> Descriptors name used to write wavefunctions into an output file
  character(len=*), public, parameter :: kpt_latt_name = "kpoints_lattice_coord"
  !> Suffix for file where there wavefunction coefficients \(\psi(t)\) are printed out
  character(len=*), public, parameter :: suffix_wavefunction_t = RTTDDFT_suffix
  !> Suffix for file where there wavefunction coefficients \(\psi(t-\Delta t)\) are printed out
  character(len=*), public, parameter :: suffix_wavefunction_t_minus_dt = '_PREVIOUS' // RTTDDFT_suffix
  !> Suffix for where there groundstate wavefunction coefficients are printed out
  character(len=*), public, parameter :: suffix_wavefunction_gnd = RTTDDFT_GND_sufix
  !> Default name of the file where timigs are printed out
  character(len=*), public, parameter :: filename_timing = 'TIMING' // RTTDDFT_suffix
  !> Default name of the file where the total energy is printed out
  character(len=*), public, parameter :: filename_etot = 'TOTENERGY' // RTTDDFT_suffix
  !> Default name of the file where the string phases are printed out
  character(len=*), public, parameter :: filename_phases = 'PHASES' // RTTDDFT_suffix
  !> Default name of the file where the polarization is printed out
  character(len=*), public, parameter :: filename_pvec = 'POLARIZATION' // RTTDDFT_suffix
end module