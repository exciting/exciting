
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modtddft
! !DESCRIPTION: 
!   Global variables for the {\tt TDDFT} implementation
!   in the {\tt EXCITING}-code.
!
!
! !REVISION HISTORY: 
!
!  Created June 2006 (SAG)
! !PUBLIC DATA MEMBERS:
  implicit none

  !---------------!
  !     NOTES     !
  !---------------!
  ! use the following units for standar output and standard error:
  ! -1: standard output
  !  0: standard error

  !---------------------------------------------------------------!
  !     muffin-tin radial mesh and angular momentum variables     !
  !---------------------------------------------------------------!
  integer, allocatable :: idxxlm(:,:)

  !------------------------------!
  !     q-point set variables    !
  !------------------------------!
  ! input type for q-point set
  character(256) :: qtype
  ! file of q-points for input
  character(256) :: qlist
  ! q-point offset
  real(8) :: vqloff(3)
  ! true if the eigenvectors for the Gamma point are to be calculated
  logical :: tq0ev
  ! true if the first q-point is the Gamma point
  logical :: tq1gamma

  !----------------------------------!
  !     G+q-vector set variables     !
  !----------------------------------!
  ! maximum |G+q| cut-off for APW functions
  real(8) gqmax
  ! number of G+q-vectors
  integer, allocatable :: ngq(:)
  ! maximum number of G+q-vectors over all q-points
  integer ngqmax
  ! index from G+q-vectors to G-vectors
  integer, allocatable :: igqig(:,:)
  ! G+q-vectors in lattice coordinates
  real(8), allocatable :: vgql(:,:,:)
  ! G+q-vectors in Cartesian coordinates
  real(8), allocatable :: vgqc(:,:,:)
  ! length of G+q-vectors
  real(8), allocatable :: gqc(:,:)
  ! (theta, phi) coordinates of G+q-vectors
  real(8), allocatable :: tpgqc(:,:,:)
  ! structure factor for the G+q-vectors
  complex(8), allocatable :: sfacgq(:,:,:)
  ! spherical harmonics of the G-vectors
  complex(8), allocatable :: ylmgq(:,:,:)

  !--------------------------------------!
  !     k-point set  variables (q=0)     !
  !--------------------------------------!
  ! k-points in lattice coordinates
  real(8), allocatable :: vkl0(:,:)

  !-------------------------!
  !     k+q-point set       !
  !-------------------------!
  ! offset for k+q-point set derived from q-point
  real(8),allocatable :: qvkloff(:,:)
  ! map from k-point index to k+q point index for same k
  integer, allocatable :: ikmapikq(:,:)

  !-----------------------------------------!
  !     G+k-vector set  variables (q=0)     !
  !-----------------------------------------!
  ! number of G+k-vectors for augmented plane waves
  integer, allocatable :: ngk0(:,:)
  ! maximum number of G+k-vectors over all k-points
  integer ngkmax0
  ! index from G+k-vectors to G-vectors
  integer, allocatable :: igkig0(:,:,:)
  ! G+k-vectors in lattice coordinates
  real(8), allocatable :: vgkl0(:,:,:,:)
  ! G+k-vectors in Cartesian coordinates
  real(8), allocatable :: vgkc0(:,:,:,:)
  ! length of G+k-vectors
  real(8), allocatable :: gkc0(:,:,:)
  ! (theta, phi) coordinates of G+k-vectors
  real(8), allocatable :: tpgkc0(:,:,:,:)
  ! structure factor for the G+k-vectors
  complex(8), allocatable :: sfacgk0(:,:,:,:)

  !---------------------------------------!
  !     Hamiltonian and APW variables     !
  !---------------------------------------!
  ! maximum nmat over all k-points (q=0)
  integer nmatmax0
  ! order of overlap and Hamiltonian matrices for each k-point (q=0)
  integer, allocatable :: nmat0(:,:)
  ! first-variational eigenvectors
  complex(8), allocatable :: evecfv(:,:,:)
  ! second variational eigenvectors
  complex(8), allocatable :: evecsv(:,:)
  ! first-variational eigenvectors (q=0)
  complex(8), allocatable :: evecfv0(:,:,:)
  ! first variational eigenvalues
  real(8), allocatable :: evalfv(:,:)
  ! second-variational eigenvalues
  real(8), allocatable :: evalsv0(:,:)
  ! matching coefficients
  complex(8), allocatable :: apwalm(:,:,:,:)
  ! matching coefficients (q=0)
  complex(8), allocatable :: apwalm0(:,:,:,:)
  ! expansion coefficients of APW functions
  complex(8), allocatable :: apwdlm(:,:,:,:)
  ! expansion coefficients of APW functions (q=0)
  complex(8), allocatable :: apwdlm0(:,:,:,:)
  ! APW coefficients for muffin-tin part of the wavefunction
  complex(8), allocatable :: wfcmt(:,:,:,:)
  ! APW coefficients for muffin-tin part of the wavefunction (q=0)
  complex(8), allocatable :: wfcmt0(:,:,:,:)

  !--------------------------------------------!
  !     eigenvalue and occupancy variables     !
  !--------------------------------------------!
  ! number of occupied valence states (valence band states)
  integer :: nstval
  ! number of unoccupied valence states (conduction band states)
  integer :: nstcon
  ! eigenvalue differences (resonant part)
  real(8), allocatable :: deou(:,:)
  ! eigenvalue differences (anti-resonant part)
  real(8), allocatable :: deuo(:,:)

  !--------------------------------------------------!
  !     matrix elements of exponential expression    !
  !--------------------------------------------------!
  ! maximum angular momentum for Rayleigh expansion of exponential
  integer :: lmaxemat
  ! (lmaxemat+1)^2
  integer :: lmmaxemat
  ! maximum angular momentum for APW functions (for matrix elements)
  integer :: lmaxapwtd
  ! (lmaxapwtd+1)^2
  integer :: lmmaxapwtd
  ! Gaunt coefficients array
  real(8), allocatable :: tdgnt(:,:,:)
  ! radial integrals coefficients (APW-APW)
  complex(8), allocatable :: intrgaa(:,:,:,:,:)
  ! radial integrals coefficients (lo-APW)
  complex(8), allocatable :: intrgloa(:,:,:,:,:)
  ! radial integrals coefficients (APW-lo)
  complex(8), allocatable :: intrgalo(:,:,:,:,:)
  ! radial integrals coefficients (lo-lo)
  complex(8), allocatable :: intrglolo(:,:,:,:,:)
  ! radial integrals (APW-APW)
  real(8), allocatable :: riaa(:,:,:,:,:,:,:)
  ! radial integrals (lo-APW)
  real(8), allocatable :: riloa(:,:,:,:,:,:)
  ! radial integrals (lo-lo)
  real(8), allocatable :: rilolo(:,:,:,:,:)
  ! helper matrix
  complex(8), allocatable :: xih(:,:)
  ! helper matrix
  complex(8), allocatable :: xihir(:,:)
  ! helper matrix
  complex(8), allocatable :: xiohalo(:,:)
  ! helper matrix
  complex(8), allocatable :: xiuhalo(:,:)
  ! helper matrix 
  complex(8), allocatable :: xiohloa(:,:)
  ! helper matrix 
  complex(8), allocatable :: xiuhloa(:,:)
  ! helper matrix
  complex(8), allocatable :: xihlolo(:,:)
  ! matrix elements array (resonant part)
  complex(8), allocatable :: xiou(:,:,:)
  ! matrix elements array (anti-resonant part)
  complex(8), allocatable :: xiuo(:,:,:)

  !---------------------------------!
  !     momentum matrix elements    !
  !---------------------------------!
  ! momentum matrix elements (resonant part)
  complex(8), allocatable :: pmou(:,:,:)
  ! momentum matrix elements (anti-resonant part)
  complex(8), allocatable :: pmuo(:,:,:)

  !------------------------------------------!
  !     response and dielectric functions    !
  !------------------------------------------!
  ! k-point step size for checkpointing patial sums of response functions
  integer :: kstepdf ! *** to be done ***
  ! type of response function (time-ordered/retarded/advanced)
  character(4) :: rsptype
  ! true if analytic continuation to the real axis is to be performed
  logical :: acont
  ! number of energy intervals
  integer :: nwdf
  ! number of energy intervals (on imaginary axis) for analytic continuation
  integer :: nwacont
  ! broadening for Kohn Sham response function
  real(8) :: brdtd
  ! true if to consider the anti-resonant part for the dielectric function
  logical :: aresdf
  ! true if only diagonal part of local field effects is considered
  logical :: lfediag
  ! true if wings of dielectric function use symmetrized momentum matr. el.
  logical :: symwings
  ! symmetrization matrix for q=0 dielectric tensor
  real(8) :: symdfq0(3,3)
  ! true on occurrance of negative 1,2-minor determinant contribution in above 
  ! symmetrization matrix
  logical :: tsymdfq0dn
  ! complex RPA macroscopic dielectric function
  real(8), allocatable :: mdfrpa(:,:,:)
  ! derivative of real part of RPA macroscopic dielectric function
  real(8), allocatable :: mdfrpad(:,:)
  ! true if a Lorentzian broadening is to be used in the default linear optics
  logical :: lorentz
  ! sampling type for Brillouin zone (0 Lorentzian broadening, 1 tetrahedron
  ! method)
  integer :: bzsampl

  !----------------------------------!
  !     angular momenta variables    !
  !----------------------------------!
  ! maximum of maximum angular momenta
  integer :: lmaxmax
  ! (lmaxmax+1)^2
  integer :: lmmaxmax

  !----------------------------!
  !     xc-kernel variables    !
  !----------------------------!
  ! muffin-tin real space exchange-correlation kernel
  complex(8), allocatable :: fxcmt(:,:,:)
  ! interstitial real space exchange-correlation kernel
  complex(8), allocatable :: fxcir(:)
  ! Fourier transform of real space exchange-correlation kernel
  complex(8), allocatable :: fxcft(:)
  ! head element of real exchange-correlation kernel
  real(8), allocatable :: fxc0(:,:)
  ! derivative of head element of real exchange-correlation kernel
  real(8), allocatable :: fxc0d(:,:)
  ! exchange-correlation kernel functional type
  integer :: fxctype
  ! exchange-correlation kernel functional description
  character(256) fxcdescr
  ! exchange-correlation kernel functional spin treatment
  integer :: fxcspin
  ! number of G-vectors for ALDA kernel (twice the G cutoff)
  integer :: ngveca
  ! alpha-parameter for the asymptotic long range part of the kernel
  ! (see [Reining PRL 2002])
  real(8) :: alphalrc
  ! alpha-parameter for the asymptotic long range part of the kernel
  ! (see [Botti PRB 2005])
  real(8) :: alphalrcdyn
  ! beta-parameter for the asymptotic long range part of the kernel
  ! (see [Botti PRB 2005])
  real(8) :: betalrcdyn

  !---------------------------!
  !     exciton variables     !
  !---------------------------!
  ! maximum number of excitons
  integer :: nexcitmax
  ! number of excitons
  integer :: nexcit(3)
  ! exciton energies
  real(8), allocatable :: excite(:,:)
  ! exciton oscillator strengths
  real(8), allocatable :: excito(:,:)

  !-----------------------!
  !     I/O variables     !
  !-----------------------!
!!$  ! file extension (to be modified to characterize/enumerate outputs)
!!$  character(256) :: tdfilext
  ! file name for resume file
  character(256) :: fnresume
  ! last value of filext
  character(256) :: filextrevert
  ! file unit for output
  integer :: unitout
  ! file units to be connected at the same time
  integer :: unit1, unit2, unit3, unit4, unit5, unit6, unit7, unit8, unit9
  ! name for status file
  character(256) :: statusnam
  ! filename for output
  character(256) :: tdfileout
  ! eigenvectors contracted with matchin coefficients
  character(256) :: fnevapw
  ! weights for Brillouin zone integration
  character(256) :: fnwtet, fnwtet_t
  ! momentum matrix elements
  character(256) :: fnpmat, fnpmat_t
  ! exponential factor matrix elements
  character(256) :: fnemat, fnemat_t
  ! exponential factor matrix elements timing
  character(256) :: fnetim
  ! Kohn-Sham response function timing
  character(256) :: fnxtim
  ! Kohn-Sham energy differences
  character(256) :: fndevalsv, fndevalsv_t
  ! Kohn-Sham response function
  character(256) :: fnchi0, fnchi0_t, fnchi0p
  ! Inverse dielectric function
  character(256) :: fnieps, fnieps_t
  ! exciton file
  character(256) :: fnexciton
  ! macroscopic dielectric function
  character(256) :: fneps
  ! loss function
  character(256) :: fnloss
  ! optical conductivity
  character(256) :: fnsigma
  ! sumrules for optics
  character(256) :: fnsumrules

  !------------------------------!
  !     parallel environment     !
  !------------------------------!
  ! maximum number of processors allowed to use
  integer, parameter :: maxproc=1000
  ! current initial q-point index
  integer :: qpari
  ! current final q-point index
  integer :: qparf
  ! current initial k-point index
  integer :: kpari
  ! current final k-point index
  integer :: kparf
  ! current initial w-point index
  integer :: wpari
  ! current final w-point index
  integer :: wparf

  !--------------------------!
  !     Timing variables     !
  !--------------------------!
  ! initial and final timings for wall clock
  integer :: systim0i, systim0f, cntrate, systimcum
  ! initial and final timings for CPU timing
  real(8) :: cputim0i, cputim0f, cputimcum
  ! muffin-tin timings
  real(8) :: cmt0,cmt1,cmt2,cmt3,cmt4
  real(8) :: cpumtaa,cpumtalo,cpumtloa,cpumtlolo

  !-----------------------------!
  !     numerical constants     !
  !-----------------------------!
  ! array of (-i)**l values
  complex(8), allocatable :: zmil(:)
  ! conversion from hartree to electron volt
  real(8), parameter :: h2ev = 27.2114

  !---------------------------------!
  !     miscellaneous variables     !
  !---------------------------------!
  ! build number
  integer, parameter :: build_tddft = 130
  ! true if energies output in eV
  logical :: tevout
  ! scaling factor for writing energies
  real(8) :: escale
  ! debugging level
  integer :: dbglev
  ! true if to append info to output file
  logical :: tappinfo
  ! task for TDDFT part
  integer :: tasktd
  ! gather option
  logical :: gather
  data gather /.false./
  ! string for messages
  character(1024) :: msg
  ! default file extension
  data msg / 'no message' /  
  ! number of times the main tddft routine was called
  integer :: calledtd
  data calledtd / 0 /
  ! true if to skip allocations of radial functions in "init1"
  logical :: skipallocs1
  data skipallocs1 /.false./
  ! true if specified vertices are included in k-path and
  ! additionally bandstructure that is not shifted to the Fermi level
  logical :: imbandstr
  data imbandstr /.false./
  ! analytic evaluation of momentum matrix elements in interstitial
  logical :: pmatira
  ! true if state is only allowed to be read from STATE.OUT file
  ! and from no other file extension
  logical :: isreadstate0
  data isreadstate0 /.false./
  ! verbose output within SCF loop
  logical :: verbscf

end module modtddft
