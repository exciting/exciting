
! Copyright (C) 2004-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
#include "maxdefinitions.inc"
Module modxs
! !DESCRIPTION:
!   Global variables for the {\tt XS} (eXcited States) implementation
!   within the {\tt EXCITING}-code.
!
! !REVISION HISTORY:
!
!  Created June 2004 (Sagmeister)
!  Modifications due to XML input, 2008-2010 (Sagmeister)
      Implicit None

      Type mtints_type
        Complex (8), Allocatable :: aa (:, :, :)
        Complex (8), Allocatable :: alo (:, :, :)
        Complex (8), Allocatable :: loa (:, :, :)
        Complex (8), Allocatable :: lolo (:, :, :)        
      end type
      Type fftmap_type
        integer, pointer :: igfft(:)
        integer :: ngrid(3),ngrtot
        
      end type

  !----------------------------!
  !     symmetry variables     !
  !----------------------------!
  ! maximum allowed number of symmetry operations (private to this module)
      Integer, Private, Parameter :: maxsymcrs = 192
  ! map to inverse crystal symmetry
      Integer :: scimap (maxsymcrs)

  !------------------------------!
  !     q-point set variables    !
  !------------------------------!
  ! total number of q-points (reduced set)
      Integer :: nqptr
  ! locations of q-points on integer grid (reduced set)
      Integer, Allocatable :: ivqr (:, :)
  ! map from non-reduced grid to reduced set (reduced set)
      Integer, Allocatable :: iqmapr (:, :, :)
  ! q-points in lattice coordinates (reduced set)
      Real (8), Allocatable :: vqlr (:, :)
  ! q-points in Cartesian coordinates (reduced set)
      Real (8), Allocatable :: vqcr (:, :)
  ! q-point weights (reduced set)
      Real (8), Allocatable :: wqptr (:)
  ! finite momentum transfer q-vector
      Real (8), Allocatable :: vqlmt (:, :)
  ! finite momentum transfer G-vector
      Integer, Allocatable :: ivgmt (:, :)
  ! index of current q-point
      Integer :: iqcu
      Data iqcu / 0 /
  ! number of crystal symmetries for the little group of q
      Integer, Allocatable :: nsymcrysq (:)
  ! map from little group of q to spacegroup
      Integer, Allocatable :: scqmap (:, :)
  ! wrapping vectors for elements of the small group of q
      Integer, Allocatable :: ivscwrapq (:, :, :)
!
  !----------------------------------!
  !     G+q-vector set variables     !
  !----------------------------------!
  ! true if |G| cutoff is used in place of |G+q| cutoff for momentum transfer calculations
      logical :: tgqmaxg
  ! G-vector grid sizes of (G+q)-vectors
      Integer :: ngridgq (3)
  ! integer grid intervals for each direction for G-vectors
      Integer :: intgqv (3, 2)
  ! number of G+q-vectors
      Integer, Allocatable :: ngq (:)
  ! maximum number of G+q-vectors over all q-points
      Integer :: ngqmax
  ! maximum number of augmentation functions
      Integer :: apwmaxsize
  ! maximum number of local orbitals 
      Integer :: lomaxsize
  ! number of local orbitals 
      Integer, allocatable :: losize(:),apwsize(:)
  ! index from G+q-vectors to G-vectors
      Integer, Allocatable :: igqig (:, :)
  ! map from integer grid to G+q-vector array
      Integer, Allocatable :: ivgigq (:, :, :, :)
  ! G+q-vectors in lattice coordinates
      Real (8), Allocatable :: vgql (:, :, :)
  ! G+q-vectors in Cartesian coordinates
      Real (8), Allocatable :: vgqc (:, :, :)
  ! length of G+q-vectors
      Real (8), Allocatable :: gqc (:, :)
  ! (theta, phi) coordinates of G+q-vectors
      Real (8), Allocatable :: tpgqc (:, :, :)
  ! structure factor for the G+q-vectors
      Complex (8), Allocatable :: sfacgq (:, :, :)
  ! spherical harmonics of the G-vectors
      Complex (8), Allocatable :: ylmgq (:, :, :)

  !---------------------------------!
  !     k-point set  variables      !
  !---------------------------------!
  ! number of k-points for q=0
      Integer :: nkpt0
  ! k-points in lattice coordinates for q=0
      Real (8), Allocatable :: vkl0 (:, :)
  ! maximum number of space group operations in stars over all k
      Integer :: nsymcrysstrmax
  ! number of space group operations for stars
      Integer, Allocatable :: nsymcrysstr (:)
  ! star of space group operations for k-points
      Integer, Allocatable :: scmapstr (:, :)
  ! star of k-point indices of non-reduced set
      Integer, Allocatable :: ikstrmapiknr (:, :)
  ! map from non-reduced k-point set to reduced one
      Integer, Allocatable :: strmap (:)
  ! map from non-reduced k-point set to associated symmetry in star
      Integer, Allocatable :: strmapsymc (:)

  !-------------------------!
  !     k+q-point set       !
  !-------------------------!
  ! offset for k+q-point set derived from q-point
      Real (8), Allocatable :: qvkloff (:, :)
  ! map from k-point index to k+q point index for same k
      Integer, Allocatable :: ikmapikq (:, :)

  !-----------------------------------------!
  !     G+k-vector set  variables (q=0)     !
  !-----------------------------------------!
  ! number of G+k-vectors for augmented plane waves
      Integer, Allocatable :: ngk0 (:, :)
  ! maximum number of G+k-vectors over all k-points
      Integer :: ngkmax0
  ! index from G+k-vectors to G-vectors
      Integer, Allocatable :: igkig0 (:, :, :)
  ! G+k-vectors in lattice coordinates
      Real (8), Allocatable :: vgkl0 (:, :, :, :)
  ! G+k-vectors in Cartesian coordinates
      Real (8), Allocatable :: vgkc0 (:, :, :, :)
  ! length of G+k-vectors
      Real (8), Allocatable :: gkc0 (:, :, :)
  ! (theta, phi) coordinates of G+k-vectors
      Real (8), Allocatable :: tpgkc0 (:, :, :, :)
  ! structure factor for the G+k-vectors
      Complex (8), Allocatable :: sfacgk0 (:, :, :, :)

  !-----------------------------------------!
  !     potential and density variables     !
  !-----------------------------------------!
  ! square root of Coulomb potential in G-space
      Real (8), Allocatable :: sptclg (:, :)
!
  !---------------------------------------!
  !     Hamiltonian and APW variables     !
  !---------------------------------------!
  ! maximum nmat over all k-points (q=0)
      Integer :: nmatmax0
  ! order of overlap and Hamiltonian matrices for each k-point (q=0)
      Integer, Allocatable :: nmat0 (:, :)
  ! first-variational eigenvectors
      Complex (8), Allocatable :: evecfv (:, :, :)
  ! second variational eigenvectors
      Complex (8), Allocatable :: evecsv (:, :)
  ! first-variational eigenvectors (q=0)
      Complex (8), Allocatable :: evecfv0 (:, :, :)
  ! first variational eigenvalues
      Real (8), Allocatable :: evalfv (:, :)
  ! second-variational eigenvalues
      Real (8), Allocatable :: evalsv0 (:, :)
  ! expansion coefficients of APW functions
      Complex (8), Allocatable :: apwcmt (:, :, :, :)
  ! expansion coefficients of APW functions (q=0)
      Complex (8), Allocatable :: apwcmt0 (:, :, :, :)

  ! expansion coefficients of local orbitals functions
      Complex (8), Allocatable :: locmt (:, :, :, :)
  ! expansion coefficients of local orbitals functions (q=0)
      Complex (8), Allocatable :: locmt0 (:, :, :, :)
      Complex (8), Allocatable :: cmtfun (:, :, :)
      Complex (8), Allocatable :: cmtfun0 (:, :, :)

  !--------------------------------------------!
  !     eigenvalue and occupancy variables     !
  !--------------------------------------------!
  ! eigenvalue differences (resonant part)
      Real (8), Allocatable :: deou (:, :)
  ! eigenvalue differences (anti-resonant part)
      Real (8), Allocatable :: deuo (:, :)
  ! occupation numbers (q=0)
      Real (8), Allocatable :: occsv0 (:, :)
  ! occupation number differences (first band combination)
      Real (8), Allocatable :: docc12 (:, :)
  ! occupation number differences (second band combination)
      Real (8), Allocatable :: docc21 (:, :)
  ! highest (at least partially) occupied state
      Integer, Allocatable :: isto0 (:), isto (:)
  ! lowest (at least partially) unoccupied state
      Integer, Allocatable :: istu0 (:), istu (:)
  ! maximum isto over k-points
      Integer :: istocc0, istocc
  ! minimum istu over k-points
      Integer :: istunocc0, istunocc
  ! number of (at least partially) occupied valence states
      Integer :: nstocc0, nstocc
  ! number of (at least partially) unoccupied valence states
      Integer :: nstunocc0, nstunocc
  ! highest (at least partially) occupied state energy
      Real (8) :: evlhpo
  ! lowest (at least partially) unoccupied state energy
      Real (8) :: evllpu
  ! lower and upper limits and numbers for band indices combinations
      Integer :: nst1, istl1, istu1, nst2, istl2, istu2
  ! lower and upper limits and numbers for band indices combinations, 2nd block
      Integer :: nst3, istl3, istu3, nst4, istl4, istu4
  ! minimum and maximum energies over k-points
      Real (8) :: evlmin, evlmax, evlmincut, evlmaxcut
  ! true if system has a Kohn-Sham gap
      Logical :: ksgap

  !--------------------------------------------------!
  !     matrix elements of exponential expression    !
  !--------------------------------------------------!
  ! (lmaxemat+1)^2
      Integer :: lmmaxemat
  ! (lmaxapwwf+1)^2
      Integer :: lmmaxapwwf
  ! Gaunt coefficients array
      Real (8), Allocatable :: xsgnt (:, :, :)
      Real (8), Allocatable :: xsgntou (:, :, :)
      Real (8), Allocatable :: xsgntuo (:, :, :)
      Real (8), Allocatable :: xsgntoo (:, :, :)
  ! radial integrals coefficients (APW-APW)
!      Complex (8), Allocatable :: intrgaa2 (:, :, :)
!      Complex (8), Allocatable :: intrgalo2 (:, :, :)
!      Complex (8), Allocatable :: intrgloa2 (:, :, :)
!      Complex (8), Allocatable :: intrglolo2 (:, :, :)
 ! radial integrals coefficients (lo-APW)
  ! radial integrals coefficients (APW-lo)
  ! radial integrals coefficients (lo-lo)
  ! radial integrals (APW-APW)
      Real (8), Allocatable :: riaa (:, :, :, :, :, :, :)
  ! radial integrals (lo-APW)
      Real (8), Allocatable :: riloa (:, :, :, :, :, :)
  ! radial integrals (lo-lo)
      Real (8), Allocatable :: rilolo (:, :, :, :, :)
  ! helper matrix
      Complex (8), Allocatable :: xih (:, :)
  ! helper matrix
!      Complex (8), Allocatable :: xihir (:, :)
  ! helper matrix
      Complex (8), Allocatable :: xiohalo (:, :)
  ! helper matrix
      Complex (8), Allocatable :: xiuhloa (:, :)
  ! matrix elements array (resonant part)
      Complex (8), Allocatable :: xiou (:, :, :)
  ! matrix elements array (anti-resonant part)
      Complex (8), Allocatable :: xiuo (:, :, :)

  !---------------------------------!
  !     momentum matrix elements    !
  !---------------------------------!
  ! radial integrals coefficients (APW-APW)
      Real (8), Allocatable :: ripaa (:, :, :, :, :, :)
  ! radial integrals coefficients (APW-lo)
      Real (8), Allocatable :: ripalo (:, :, :, :, :, :)
  ! radial integrals coefficients (lo-APW)
      Real (8), Allocatable :: riploa (:, :, :, :, :, :)
  ! radial integrals coefficients (lo-lo)
      Real (8), Allocatable :: riplolo (:, :, :, :, :, :)
  ! momentum matrix elements (resonant part)
      Complex (8), Allocatable :: pmou (:, :, :)
  ! momentum matrix elements (anti-resonant part)
      Complex (8), Allocatable :: pmuo (:, :, :)

  !------------------------------------------!
  !     response and dielectric functions    !
  !------------------------------------------!
  ! factor for time-ordering
      Real (8) :: tordf
  ! number of energy intervals
      Integer :: nwdf
  ! symmetrization tensor
      Real (8) :: symt2 (3, 3, 3, 3)
  ! sampling type for Brillouin zone (0 Lorentzian broadening, 1 tetrahedron
  ! method)
      Integer :: bzsampl

  !----------------------------!
  !     xc-kernel variables    !
  !----------------------------!
  ! factor for time-ordering
      Real (8) :: torfxc
  ! muffin-tin real space exchange-correlation kernel
      Complex (8), Allocatable :: fxcmt (:, :, :)
  ! interstitial real space exchange-correlation kernel
      Complex (8), Allocatable :: fxcir (:)
  ! exchange-correlation kernel functional description
      Character (256) :: fxcdescr
  ! exchange-correlation kernel functional spin treatment
      Integer :: fxcspin

  !---------------------------!
  !     exciton variables     !
  !---------------------------!
  ! number of excitons
      Integer :: nexcit (3)
  ! exciton energies
      Real (8), Allocatable :: excite (:, :)
  ! exciton oscillator strengths
      Real (8), Allocatable :: excito (:, :)

  !-----------------------------!
  !     screening variables     !
  !-----------------------------!
  ! true if one of the screening tasks is executed
      Logical :: tscreen
  ! dielectric tensor in the RPA
      Complex (8) :: dielten (3, 3)
  ! dielectric tensor in the independent particle approximation
      Complex (8) :: dielten0 (3, 3)
  ! (lmaxdielt+1)^2
      Integer :: lmmaxdielt
  ! true if Lebedev Laikov meshes are to be used
      Logical :: tleblaik

  !------------------------------------------!
  !     Bethe-Salpeter (kernel) variables    !
  !------------------------------------------!
      Logical :: tfxcbse
  ! number of states below Fermi energy (Coulomb - and exchange term)
      Integer :: nbfce
  ! number of states above Fermi energy (Coulomb - and exchange term)
      Integer :: nafce
  ! number of states below Fermi energy
      Integer :: nbfbse
  ! number of states above Fermi energy
      Integer :: nafbse
  ! diagonal of BSE kernel (mean value, lower, upper limit and range)
      Complex (8) :: bsed, bsedl, bsedu, bsedd
  ! BSE matrix sizes
	  Integer :: sta1, sto1, sta2, sto2
  !-----------------------!
  !     I/O variables     !
  !-----------------------!
  ! file name for resume file
      Character (256) :: fnresume
  ! last value of filext
      Character (256) :: filextrevert
  ! file unit for output
      Integer :: unitout
  ! file units to be connected at the same time
      Integer :: unit1, unit2, unit3, unit4, unit5, unit6, unit7, &
     & unit8, unit9
  ! filename for output
      Character (256) :: xsfileout
  ! weights for Brillouin zone integration
      Character (256) :: fnwtet
  ! momentum matrix elements
      Character (256) :: fnpmat, fnpmat_t
  ! exponential factor matrix elements
      Character (256) :: fnemat, fnemat_t
  ! exponential factor matrix elements timing
      Character (256) :: fnetim
  ! Kohn-Sham response function timing
      Character (256) :: fnxtim
  ! Kohn-Sham response function
      Character (256) :: fnchi0, fnchi0_t
  ! macroscopic dielectric function
      Character (256) :: fneps
  ! loss function
      Character (256) :: fnloss
  ! optical conductivity
      Character (256) :: fnsigma
  ! MOKE spectrum
      Character (256) :: fnmoke
  ! sumrules for optics
      Character (256) :: fnsumrules

  !--------------------------!
  !     backup variables     !
  !--------------------------!
  ! filename extension
      Character (256) :: filext_b
  ! nosym is .true. if no symmetry information should be used
      Logical :: nosym_b
  ! smallest muffin-tin radius times gkmax
      Real (8) :: rgkmax_b
  ! number of empty states
      Integer :: nempty_b
  ! reducek is .true. if k-points are to be reduced (with crystal symmetries)
      Logical :: reducek_b
  ! k-point grid sizes
      Integer :: ngridk_b (3)
  ! k-point offset
      Real (8) :: vkloff_b (3)
  ! q-point grid sizes
      Integer :: ngridq_b (3)
  ! reducek is .true. if q-points are to be reduced (with crystal symmetries)
      Logical :: reduceq_b
  ! type of matrix element generation (band-combinations)
      Integer :: emattype_b
      Real (8) :: swidth_b
      Integer :: lmaxapw_b
      Integer :: lmaxmat_b
  ! maximum self-consistency steps
      Integer :: maxscl_b
  ! bfieldc
      Real (8) :: bfieldc_b (3)

  !------------------------------!
  !     parallel environment     !
  !------------------------------!
  ! maximum number of processors allowed to use
      Integer, Parameter :: maxproc = 1000
  ! parallelization type (values are 'q', 'k', 'w')
      Character (1) :: partype
  ! current initial q-point index
      Integer :: qpari
  ! current final q-point index
      Integer :: qparf
  ! current initial k-point index
      Integer :: kpari
  ! current final k-point index
      Integer :: kparf
  ! current initial (k,kp) pair index
      Integer :: ppari
  ! current final (k,kp) pair index
      Integer :: pparf
  ! current initial w-point index
      Integer :: wpari
  ! current final w-point index
      Integer :: wparf

  !--------------------------!
  !     Timing variables     !
  !--------------------------!
  ! initial and final timings for wall clock
      Integer :: systim0i, systim0f, cntrate, systimcum
  ! initial and final timings for CPU timing
      Real (8) :: cputim0i, cputim0f
  ! muffin-tin timings
!      Real (8) :: cmt0, cmt1, cmt2, cmt3, cmt4
      Real (8) :: cpumtaa, cpumtalo, cpumtloa, cpumtlolo

  !---------------------------------!
  !     miscellaneous variables     !
  !---------------------------------!
  ! spherical covering set in Cartesian coordinates
      Real (8), Allocatable :: sphcov (:, :)
  ! spherical covering set in tetha/phi angles
      Real (8), Allocatable :: sphcovtp (:, :)
  ! scaling factor for writing energies
      Real (8) :: escale
  ! string for messages
      Character (1024) :: msg
      Data msg / 'no message' /
  ! number of times the main excited states routine was called
      Integer :: calledxs
      Data calledxs / 0 /
  ! true if only symmetries are recalculated in "init0"
      Logical :: init0symonly
      Data init0symonly / .False. /
  ! true if to skip allocations of radial functions in "init1"
      Logical :: init1norealloc
      Data init1norealloc / .False. /
  ! true if state (density and potential) is only allowed to be read from
  ! STATE.OUT file (no other file extension allowed)
      Logical :: isreadstate0
      Data isreadstate0 / .False. /
  ! true if calculation of matrix elements of exponential expression
  ! is needed
      Logical :: temat
      Data temat / .True. /
  ! set true if HF-Hybrids are used as starting point
      Logical :: hybridhf
      Data  hybridhf / .false. /
  ! if set to true KS eingenvalues and eigenvectors are not recalculated
      Logical :: skipgnd
!STK: include variables necessary for double grid computations

  !---------------------------------------------!
  !     sub-grid k-vectors for double grids     !
  !---------------------------------------------!
  ! number of sub k-points
      Integer :: nksubpt
  ! sub kpt that is worked on
      Integer :: iksubpt
  ! sub k-points in lattice coordinates
      Real (8), Allocatable :: vksubl (:, :)
  ! sub k-points in cartesian coordinates
      Real (8), Allocatable :: vksubc (:, :)
  ! k-point weights
      Real (8), Allocatable :: wksubpt (:)
  ! locations of k-points on integer grid
      Integer, Allocatable :: ivksub (:, :)
  ! map from non-reduced grid to reduced set
      Integer, Allocatable :: iksubmap (:, :, :)

  !---------------------------------------!
  !     other double grid variables       !
  !---------------------------------------!
  ! are we doing a double grid run?
      Logical :: dgrid
  ! backup input%xs%screening%do 
      Character (11) :: doscreen0
  ! backup for XS vkloff
      Real (8) :: vkloff_xs_b(3)


End Module modxs
