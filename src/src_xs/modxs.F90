!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module modxs
! !DESCRIPTION:
!   Global variables for the {\tt XS} (eXcited States) implementation
!   in the {\tt EXCITING}-code.
!
! !REVISION HISTORY:
!
!  Created June 2004 (Sagmeister)
      Implicit None
!
  !----------------------------!
  !     symmetry variables     !
  !----------------------------!
  ! maximum allowed number of symmetry operations (private to this module)
      Integer, Private, Parameter :: maxsymcrs = 192
  ! true if only symmorphic space-group operations are to be considered
  ! allow only symmetries without non-primitive translations
!replaced by inputstructure  logical :: symmorph
  ! map to inverse crystal symmetry
      Integer :: scimap (maxsymcrs)
!
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
  ! number of Q-points for momentum transfer
  !integer ::size(input%xs%qpointset%qpoint,2)
  ! finite momentum transfer G+q-vector
 ! real(8), allocatable :: input%xs%qpointset%qpoint(:,:)
  ! finite momentum transfer q-vector
      Real (8), Allocatable :: vqlmt (:, :)
  ! finite momentum transfer G-vector
      Integer, Allocatable :: ivgmt (:, :)
  ! treatment of macroscopic dielectric function for Q-point outside of
  ! Brillouin zone
!replaced by inputstructure  integer :: mdfqtype
  ! index of current q-point
      Integer :: iqcu
      Data iqcu / 0 /
!
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
  ! G-vector grid sizes of (G+q)-vectors
      Integer :: ngridgq (3)
  ! integer grid intervals for each direction for G-vectors
      Integer :: intgqv (3, 2)
  ! maximum |G+q| cut-off for APW functions
!replaced by inputstructure  real(8)::gqmax
  ! number of G+q-vectors
      Integer, Allocatable :: ngq (:)
  ! maximum number of G+q-vectors over all q-points
      Integer :: ngqmax
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
!
  !---------------------------------!
  !     k-point set  variables      !
  !---------------------------------!
  ! number of k-points for q=0
      Integer :: nkpt0
  ! k-points in lattice coordinates for q=0
      Real (8), Allocatable :: vkl0 (:, :)
!
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
!
!
  !-------------------------!
  !     k+q-point set       !
  !-------------------------!
  ! offset for k+q-point set derived from q-point
      Real (8), Allocatable :: qvkloff (:, :)
  ! map from k-point index to k+q point index for same k
      Integer, Allocatable :: ikmapikq (:, :)
!
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
!
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
!
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
!
  !--------------------------------------------------!
  !     matrix elements of exponential expression    !
  !--------------------------------------------------!
  ! fast method to calculate APW-lo, lo-APW and lo-lo parts in MT
!replaced by inputstructure  logical :: fastemat
  ! type of matrix element generation (band-combinations)
!replaced by inputstructure  integer :: emattype
  ! maximum angular momentum for Rayleigh expansion of exponential
!replaced by inputstructure  integer :: lmaxemat
  ! (lmaxemat+1)^2
      Integer :: lmmaxemat
  ! maximum angular momentum for APW functions (for matrix elements)
!replaced by inputstructure  integer :: lmaxapwwf
  ! (lmaxapwwf+1)^2
      Integer :: lmmaxapwwf
  ! Gaunt coefficients array
      Real (8), Allocatable :: xsgnt (:, :, :)
  ! radial integrals coefficients (APW-APW)
      Complex (8), Allocatable :: intrgaa (:, :, :, :, :)
  ! radial integrals coefficients (lo-APW)
      Complex (8), Allocatable :: intrgloa (:, :, :, :, :)
  ! radial integrals coefficients (APW-lo)
      Complex (8), Allocatable :: intrgalo (:, :, :, :, :)
  ! radial integrals coefficients (lo-lo)
      Complex (8), Allocatable :: intrglolo (:, :, :, :, :)
  ! radial integrals (APW-APW)
      Real (8), Allocatable :: riaa (:, :, :, :, :, :, :)
  ! radial integrals (lo-APW)
      Real (8), Allocatable :: riloa (:, :, :, :, :, :)
  ! radial integrals (lo-lo)
      Real (8), Allocatable :: rilolo (:, :, :, :, :)
  ! helper matrix
      Complex (8), Allocatable :: xih (:, :)
  ! helper matrix
      Complex (8), Allocatable :: xihir (:, :)
  ! helper matrix
      Complex (8), Allocatable :: xiohalo (:, :)
  ! helper matrix
      Complex (8), Allocatable :: xiuhloa (:, :)
  ! matrix elements array (resonant part)
      Complex (8), Allocatable :: xiou (:, :, :)
  ! matrix elements array (anti-resonant part)
      Complex (8), Allocatable :: xiuo (:, :, :)
!
  !---------------------------------!
  !     momentum matrix elements    !
  !---------------------------------!
  ! fast method to calculate matrix elements
!replaced by inputstructure  logical :: fastpmat
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
!
  !------------------------------------------!
  !     response and dielectric functions    !
  !------------------------------------------!
  ! time ordering of response function (time-ordered/retarded)
!replaced by inputstructure  character(32) :: torddf
  ! factor for time-ordering
      Real (8) :: tordf
  ! true if analytic continuation to the real axis is to be performed
!replaced by inputstructure  logical :: acont
  ! number of energy intervals
      Integer :: nwdf
  ! number of energy intervals (on imaginary axis) for analytic continuation
!replaced by inputstructure  integer :: nwacont
  ! broadening for Kohn Sham response function
!replaced by inputstructure  real(8) :: broad
  ! true if Lindhard like function is calculated (trivial matrix elements)
!replaced by inputstructure  logical :: lindhard
  ! true if to consider the anti-resonant part for the dielectric function
!replaced by inputstructure  logical :: aresdf
  ! true if only diagonal part of xc-kernel is used
!replaced by inputstructure  logical :: kerndiag
  ! true if off-diagonal tensor components of dielectric function are calculated
!replaced by inputstructure  logical :: dfoffdiag
  ! symmetrization tensor
      Real (8) :: symt2 (3, 3, 3, 3)
  ! true if tetrahedron method is used for dielectric function/matrix
!replaced by inputstructure  logical :: tetradf
  ! sampling type for Brillouin zone (0 Lorentzian broadening, 1 tetrahedron
  ! method)
      Integer :: bzsampl
  ! choice of weights and nodes for tetrahedron method and non-zero Q-point
!replaced by inputstructure  integer :: tetraqweights
  ! number of band transitions for analysis
      Integer :: ndftrans
  ! k-point and band combination analysis
      Integer, Allocatable :: dftrans (:, :)
  ! smallest energy difference for which the inverse square will be considered
!replaced by inputstructure  real(8) :: epsdfde
  ! cutoff energy for dielectric function
!replaced by inputstructure  real(8) :: emaxdf
!
  !----------------------------!
  !     xc-kernel variables    !
  !----------------------------!
  ! time ordering of xc-kernel function (time-ordered/retarded)
!replaced by inputstructure  character(32) :: tordfxc
  ! factor for time-ordering
      Real (8) :: torfxc
  ! true if to consider the anti-resonant part
!replaced by inputstructure  logical :: aresfxc
  ! maximum angular momentum for Rayleigh expansion of exponential in
  ! ALDA-kernel
!replaced by inputstructure  integer :: lmaxalda
  ! muffin-tin real space exchange-correlation kernel
      Complex (8), Allocatable :: fxcmt (:, :, :)
  ! interstitial real space exchange-correlation kernel
      Complex (8), Allocatable :: fxcir (:)
  ! exchange-correlation kernel functional type
!replaced by inputstructure  integer :: fxctype
  ! exchange-correlation kernel functional description
      Character (256) :: fxcdescr
  ! exchange-correlation kernel functional spin treatment
      Integer :: fxcspin
  ! alpha-parameter for the asymptotic long range part of the kernel
  ! (see [Reining PRL 2002])
!replaced by inputstructure  real(8) :: alphalrc
  ! alpha-parameter for the asymptotic long range part of the kernel
  ! (see [Botti PRB 2005])
!replaced by inputstructure  real(8) :: alphalrcdyn
  ! beta-parameter for the asymptotic long range part of the kernel
  ! (see [Botti PRB 2005])
!replaced by inputstructure  real(8) :: betalrcdyn
  ! split parameter for degeneracy in energy differences of BSE-kernel
!replaced by inputstructure  real(8) :: fxcbsesplit
!
  !---------------------------!
  !     exciton variables     !
  !---------------------------!
  ! maximum number of excitons
!replaced by inputstructure  integer :: nexcitmax
  ! number of excitons
      Integer :: nexcit (3)
  ! exciton energies
      Real (8), Allocatable :: excite (:, :)
  ! exciton oscillator strengths
      Real (8), Allocatable :: excito (:, :)
!
  !-----------------------------!
  !     screening variables     !
  !-----------------------------!
  ! true if one of the screening tasks is executed
      Logical :: tscreen
  ! true if q-point set is taken from first Brillouin zone
!replaced by inputstructure  logical :: fbzq
  ! screening type: can be either "full", "diag", "noinvdiag" or "constant"
!replaced by inputstructure  character(32) :: screentype
  ! nosym is .true. if no symmetry information should be used
!replaced by inputstructure  logical::nosymscr
  ! reducek is .true. if k-points are to be reduced (with crystal symmetries)
!replaced by inputstructure  logical::reducekscr
  ! k-point grid sizes
!replaced by inputstructure  integer :: ngridkscr(3)
  ! k-point offset
!replaced by inputstructure  real(8) :: vkloffscr(3)
  ! smallest muffin-tin radius times gkmax
!replaced by inputstructure  real(8) :: rgkmaxscr
  ! number of empty states
!replaced by inputstructure  integer :: nemptyscr
  ! Hermitian treatment
!replaced by inputstructure  integer :: scrherm
  ! dielectric tensor in the RPA
      Complex (8) :: dielten (3, 3)
  ! dielectric tensor in the independent particle approximation
      Complex (8) :: dielten0 (3, 3)
  ! averaging type for singular term in screenend Coulomb interaction
      Character (256) :: sciavtype
  ! average of body for screened Coulomb interaction at Gamma-point
!replaced by inputstructure  logical :: sciavbd
  ! average of head, wings and body for screened Coulomb interaction at
  ! non-zero q-point
!replaced by inputstructure!replaced by inputstructure!replaced by inputstructure  logical :: sciavqhd, sciavqwg, sciavqbd
  ! maximum angular momentum for angular average of dielectric tensor
!replaced by inputstructure  integer :: lmaxdielt
  ! (lmaxdielt+1)^2
      Integer :: lmmaxdielt
  ! number of points for Lebedev Laikov meshes
!replaced by inputstructure  integer :: nleblaik
  ! true if Lebedev Laikov meshes are to be used
      Logical :: tleblaik
!
  !------------------------------------------!
  !     Bethe-Salpeter (kernel) variables    !
  !------------------------------------------!
  ! type of BSE-Hamiltonian
!replaced by inputstructure  character(32) :: bsetype
  ! true if effective singular part of direct term of BSE Hamiltonian is to be used
!replaced by inputstructure  logical :: bsedirsing
  ! nosym is .true. if no symmetry information should be used
!replaced by inputstructure  logical::nosymbse
  ! reducek is .true. if k-points are to be reduced (with crystal symmetries)
!replaced by inputstructure  logical::reducekbse
  ! k-point offset
!replaced by inputstructure  real(8) :: vkloffbse(3)
  ! smallest muffin-tin radius times gkmax
!replaced by inputstructure  real(8) :: rgkmaxbse
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
!
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
  ! sumrules for optics
      Character (256) :: fnsumrules
!
!
  !------------------------------------------!
  !     xs-parameters related to GS ones     !
  !------------------------------------------!
!replaced by inputstructure  logical :: nosymxs
!replaced by inputstructure  integer :: ngridkxs(3)
!replaced by inputstructure  real(8) :: vkloffxs(3)
!replaced by inputstructure  logical :: reducekxs
!replaced by inputstructure  integer :: ngridqxs(3)
!replaced by inputstructure  logical :: reduceqxs
!replaced by inputstructure  real(8) :: rgkmaxxs
!replaced by inputstructure  real(8) :: swidthxs
!replaced by inputstructure  integer :: lmaxapwxs
!replaced by inputstructure  integer :: lmaxmatxs
!replaced by inputstructure  integer :: nemptyxs
!
!
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
!
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
!
  !--------------------------!
  !     Timing variables     !
  !--------------------------!
  ! initial and final timings for wall clock
      Integer :: systim0i, systim0f, cntrate, systimcum
  ! initial and final timings for CPU timing
      Real (8) :: cputim0i, cputim0f, cputimcum
  ! muffin-tin timings
      Real (8) :: cmt0, cmt1, cmt2, cmt3, cmt4
      Real (8) :: cpumtaa, cpumtalo, cpumtloa, cpumtlolo
!
  !-----------------------------!
  !     numerical constants     !
  !-----------------------------!
  ! Kronecker delta
      Integer, Parameter :: krondelta (3, 3) = reshape ( (/ 1, 0, 0, 0, &
     & 1, 0, 0, 0, 1 /), (/ 3, 3 /))
  ! conversion from hartree to electron volt
      Real (8), Parameter :: h2ev = 27.2114d0
!
  !---------------------------------!
  !     miscellaneous variables     !
  !---------------------------------!
  ! spherical covering set in Cartesian coordinates
      Real (8), Allocatable :: sphcov (:, :)
  ! spherical covering set in tetha/phi angles
      Real (8), Allocatable :: sphcovtp (:, :)
  ! xs code version
      Integer :: versionxs (2)
  ! true if energies output in eV
!replaced by inputstructure  logical :: tevout
  ! scaling factor for writing energies
      Real (8) :: escale
  ! debugging level
!replaced by inputstructure  integer :: dbglev
  ! true if to append info to output file
!replaced by inputstructure  logical :: tappinfo
  ! gather option
!replaced by inputstructure  logical :: gather
!  data gather /.false./
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
!
End Module modxs
!
