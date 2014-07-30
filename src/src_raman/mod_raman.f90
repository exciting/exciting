! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
!
! parameters
module raman_params
! source: http://physics.nist.gov/constants
real(8), parameter :: pi      =      3.14159265359d0
real(8), parameter :: faua    =      0.529177211d0     ! 1 Bohr = faua Angstrom
real(8), parameter :: fau3cm3 =      1.d-24*faua**3    ! 1 Bohr^3 = fau3cm3 cm^3
real(8), parameter :: fryev   =     13.60569253d0      ! 1 Ry = fryev eV
real(8), parameter :: fhaev   =     27.21138505d0      ! 1 Ha = fhaev eV
real(8), parameter :: fevry   =      1.d0/fryev        ! 1 eV = fevry Ry
real(8), parameter :: fevha   =      1.d0/fhaev        ! 1 eV = fevha Ha
real(8), parameter :: frywn   =      1.09737315685d5   ! 1 Ry = frywn cm^-1
real(8), parameter :: fhawn   =      2.194746313708d5  ! 1 Ha = fhawn cm^-1
real(8), parameter :: fwnry   =      1.d0/frywn        ! 1 cm^-1 = fwnry Ry
real(8), parameter :: fwnha   =      1.d0/fhawn        ! 1 cm^-1 = fwnha Ha
real(8), parameter :: fkwn    =      0.69503476d0      ! 1 K = fkwn cm^-1
real(8), parameter :: fwnk    =      1.d0/fkwn         ! 1 cm^-1 = fwnk K
real(8), parameter :: fryk    =      frywn*fwnk        ! 1 Ry = fryk K
real(8), parameter :: fhak    =      fhawn*fwnk        ! 1 Ha = fhak K
real(8), parameter :: fkry    =      1.0d0/fryk        ! 1 K = fkry Ry
real(8), parameter :: fwnmev  =      0.1239841930d0    ! 1 cm^-1 = fwnmev meV
real(8), parameter :: fmevwn  =      1.d0/fwnmev       ! 1 meV = fmevwn cm^-1
real(8), parameter :: frnmha  =      1.5198298460045d-7*299792458   ! 1/nm = frnmha Ha
real(8), parameter :: fharnm  =      1.d0/frnmha       ! 1 Ha = fharnm 1/nm
real(8), parameter :: famuau  =      1822.88848        ! 1 amu = famuau atomic units of mass
! these are derived later on from the above
real(8) :: fener,flang
! speed of light in vacuum
real(8), parameter :: speed_of_light = 299792458.d0    ! m s^-1
! threshold for symmetry analysis
real(8), parameter :: eps     =      1.0d-5
! factorial
real(8), dimension(6), parameter :: factorial = (/ 1.d0, 2.d0, 6.d0, 24.d0, 120.d0, 720.d0 /)
! complex constants
complex(8), parameter :: zzero = (0.d0, 0.d0)
complex(8), parameter :: zone  = (1.d0, 0.d0)
complex(8), parameter :: ztwo  = (2.d0, 0.d0)
end module raman_params
!
!
! eigenvalues, transitions matrix elements, eigen functions
module raman_ew
use raman_params
real(8), allocatable :: eigen(:)
real(8), allocatable :: transme1(:),transme2(:),transme3(:),transme4(:),transme5(:),transme6(:)
real(8), allocatable :: z1(:,:),z2(:,:)
real(8), allocatable :: e1(:),e2(:),e3(:),de(:)
real(8) :: fgew
integer, allocatable :: indexi(:)
end module raman_ew
!
!
!
! matrices T and R
module raman_trmat
use raman_params
real(8), allocatable :: T(:,:),R(:,:)
end module raman_trmat
!
!
! coefficients for potential (a's) and derivatives for the diagonal elements of the dielectric tensor
! these are potentially changed from input values during execution 
module raman_coeff
use raman_params
! potential coefficients
real(8) :: A0,A1,A2,A3,A4,A5,A6
! number of unit cells in scattering volume
integer :: ncell
! potential raw data, energies
real(8), allocatable :: potinx(:, :),potiny(:, :)
! number of steps per mode
integer :: numpot
! dielec raw data
complex (8), allocatable :: df(:, :, :, :, :)
! dielec coefficients
complex (8), dimension (3, 3) :: deq, d2eq, d3eq, d4eq, d5eq, d6eq
! off-diagonal optical components
logical :: offdiag
! name of components
character(2), dimension(3, 3) :: comp = reshape ( (/'11', '21', '31', '12', '22', '32', '13', '23', '33'/), (/ 3, 3 /) )
end module raman_coeff
!
!
!
! finite element method matrices
module raman_fij
use raman_params
real(8) :: FI(4,4,0:6),DFI(4,4)
real(8) :: sfact
end module raman_fij
!
!
!
! interval
module raman_inter
use raman_params
! low value of every interval
real(8), allocatable :: xa(:)
! fine mesh for interval with maxp knots
real(8), allocatable :: xpot(:)
! potential values at maxp
real(8), allocatable :: pot(:)
! interval height
real(8) :: h
real(8), allocatable :: b0(:),b1(:),b2(:),b3(:),b4(:)
end module raman_inter
!
!
!
! Raman relevant symmetry variables
module raman_symmetry
use raman_params
! number of crystal SOPs, for equilibirum geometry
integer :: numsop
! rotational parts of crystal SOPs, real, for equilibirum geometry
real(8) :: sopmat(3, 3, 48)
! rotational parts of crystal SOPs, cartesian, for equilibirum geometry
real(8) :: sopmatc(3, 3, 48)
! character table
! note that due to the diagonalization the eigenvector representing IREP i is in column i, instead of in row i
complex (8),allocatable :: charact(:,:)
! number of classes
integer :: cl
! class of a give SOP
integer :: class(48)
! elements in a class
integer :: elem_cl(48)
! Raman active IREPs
logical :: raman_active(48)
! number of times an IREP occurs in the phonon symmetries
integer :: vib_ireps(48)
! effect of SOPs on atoms
integer, allocatable :: atom_sop(:, :)
! name of IREPs
character(4) :: irep_ch(48)
! is construction of symmetry vectors meaningful?
logical :: sym_out
! number of group of atoms that mix by SOPs
integer :: gr
! array containing atom number of group (index1) and group member (index2)
integer, allocatable :: gr_atoms(:, :)
! number of atoms of the groups
integer, allocatable :: gr_atoms_no(:)
! is a phonon mode belonging to this IREP?
logical :: vib_mode(48)
! max number of non-orthogonal vectors in symvec
integer :: no_vec_nonorth
end module raman_symmetry
!
