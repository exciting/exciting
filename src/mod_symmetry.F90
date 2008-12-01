#include "maxdefinitions.inc"
module mod_symmetry
!----------------------------!
!     symmetry variables     !
!----------------------------!
! nosym is .true. if no symmetry information should be used
logical nosym
! number of Bravais lattice point group symmetries
integer nsymlat
! Bravais lattice point group symmetries
integer symlat(3,3,48)
! determinants of lattice symmetry matrices (1 or -1)
integer symlatd(48)
! index to inverses of the lattice symmetries
integer isymlat(48)
! lattice point group symmetries in Cartesian coordinates
real(8) symlatc(3,3,48)
! tshift is .true. if atomic basis is allowed to be shifted
logical tshift
! maximum of symmetries allowed
integer, parameter :: maxsymcrys=192
! number of crystal symmetries
integer nsymcrys
! crystal symmetry translation vector in lattice coordinates
real(8) vtlsymc(3,maxsymcrys)
! spatial rotation element in lattice point group for each crystal symmetry
integer lsplsymc(maxsymcrys)
! global spin rotation element in lattice point group for each crystal symmetry
integer lspnsymc(maxsymcrys)
! equivalent atom index for each crystal symmetry
integer, allocatable :: ieqatom(:,:,:)
! eqatoms(ia,ja,is) is .true. if atoms ia and ja are equivalent
logical, allocatable :: eqatoms(:,:,:)
! number of site symmetries
integer, allocatable :: nsymsite(:)
! site symmetry spatial rotation element in lattice point group
integer, allocatable :: lsplsyms(:,:)
! site symmetry global spin rotation element in lattice point group
integer, allocatable :: lspnsyms(:,:)
end module

