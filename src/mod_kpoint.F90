

#include "maxdefinitions.inc"
module mod_kpoint
!-------------------------------!
!     k-point set variables     !
!-------------------------------!
! autokpt is .true. if the k-point set is determined automatically
!replaced by inputstructurelogical::autokpt
! radius of sphere used to determine k-point density when autokpt is .true.
!replaced by inputstructurereal(8)::radkpt
! k-point grid sizes
!replaced by inputstructureinteger::ngridk(3)
! total number of k-points
integer::nkpt
! k-point offset
!replaced by inputstructurereal(8)::vkloff(3)
! reducek is .true. if k-points are to be reduced (with crystal symmetries)
!replaced by inputstructurelogical::reducek
! locations of k-points on integer grid
integer, allocatable :: ivk(:, :)
! k-points in lattice coordinates
real(8), allocatable :: vkl(:, :)
! k-points in Cartesian coordinates
real(8), allocatable :: vkc(:, :)
! k-point weights
real(8), allocatable :: wkpt(:)
! map from non-reduced grid to reduced set
integer, allocatable :: ikmap(:, :, :)
! total number of non-reduced k-points
integer::nkptnr
! locations of non-reduced k-points on integer grid
integer, allocatable :: ivknr(:, :)
! non-reduced k-points in lattice coordinates
real(8), allocatable :: vklnr(:, :)
! non-reduced k-points in Cartesian coordinates
real(8), allocatable :: vkcnr(:, :)
! non-reduced k-point weights
real(8), allocatable :: wkptnr(:)
! map from non-reduced grid to non-reduced set
integer, allocatable :: ikmapnr(:, :, :)
! k-point at which to determine effective mass tensor
!replaced by inputstructurereal(8)::vklem(3)
! displacement size for computing the effective mass tensor
!replaced by inputstructurereal(8)::deltaem
! number of displacements in each direction
!replaced by inputstructureinteger::ndspem
end module

