

#include "maxdefinitions.inc"
module mod_eigenvalue_occupancy
use modinput
!--------------------------------------------!
!     eigenvalue and occupancy variables     !
!--------------------------------------------!
! number of empty states
!replaced by inputstructureinteger::nempty
! number of first-variational states
integer::nstfv
! number of second-variational states
integer::nstsv
! smearing type
!replaced by inputstructureinteger::stype
! smearing function description
character(256)::sdescr
! smearing width
!replaced by inputstructurereal(8)::swidth
! maximum allowed occupancy (1 or 2)
real(8)::occmax
! convergence tolerance for occupancies
!replaced by inputstructurereal(8)::epsocc
! second-variational occupation number array
real(8), allocatable :: occsv(:, :)
! Fermi energy for second-variational states
real(8)::efermi
! density of states at the Fermi energy
real(8)::fermidos
! error tolerance for the first-variational eigenvalues
!replaced by inputstructurereal(8)::evaltol
! minimum allowed eigenvalue
!replaced by inputstructurereal(8)::evalmin
! second-variational eigenvalues
real(8), allocatable :: evalsv(:, :)
! tevecsv is .true. if second-variational eigenvectors are calculated
!replaced by inputstructurelogical::tevecsv
! maximum number of k-point and states indices in user-defined list
integer, parameter :: maxkst=20
! number of k-point and states indices in user-defined list
integer::nkstlist
! user-defined list of k-point and state indices
integer::kstlist(3, maxkst)
end module

