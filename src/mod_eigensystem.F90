
#include "maxdefinitions.inc"
module mod_eigensystem
!-------------------------------------------!
!     overlap and Hamiltonian variables     !
!-------------------------------------------!
! order of overlap and Hamiltonian matrices for each k-point
integer, allocatable :: nmat(:, :)
! maximum nmat over all k-points
integer::nmatmax
! size of packed matrices
integer, allocatable :: npmat(:, :)
! index to the position of the local-orbitals in the H and O matrices
integer, allocatable :: idxlo(:, :, :)
! APW-local-orbital overlap integrals
real(8), allocatable :: oalo(:, :, :)
! local-orbital-local-orbital overlap integrals
real(8), allocatable :: ololo(:, :, :)
! APW-APW Hamiltonian integrals
real(8), allocatable :: haa(:, :, :, :, :, :)
! local-orbital-APW Hamiltonian integrals
real(8), allocatable :: hloa(:, :, :, :, :)
! local-orbital-local-orbital Hamiltonian integrals
real(8), allocatable :: hlolo(:, :, :, :)
! complex Gaunt coefficient array
complex(8), allocatable :: gntyry(:, :, :)
! tseqit is .true. if the first-variational secular equation is to be solved
! iteratively
logical::tseqit
! number of secular equation iterations per self-consistent loop
integer::nseqit
! iterative solver step length
real(8)::tauseq
end module
