

#include "maxdefinitions.inc"
module mod_LDA_LU
!-------------------------!
!     LDA+U variables     !
!-------------------------!
! type of LDA+U to use (0: none)
integer::ldapu
! maximum angular momentum
integer, parameter :: lmaxlu=3
integer, parameter :: lmmaxlu=(lmaxlu+1)**2
! angular momentum for each species
integer::llu(_MAXSPECIES_)
! U and J values for each species
real(8)::ujlu(2, _MAXSPECIES_)
! LDA+U density matrix
complex(8), allocatable :: dmatlu(:, :, :, :, :)
! LDA+U potential matrix in (l,m) basis
complex(8), allocatable :: vmatlu(:, :, :, :, :)
! LDA+U energy for each atom
real(8), allocatable :: engyalu(:)
! interpolation constant alpha for each atom (PRB 67, 153106 (2003))
real(8), allocatable :: alphalu(:)
! energy from the LDA+U correction
real(8)::engylu
end module

