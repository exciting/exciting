

#include "maxdefinitions.inc"
module mod_RDMFT
!-------------------------------------------------------------!
!     reduced density matrix functional (RDMFT) variables     !
!-------------------------------------------------------------!
! non-local matrix elements for varying occupation numbers
real(8), allocatable :: vnlrdm(:, :, :, :)
! Coulomb potential matrix elements
complex(8), allocatable :: vclmat(:, :, :)
! derivative of kinetic energy w.r.t. natural orbital coefficients
complex(8), allocatable :: dkdc(:, :, :)
! step size for occupation numbers
real(8)::taurdmn
! step size for natural orbital coefficients
real(8)::taurdmc
! xc functional
integer::rdmxctype
! maximum number of self-consistent loops
integer::rdmmaxscl
! maximum number of iterations for occupation number optimisation
integer::maxitn
! maximum number of iteration for natural orbital optimisation
integer::maxitc
! exponent for the functional
real(8)::rdmalpha
! temperature
real(8)::rdmtemp
! entropy
real(8)::rdmentrpy
end module
