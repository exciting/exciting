

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
!replaced by inputstructurereal(8)::taurdmn
! step size for natural orbital coefficients
!replaced by inputstructurereal(8)::taurdmc
! xc functional
!replaced by inputstructureinteger::rdmxctype
! maximum number of self-consistent loops
!replaced by inputstructureinteger::rdmmaxscl
! maximum number of iterations for occupation number optimisation
!replaced by inputstructureinteger::maxitn
! maximum number of iteration for natural orbital optimisation
!replaced by inputstructureinteger::maxitc
! exponent for the functional
!replaced by inputstructurereal(8)::rdmalpha
! temperature
!replaced by inputstructurereal(8)::rdmtemp
! entropy
real(8)::rdmentrpy
end module

