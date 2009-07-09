

#include "maxdefinitions.inc"
module mod_convergence
!-------------------------------!
!     convergence variables     !
!-------------------------------!
! maximum number of self-consistent loops
!replaced by inputstructureinteger::maxscl
! current self-consistent loop number
integer::iscl
! effective potential convergence tolerance
!replaced by inputstructurereal(8)::epspot
! energy convergence tolerance
!replaced by inputstructurereal(8)::epsengy
! force convergence tolerance
!replaced by inputstructurereal(8)::epsforce
!curent convergence
real(8)::currentconvergence
end module

