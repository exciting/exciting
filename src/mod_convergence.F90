

#include "maxdefinitions.inc"
module mod_convergence
!-------------------------------!
!     convergence variables     !
!-------------------------------!
! maximum number of self-consistent loops
integer::maxscl
! current self-consistent loop number
integer::iscl
! effective potential convergence tolerance
real(8)::epspot
! energy convergence tolerance
real(8)::epsengy
! force convergence tolerance
real(8)::epsforce
!curent convergence
real(8)::currentconvergence
end module
