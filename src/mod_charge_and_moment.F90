

#include "maxdefinitions.inc"
module mod_charge_and_moment
!-------------------------------------!
!     charge and moment variables     !
!-------------------------------------!
! tolerance for error in total charge
real(8)::epschg
! total nuclear charge
real(8)::chgzn
! total core charge
real(8)::chgcr
! core leakage charge
real(8)::chgcrlk
! total valence charge
real(8)::chgval
! excess charge
real(8)::chgexs
! total charge
real(8)::chgtot
! calculated total charge
real(8)::chgcalc
! interstitial region charge
real(8)::chgir
! muffin-tin charges
real(8), allocatable :: chgmt(:)
! total muffin-tin charge
real(8)::chgmttot
! effective Wigner radius
real(8)::rwigner
! total moment
real(8)::momtot(3)
! interstitial region moment
real(8)::momir(3)
! muffin-tin moments
real(8), allocatable :: mommt(:, :)
! total muffin-tin moment
real(8)::mommttot(3)
end module
