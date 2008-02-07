
subroutine hamiltonandoverlapsetupnotpacked(n,ngp,apwalm,igpig,vgpc)
use modmain, only:ngkmax,apwordmax,lmmaxapw,natmtot
use diisinterfaces
use modfvsystem
implicit none
integer, intent(in)::n,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)



integer ::ip, icolumn,irow
ohrank=n


end subroutine 