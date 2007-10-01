real(8) function fopt(param)
use modeos
implicit none
! arguments
real(8), intent(in) :: param(*)
! local variables
integer ipt
real(8) rsum
! external functions
real(8) eveos
external eveos
rsum=0.d0
do ipt=1,nevpt
  rsum=rsum+(eveos(etype,param,vpt(ipt))-ept(ipt))**2
end do
fopt=rsum
return
end function
