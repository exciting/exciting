real(8) function pveos(etype,param,v)
! pressure-volume equation of state function
implicit none
! arguments
integer, intent(in) :: etype
real(8), intent(in) :: param(*)
real(8), intent(in) :: v
! local variables
real(8) vm,vp,pm,pp,dv
! external functions
real(8) eveos
external eveos
! use central differences
dv=1.d-3
vm=v-dv
vp=v+dv
pm=eveos(etype,param,vm)
pp=eveos(etype,param,vp)
pveos=-(pp-pm)/(2.d0*dv)
return
end function
