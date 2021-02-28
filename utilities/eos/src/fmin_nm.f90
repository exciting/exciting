
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function fmin_nm(x)
use modmain
implicit none
! arguments
real(8), intent(in) :: x
! local variables
integer i
real(8) sum
! external functions
real(8) eveos
external eveos
sum=0.d0
do i=1,nevpt
  sum=sum+(eveos(etype,x,vpt(i))-ept(i))**2
end do
fmin_nm=sum
return
end function

