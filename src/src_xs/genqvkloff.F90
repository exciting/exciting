


! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine genqvkloff(vq, voff)
  use modmain
use modinput
  use modxs
  implicit none
  ! arguments
  real(8), intent(in) :: vq(3)
  real(8), intent(out) :: voff(3)
  ! local variables
  real(8) :: v1(3)
  integer :: iv(3)
  if (any(input%groundstate%vkloff/dble(input%groundstate%ngridk)+vq.ge.1.d0)) then
     ! vector is outside Brillouine zone
     v1=input%groundstate%vkloff/dble(input%groundstate%ngridk)+vq
     call mapkto01(v1)
     voff=v1*dble(input%groundstate%ngridk)
     if (any(v1*dble(input%groundstate%ngridk).ge.1.d0)) then
	v1=v1*dble(input%groundstate%ngridk)
	call mapkto01(v1)
	voff=v1
     end if
  else if (any(input%groundstate%vkloff+vq*dble(input%groundstate%ngridk).ge.1.d0)) then
     ! vector is inside Brillouine zone but outside k-point spacing
     v1=input%groundstate%vkloff+vq*dble(input%groundstate%ngridk)
     call mapkto01(v1)
     voff=v1
  else
     ! vector is inside k-point spacing
     voff=input%groundstate%vkloff+vq*input%groundstate%ngridk
  end if
  ! treatment of values close to zero or one
  call r3frac(input%structure%epslat, voff, iv)
end subroutine genqvkloff


subroutine mapkto01(v)
  implicit none
  ! arguments
  real(8), intent(inout) :: v(3)
  ! local variables
  !integer :: id(3)
  integer(8) :: v2(3), v3(3)
  real(8), parameter :: fac=1.d15
  !call r3frac(epslat,v,id)
  v2=dint(v)
  v3=dint(fac*v)
  v3=v3-v2*dint(fac)
  v=dble(v3/dint(fac))
  where(v.lt.0.d0) v=v+1.d0
end subroutine mapkto01
