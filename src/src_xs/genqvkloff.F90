
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genqvkloff(vq,voff)
  use modmain
  use modxs
  implicit none
  ! arguments
  real(8), intent(in) :: vq(3)
  real(8), intent(out) :: voff(3)
  ! local variables
  real(8) :: v2(3)
  if (any(vkloff/dble(ngridk)+vq.ge.1.d0)) then
     ! vector is outside Brillouine zone
     v2=vkloff/dble(ngridk)+vq
     call mapkto01(v2)
     voff=v2*dble(ngridk)
     if (any(v2*dble(ngridk).ge.1.d0)) then
        v2=v2*dble(ngridk)
        call mapkto01(v2)
        voff=v2
     end if
  else if (any(vkloff+vq*dble(ngridk).ge.1.d0)) then
     ! vector is inside Brillouine zone but outside k-point spacing
     v2=vkloff+vq*dble(ngridk)
     call mapkto01(v2)
     voff=v2
  else
     ! vector is inside k-point spacing
     voff=vkloff+vq*ngridk
  end if
end subroutine genqvkloff

subroutine mapkto01(v)
  implicit none
  ! arguments
  real(8), intent(inout) :: v(3)
  ! local variables
  !integer :: id(3)
  integer(8) :: v2(3),v3(3)
  real(8),parameter :: fac=1.d15
  !call r3frac(epslat,v,id)
  v2=dint(v)
  v3=dint(fac*v)
  v3=v3-v2*dint(fac)
  v=dble(v3/dint(fac))
  where(v.lt.0.d0) v=v+1.d0

end subroutine mapkto01
