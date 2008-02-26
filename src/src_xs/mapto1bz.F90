
! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mapto1bz(vl,vl1bz,iv)
  use modmain
  real(8), intent(in) :: vl(3)
  real(8), intent(out) :: vl1bz(3)
  integer, intent(out) :: iv(3)
  ! local variables
  real(8) :: v0(3),v1(3)
  ! map the q-vector into the first Brillouin zone
  t1=1.d8
  do i1=-1,1
     do i2=-1,1
        do i3=-1,1
           v0=vl+dble((/i1,i2,i3/))
           v1=matmul(bvec,v0)
           t2=v1(1)**2+v1(2)**2+v1(3)**2
           ! favour positive coordinates
           if (t2.lt.(t1+1.d-8)) then
              t1=t2
              vl1bz(:)=v0(:)
              iv=-(/i1,i2,i3/)
           end if
        end do
     end do
  end do
end subroutine mapto1bz
