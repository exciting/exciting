
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potnucl(ptnucl,nr,r,zn,vn)
implicit none
! arguments
logical, intent(in) :: ptnucl
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: zn
real(8), intent(out) :: vn(nr)
! local variables
integer ir
! nuclear radius constant in Bohr
real(8), parameter :: r0=1.25d-15/0.52917720859d-10
real(8) rn,t1,t2
if (zn.eq.0.d0) then
  vn(:)=0.d0
  return
end if
if (ptnucl) then
! nucleus is taken to be a point particle
  do ir=1,nr
    vn(ir)=zn/r(ir)
  end do
else
! nucleus has a finite radius approximated by r0*A^(1/3)
  rn=r0*abs(zn)**(1.d0/3.d0)
  t1=zn/(2.d0*rn**3)
  t2=3.d0*rn**2
  do ir=1,nr
    if (r(ir).lt.rn) then
      vn(ir)=t1*(t2-r(ir)**2)
    else
      vn(ir)=zn/r(ir)
    end if
  end do
end if
return
end subroutine

