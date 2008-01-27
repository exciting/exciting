
! Copyright (C) 2002-2005 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine linoptk(ik,i1,i2,sc,d,delta,pmat,e,f,pmatint)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
integer, intent(in) :: i1
integer, intent(in) :: i2
real(8), intent(in) :: sc(3,3,nsymcrys)
real(8), intent(in) :: d(nsymcrys)
real(8), intent(in) :: delta(nstsv,nstsv)
complex(8), intent(in) :: pmat(3,nstsv,nstsv)
real(8), intent(inout) :: e(nstsv*nstsv)
real(8), intent(inout) :: f(nstsv*nstsv)
real(8), intent(out) :: pmatint(nstsv)
! local variables
integer ist,jst,isym,i,m
real(8) sum
complex(8) zt1(3)
! get the eigenvalues and occupancies from file
call getevalsv(vkl(1,ik),evalsv(1,ik))
call getoccsv(vkl(1,ik),occsv(1,ik))
m=0
do ist=1,nstsv
  do jst=1,nstsv
    m=m+1
! symmetrise the matrix elements
    sum=0.d0
    do isym=1,nsymcrys
      zt1(:)=0.d0
      if (i1.eq.i2) then
        do i=1,3
          zt1(i1)=zt1(i1)+sc(i,i1,isym)*pmat(i,ist,jst)
        end do
      else
        do i=1,3
          zt1(i1)=zt1(i1)+sc(i,i1,isym)*pmat(i,ist,jst)
          zt1(i2)=zt1(i2)+sc(i,i2,isym)*pmat(i,ist,jst)
        end do
      end if
! calculate the desired dielectric tensor components
      if (i1.eq.i2) then
        sum=sum+dble(zt1(i1)*conjg(zt1(i2)))
      else
        if (d(isym).lt.0.d0) then
          sum=sum+aimag(zt1(i2)*conjg(zt1(i1)))
        else
          sum=sum+aimag(zt1(i1)*conjg(zt1(i2)))
        endif
      end if
! end of symmetrisation
    end do
    sum=sum/dble(nsymcrys)
    e(m)=evalsv(ist,ik)-evalsv(jst,ik)
! store the matrix elements for intraband transitions
    if (ist.eq.jst) pmatint(ist)=sum
! scissors correction
    if (evalsv(ist,ik).gt.efermi) e(m)=e(m)+scissor
    if (evalsv(jst,ik).gt.efermi) e(m)=e(m)-scissor
! generalised DFT correction
    if (usegdft) e(m)=e(m)+delta(ist,jst)
    f(m)=(occsv(ist,ik)-occsv(jst,ik))*sum
  end do
end do
return
end subroutine

