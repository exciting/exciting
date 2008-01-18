
! Copyright (C) 2002-2005 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine linoptkpq(ik,i1,i2,qmat,e,f)
use modmain
use modxs
use m_genfilname
implicit none
! arguments
integer, intent(in) :: ik
integer, intent(in) :: i1
integer, intent(in) :: i2
complex(8), intent(in) :: qmat(nstsv,nstsv)
real(8), intent(inout) :: e(nstsv*nstsv)
real(8), intent(inout) :: f(nstsv*nstsv)
! local variables
integer, parameter :: iq=1
integer ist,jst,isym,i,m
real(8) sum
complex(8) zt1(3)

integer :: ikq

!k+q-point
ikq=ikmapikq(ik,iq)
call getevalsv(vkl(1,ikq),evalsv(1,ikq))
call getoccsv(vkl(1,ikq),occsv(1,ikq))

!k-point
call genfilname(iq=0,setfilext=.true.)
call getevalsv0(vkl0(1,ik),evalsv0(1,ik))
call getoccsv0(vkl0(1,ik),occsv0(1,ik))
call genfilname(iq=iq,setfilext=.true.)



m=0
do ist=1,nstsv
  do jst=1,nstsv
    m=m+1


!!$! symmetrise the matrix elements
!!$    sum=0.d0
!!$    do isym=1,nsymcrys
!!$      zt1(:)=0.d0
!!$      if (i1.eq.i2) then
!!$        do i=1,3
!!$          zt1(i1)=zt1(i1)+sc(i,i1,isym)*pmat(i,ist,jst)
!!$        end do
!!$      else
!!$        do i=1,3
!!$          zt1(i1)=zt1(i1)+sc(i,i1,isym)*pmat(i,ist,jst)
!!$          zt1(i2)=zt1(i2)+sc(i,i2,isym)*pmat(i,ist,jst)
!!$        end do
!!$      end if
!!$! calculate the desired dielectric tensor components
!!$      if (i1.eq.i2) then
!!$        sum=sum+dble(zt1(i1)*conjg(zt1(i2)))
!!$      else
!!$        if (d(isym).lt.0.d0) then
!!$          sum=sum+aimag(zt1(i2)*conjg(zt1(i1)))
!!$        else
!!$          sum=sum+aimag(zt1(i1)*conjg(zt1(i2)))
!!$        endif
!!$      end if
!!$! end of symmetrisation
!!$    end do
!!$    sum=sum/dble(nsymcrys)

    


    e(m)=evalsv0(ist,ik)-evalsv(jst,ikq)
! scissors correction
    if (evalsv0(ist,ik).gt.efermi) e(m)=e(m)+scissor
    if (evalsv(jst,ikq).gt.efermi) e(m)=e(m)-scissor
    f(m)=(occsv0(ist,ik)-occsv(jst,ikq))* &
         abs(qmat(ist,jst))**2/gqc(1,1)**2
  end do
end do
return
end subroutine

