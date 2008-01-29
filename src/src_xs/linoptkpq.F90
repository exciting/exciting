
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine linoptkpq(ik,qmat,e,f)
  use modmain
  use modxs
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: ik
  complex(8), intent(in) :: qmat(nstsv,nstsv,ngq(iq))
  real(8), intent(inout) :: e(nstsv*nstsv)
  real(8), intent(inout) :: f(nstsv*nstsv)
  ! local variables
  integer, parameter :: iq=1
  integer ist,jst,m
  integer :: ikq
  ! eigenvalues and occupancies for k+q-point
  ikq=ikmapikq(ik,iq)
  call getevalsv(vkl(1,ikq),evalsv(1,ikq))
  call getoccsv(vkl(1,ikq),occsv(1,ikq))
  ! eigenvalues and occupancies for k-point
  call genfilname(iq=0,setfilext=.true.)
  call getevalsv0(vkl0(1,ik),evalsv0(1,ik))
  call getoccsv0(vkl0(1,ik),occsv0(1,ik))
  call genfilname(iq=iq,setfilext=.true.)
  m=0
  do ist=1,nstsv
     do jst=1,nstsv
        m=m+1
        e(m)=evalsv0(ist,ik)-evalsv(jst,ikq)
        ! scissors correction
        if (evalsv0(ist,ik).gt.efermi) e(m)=e(m)+scissor
        if (evalsv(jst,ikq).gt.efermi) e(m)=e(m)-scissor
        f(m)=(occsv0(ist,ik)-occsv(jst,ikq))* &
             abs(qmat(ist,jst,1))**2/gqc(1,1)**2
        if ((.not.intraband).and.(ist.eq.jst)) f(m)=0.d0
     end do
  end do
end subroutine linoptkpq
