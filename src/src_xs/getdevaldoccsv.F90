
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getdevaldoccsv(iq,ik,ikq,l1,u1,l2,u2,devalsv,doccsv,scissv)
  ! tdsave0 has to be called in advance.
  use modmain
  use modxs
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: iq,ik,ikq,l1,u1,l2,u2
  real(8), intent(out) :: devalsv(u1-l1+1,u2-l2+1)
  real(8), intent(out) :: doccsv(u1-l1+1,u2-l2+1)
  real(8), intent(out) :: scissv(u1-l1+1,u2-l2+1)
  ! local variables
  integer :: ist,jst
  real(8), allocatable :: e0(:),e(:),o0(:),o(:)
  ! check output array sizes
  allocate(e0(nstsv),e(nstsv),o0(nstsv),o(nstsv))
  ! eigenvalues and occupancies for k+q-point
!!!  call genfilname(iqmt=iq,setfilext=.true.)
  call getevalsv(vkl(1,ikq),e)
  call getoccsv(vkl(1,ikq),o)
  ! eigenvalues and occupancies for k-point
!!!  call genfilname(iqmt=0,setfilext=.true.)
  call getevalsv0(vkl0(1,ik),e0)
  call getoccsv0(vkl0(1,ik),o0)
!!!  call genfilname(iqmt=iq,setfilext=.true.)
  ! scissors correction
  scissv(:,:)=0.d0
  do ist=l1,u1
     do jst=l2,u2
        devalsv(ist-l1+1,jst-l2+1)=e0(ist)-e(jst)
        doccsv(ist-l1+1,jst-l2+1)=o0(ist)-o(jst)
        if ((e0(ist).le.efermi).and.(e(jst).gt.efermi)) &
             scissv(ist-l1+1,jst-l2+1)=-scissor
        if ((e0(ist).gt.efermi).and.(e(jst).le.efermi)) &
             scissv(ist-l1+1,jst-l2+1)=scissor
     end do
  end do
  deallocate(e0,e,o0,o)
end subroutine getdevaldoccsv
