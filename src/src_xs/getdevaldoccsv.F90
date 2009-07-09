


! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine getdevaldoccsv(iq, ik, ikq, l1, u1, l2, u2, devalsv, doccsv, scissv)
  ! xssave0 has to be called in advance.
use modinput
  use modmain
  use modxs
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(in) :: iq, ik, ikq, l1, u1, l2, u2
  real(8), intent(out) :: devalsv(u1-l1+1, u2-l2+1)
  real(8), intent(out) :: doccsv(u1-l1+1, u2-l2+1)
  real(8), intent(out) :: scissv(u1-l1+1, u2-l2+1)
  ! local variables
  integer :: ist, jst, iqt
  real(8), allocatable :: e0(:), e(:), o0(:), o(:)
  iqt=iq
  allocate(e0(nstsv), e(nstsv), o0(nstsv), o(nstsv))
  ! eigenvalues and occupancies for k+q-point
  call getevalsv(vkl(1, ikq), e)
  call getoccsv(vkl(1, ikq), o)
  ! eigenvalues and occupancies for k-point
  call getevalsv0(vkl0(1, ik), e0)
  call getoccsv0(vkl0(1, ik), o0)
  ! scissors correction
  scissv(:, :)=0.d0
  do ist=l1, u1
     do jst=l2, u2
        ! neglect contributions above cutoff
	if ((e0(ist).lt.input%xs%emaxdf).and.(e(jst).lt.input%xs%emaxdf)) then
	   devalsv(ist-l1+1, jst-l2+1)=e0(ist)-e(jst)
	   doccsv(ist-l1+1, jst-l2+1)=o0(ist)-o(jst)
	   if ((e0(ist).le.efermi).and.(e(jst).gt.efermi)) &
		scissv(ist - l1 + 1, jst - l2 + 1) =- input%properties%bandstructure%scissor
	   if ((e0(ist).gt.efermi).and.(e(jst).le.efermi)) &
		scissv(ist - l1 + 1, jst - l2 + 1) = input%properties%bandstructure%scissor
	else
	  if (input%xs%dbglev.gt.1) then
	   write(*, '("Info(getdevaldoccsv): cutoff applied: iq, ik, ist, jst, &
		&energies:", 4i6, 2g18.10)') iq, ik, ist, jst, e0(ist), e(jst)
	  end if
           ! set energy difference to arbitrary number
	   devalsv(ist-l1+1, jst-l2+1)=1.d0
           ! set occupation number difference to zero as a factor in the
           ! terms of the dielectric function
	   doccsv(ist-l1+1, jst-l2+1)=0.d0
	end if
     end do
  end do
  deallocate(e0, e, o0, o)
end subroutine getdevaldoccsv
