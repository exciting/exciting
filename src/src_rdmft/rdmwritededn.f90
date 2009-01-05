
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine rdmwritededn(dedn)
! writes derivative of total energy w.r.t. occupancies to file
use modmain
implicit none
! arguments
real(8), intent(in) :: dedn(nstsv,nkpt)
! local variables
integer ik,ist
open(50,file='RDM_DEDN.OUT',action='WRITE',form='FORMATTED')
write(50,'(I6," : nkpt")') nkpt
write(50,'(I6," : nstsv")') nstsv
do ik=1,nkpt
  write(50,*)
  write(50,'(I6,3G18.10," : k-point, vkl")') ik,vkl(:,ik)
  write(50,'("     (state, occupancy and derivative below)")')
  do ist=1,nstsv
    write(50,'(I6,4G18.10)') ist,occsv(ist,ik),-dedn(ist,ik)
  end do
end do
close(50)
return
end subroutine

