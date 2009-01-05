
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writegamma(gq)
use modmain
implicit none
! arguments
real(8), intent(in) :: gq(3*natmtot,nqpt)
! local variables
integer iq,i
open(50,file='GAMMAQ.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'(I4," : total number of atoms")') natmtot
write(50,'(I6," : number of q-points")') nqpt
write(50,*)
do iq=1,nqpt
  write(50,'(I6," : q-point")') iq
  write(50,'(3G18.10," : q-vector (lattice coordinates)")') vql(:,iq)
  write(50,'(3G18.10," : q-vector (Cartesian coordinates)")') vqc(:,iq)
  do i=1,3*natmtot
    write(50,'(I4,G18.10)') i,gq(i,iq)
  end do
  write(50,*)
end do
close(50)
write(*,*)
write(*,'("Info(writegamma):")')
write(*,'(" wrote phonon linewidths for all q-points to GAMMAQ.OUT")')
write(*,*)
return
end subroutine

