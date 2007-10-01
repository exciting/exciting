
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepvnl(vnlcv,vnlvv)
use modmain
implicit none
! arguments
complex(8), intent(out) :: vnlcv(ncrmax,natmtot,nstsv,nkpt)
complex(8), intent(out) :: vnlvv(nstsv,nstsv,nkpt)
! local variables
integer ik
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkpt
!$OMP CRITICAL
  write(*,'("Info(oepvnl): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
  call oepvnlk(ik,vnlcv(1,1,1,ik),vnlvv(1,1,ik))
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

