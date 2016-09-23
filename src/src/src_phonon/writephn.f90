
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writephn
  use modmain
  use modmpi
  implicit none
! writeout only in master process
  if (rank .ne. 0) goto 10
! initialise universal variables
  Call init0
  Call init2
  call writephnlist(nphwrt,vqlwrt,.true.,"PHONON.OUT")
10 continue
#ifdef MPI
  call MPI_Barrier(MPI_Comm_World, ierr)
#endif
  return
end subroutine
