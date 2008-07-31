
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine terminate
  use modmpi
  implicit none
  ! abort MPI if necessary
#ifdef MPI
  call mpi_abort(mpi_comm_world, ierr)
  if (ierr.eq.0) then
     write(*,'(a)') 'MPI abort clean. STOP'
  else
     write(*,'(a)') 'MPI abort unclean - zombie processes might remain! STOP'
  end if
#endif
#ifndef MPI
  write(*,'(a)') 'Abort clean. STOP'
#endif
  ! stop
  stop
end subroutine terminate

subroutine terminate_inqr(str)
  use m_getunit
  implicit none
  ! arguments
  character(*) :: str
  ! local variables
  logical :: exis
  integer :: un
  ! check for terminating
  inquire(file='TERMINATE',exist=exis)
  if (exis) then
     write(*,'(a)') 'Error: user termination of program in routine: '// &
          trim(str)
     call getunit(un)
     open(un,file='TERMINATE',action='write')
     close(un,status='delete')
     call terminate
  end if
end subroutine terminate_inqr
