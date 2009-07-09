


! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine terminate
  use modmpi
  implicit none
  ! abort MPI if necessary
#ifdef MPI
  call mpi_abort(mpi_comm_world, 1, ierr)
  if (ierr.eq.0) then
     write(*, '(a)') 'MPI abort'
  else
     write(*, '(a)') 'MPI abort with errors - zombie processes might remain!'
  end if
#endif
#ifndef MPI
  write(*, '(a)') 'Abort'
#endif
  ! stop program
  stop
end subroutine terminate
