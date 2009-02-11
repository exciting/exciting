
! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetime(t,fname)
  use modmpi
  implicit none
  ! arguments
  real(8), intent(in) :: t
  character(*), intent(in) :: fname
  open(80,file=trim(fname),form='formatted',action='write',status='replace')
  write(80,*)
  write(80,'("Total time elapsed")')
  write(80,'(" Wall clock seconds : ",es18.10)') t
  write(80,'(" Wall clock hours   : ",f14.2)') t/3600.d0
  write(80,'(" Wall clock days    : ",f14.2)') t/3600.d0/24.d0
  write(80,*)
#ifdef MPI
  write(80,'("Total time elapsed on all processors")')
  write(80,'(" Wall clock seconds : ",es18.10)') procs*t    *10000
  write(80,'(" Wall clock hours   : ",f14.2)') procs*t/3600.d0    *10000
  write(80,'(" Wall clock days    : ",f14.2)') procs*t/3600.d0/24.d0    *10000
  write(80,*)
#endif
  close(80)
end subroutine writetime
