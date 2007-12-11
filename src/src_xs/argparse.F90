
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine argparse()
  use modmain
  use modmpi
  implicit none
  ! local variables
  character(256) :: str
  
#ifndef MPI
  ! read rank
  call getarg(1,str)
  if (trim(str).eq.'') then
     rank=0
  else
     str=adjustl(str)
     read(str,*) rank
  end if

  ! read procs
  call getarg(2,str)
  if (trim(str).eq.'') then
     procs=1
  else
     str=adjustl(str)
     read(str,*) procs
  end if
#endif

end subroutine argparse
