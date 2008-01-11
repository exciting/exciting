
! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genparidxran(typ)
  use modmain
  use modmpi
  use modxs
  implicit none
  ! arguments
  character(1), intent(in) :: typ
  select case (typ)
     case('w')
        wpari=firstofset(rank,nwdf)
        wparf=lastofset(rank,nwdf)
        qpari=1
        qparf=nqpt
        kpari=1
        kparf=nkpt
     case('q')
        wpari=1
        wparf=nwdf
        qpari=firstofset(rank,nqpt)
        qparf=lastofset(rank,nqpt)        
        kpari=1
        kparf=nkpt
     case('k')
        wpari=1
        wparf=nwdf
        qpari=1
        qparf=nqpt
        kpari=firstofset(rank,nkpt)
        kparf=lastofset(rank,nkpt)
     case default
        write(*,*)
        write(*,'("Error(genparidxran): unknown parallelization type: ",a)') typ
        write(*,*)
        call terminate
  end select
end subroutine genparidxran
