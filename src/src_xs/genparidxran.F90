! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
subroutine genparidxran(typ, n)
  use modmpi, only: procs, rank, firstofset, lastofset
  use mod_qpoint, only: nqpt
  use mod_kpoint, only: nkpt
  use modxs, only: wpari, wparf, qpari, qparf,&
                 & kpari, kparf, nwdf, ppari,&
                 & pparf, partype

  implicit none

  ! Arguments
  character(1), intent(in) :: typ
  integer, intent(in) :: n

  ! Local variables
  integer :: np

  ! Check if number of processors is greater than set
  if(procs .gt. n) then
    write(*,*)
    write(*, '("Error(genparidxran): number of processors exceeds size of set")')
    write(*, '(" parallelization type : ", a)') typ
    write(*, '(" size of set	       : ", i6)') n
    write(*, '(" number of processors : ", i6)') procs
    write(*,*)
    call terminate
  end if

  ! Default values
  wpari = 1
  wparf = nwdf
  qpari = 1
  qparf = nqpt
  kpari = 1
  kparf = nkpt

  ! Number of(k,kp) pairs
  np = nkpt * (nkpt+1) / 2
  ppari = 1
  pparf = np

  select case(typ)
    case('w')
      wpari = firstofset(rank, n)
      wparf = lastofset(rank, n)
    case('q')
      qpari = firstofset(rank, n)
      qparf = lastofset(rank, n)
    case('k')
      kpari = firstofset(rank, n)
      kparf = lastofset(rank, n)
    case('p')
      ppari = firstofset(rank, n)
      pparf = lastofset(rank, n)
    case default
      write(*,*)
      write(*, '("Error(genparidxran):&
        & unknown parallelization type: ", a)') typ
      write(*,*)
      call terminate
  end select

  partype = typ

end subroutine genparidxran
