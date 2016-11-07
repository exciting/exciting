! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!

!BOP
! !ROUTINE: genparidxran
! !INTERFACE:
subroutine genparidxran(typ, n)
! !USES:
  use modmpi
  use mod_qpoint, only: nqpt
  use mod_kpoint, only: nkpt
  use modxs, only: wpari, wparf, qpari, qparf,&
                 & kpari, kparf, nwdf, ppari,&
                 & pparf, partype, unitout
! !INPUT/OUTPUT PARAMETER:
!   IN:
!   character(1) :: typ  ! Physical meaning of index
!   integer(4) :: n      ! Number of elements to distribute
!   Module IN:
!   integer(4) :: nwdf   ! Number of frequencies for the construction of 
!                        ! the dielectric matrix in RPA
!   integer(4) :: nkpt   ! Number of k-points
!   integer(4) :: nqpt   ! Number of q-points 
!
! !DESCRIPTION:
!   The routine sets loop indices for each calling {\tt MPI} rank
!   for k-points, q-points, p-point and $\omega$-points in modxs.
!
! !REVISION HISTORY:
!   Added to documentation scheme. 2016 (Aurich)
!   Changed to also work for the case that there
!   are more processes than elements. (Aurich)
!EOP
!BOC

  implicit none

  ! Arguments
  character(1), intent(in) :: typ
  integer, intent(in) :: n

  ! Local variables
  integer :: np

  ! Warn if number of processors is greater than set
  if(procs .gt. n) then
    if(rank == 0) then 
      write(unitout, '("Warning(genparidxran):&
        & number of processors exceeds size of set")')
      write(unitout, '("  parallelization type : ", a)') typ
      write(unitout, '("  size of set	       : ", i6)') n
      write(unitout, '("  number of processors : ", i6)') procs
    end if
  end if

  ! Default values
  wpari = 1
  wparf = nwdf
  qpari = 1
  qparf = nqpt
  kpari = 1
  kparf = nkpt

  ! Number of (k,kp) pairs, where ik'>=ik
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
!EOC
