! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_exccoulintlauncher
! !INTERFACE:
subroutine b_exccoulintlauncher
! !USES:
  use modmpi
  use modxs, only: unitout
  use modinput, only: input
  use modbse
! !DESCRIPTION:
!   Launches the calculation of the exchange term of the Bethe-Salpeter Hamiltonian.
!
! !REVISION HISTORY:
!   Created. 2016 (Aurich)
!EOP
!BOC      

  implicit none

  logical :: fra, fti
  integer(4) :: iqmt

  ! Also calculate coupling blocks
  if(input%xs%bse%coupling == .true.) then 
    fra = .true.
  else
    fra = .false.
  end if

  ! Use time inverted anti-resonant basis
  if(input%xs%bse%ti == .true.) then 
    fti = .true.
  else
    fti = .false.
  end if

  ! Note: Only iqmt=0 is supported 
  iqmt = 0

  ! RR block
  if(mpiglobal%rank == 0) then
    write(unitout, '("Info(b_exccoulintlauncher):&
      & Calculating RR block of V for qmt=0")')
  end if
  call b_exccoulint(iqmt, .false., .false.)
  call barrier(mpiglobal)

  ! RA block
  if(fra) then 
    if(mpiglobal%rank == 0) then
      write(unitout, '("Info(b_exccoulintlauncher):&
        & Calculating RA block of V for qmt=0")')
    end if
    if(fti) then
      write(unitout, '("Info(b_exccoulintlauncher):&
        & RR = RA^{ti} no further calculation needed.")')
    else
      call b_exccoulint(iqmt, .true., fti)
    end if
    call barrier(mpiglobal)
  end if

end subroutine b_exccoulintlauncher
!EOC

