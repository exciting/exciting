! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_scrcoulintlauncher
! !INTERFACE:
subroutine b_scrcoulintlauncher
! !USES:
  use modmpi
  use modxs, only: unitout
  use modinput, only: input
  use modbse
! !DESCRIPTION:
!   Launches the calculation of the direct term of the Bethe-Salpeter Hamiltonian.
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

  ! Note: Only iqmt=0 is supported currently
  iqmt = 0

  ! RR block
  if(mpiglobal%rank == 0) then
    write(unitout, '("Info(b_scrcoulintlauncer):&
      & Calculating RR block of W for qmt=0")')
  end if
  ! b_scrcoulint(iqmt, fra=.false., fti=.false.)
  call b_scrcoulint(iqmt, .false., .false.)
  call barrier(mpiglobal)

  ! RA block
  if(fra) then 
    if(mpiglobal%rank == 0) then
      write(unitout, '("Info(b_scrcoulintlauncer):&
        & Calculating RA block of W for qmt=0")')
      if(fti) then
        write(unitout, '("Info(b_scrcoulintlauncer):&
          & Using time inverted anti-resonant basis")')
      end if
    end if
    call b_scrcoulint(iqmt, .true., fti)
    call barrier(mpiglobal)
  end if

end subroutine b_scrcoulintlauncher
!EOC

