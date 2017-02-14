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

  logical :: fcoup, fti
  integer(4) :: iqmt, iqmti, iqmtf
  real(8) :: vqmt(3)

  ! Also calculate coupling blocks
  if(input%xs%bse%coupling == .true.) then 
    fcoup = .true.
  else
    fcoup = .false.
  end if

  ! Use time inverted anti-resonant basis
  if(input%xs%bse%ti == .true.) then 
    fti = .true.
  else
    fti = .false.
  end if

  ! Testing 
  iqmti = 1
  iqmtf = size(input%xs%qpointset%qpoint, 2)

  do iqmt = iqmti, iqmtf

    vqmt(:) = input%xs%qpointset%qpoint(:, iqmt)

    if(mpiglobal%rank == 0) then
      write(unitout, *)
      call printline(unitout, "+")
      write(unitout, '("Info(b_scrcoulintlauncher):", a)')&
        & " Calculating screened Coulomb interaction matrix W"
      write(unitout, '("Info(b_scrcoulintlauncher):", a, i3)')&
        & " Momentum tranfer list index: iqmt=", iqmt
      write(unitout, '("Info(b_scrcoulintlauncher):", a, 3f8.3)')&
        & " Momentum tranfer: vqmtl=", vqmt(1:3)
      call printline(unitout, "+")
      write(unitout,*)
    end if

    ! RR block
    if(mpiglobal%rank == 0) then
      write(unitout, '("Info(b_scrcoulintlauncher):&
        & Calculating RR block of W")')
      write(unitout,*)
    end if
    ! b_scrcoulint(iqmt, fra=.false., fti=.false.)
    call b_scrcoulint(iqmt, .false., fti)
    call barrier(mpiglobal)

    ! RA block
    if(fcoup) then 
      if(mpiglobal%rank == 0) then
        call printline(unitout, "-")
        write(unitout, '("Info(b_scrcoulintlauncher):&
          & Calculating RA block of W")')
        if(fti) then
          write(unitout, '("Info(b_scrcoulintlauncher):&
            & Using time inverted anti-resonant basis")')
        end if
        write(unitout,*)
      end if
      call b_scrcoulint(iqmt, .true., fti)
      call barrier(mpiglobal)
    end if

    if(mpiglobal%rank == 0) then
      call printline(unitout, "+")
      write(unitout, '("Info(b_scrcoulintlauncher): Screened coulomb interaction&
        & finished for iqmt=",i4)') iqmt
      call printline(unitout, "+")
    end if

  end do

end subroutine b_scrcoulintlauncher
!EOC

