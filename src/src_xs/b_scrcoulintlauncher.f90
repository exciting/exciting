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
  character(256) :: casestring
  character(*), parameter :: thisname = "b_scrcoulintlauncher"

  write(*,*) "b_scrcoulintlauncher here at rank", rank

  ! Calculate RR, RA or RR and RA blocks
  casestring = input%xs%bse%blocks

  ! Also calculate coupling blocks
  if(input%xs%bse%coupling .eqv. .true.) then 
    fcoup = .true.
  else
    fcoup = .false.
    if(trim(casestring) /= "rr") then 
      ! silently set default from default of "both" to "rr"
      casestring="rr"
    end if
  end if

  ! Use time inverted anti-resonant basis
  if(input%xs%bse%ti .eqv. .true.) then 
    fti = .true.
  else
    fti = .false.
  end if

  ! Which Q points to consider 
  iqmti = 1
  iqmtf = size(input%xs%qpointset%qpoint, 2)
  if(input%xs%bse%iqmtrange(1) /= -1) then 
    iqmti=input%xs%bse%iqmtrange(1)
    iqmtf=input%xs%bse%iqmtrange(2)
  end if

  call printline(unitout, "+")
  write(unitout, '("Info(",a,"):", a)') trim(thisname),&
    & " Setting up screened interaction matrix."
  write(unitout, '("Info(",a,"):", a, i3, a, i3)') trim(thisname),&
    & " Using momentum transfer vectors from list : ", iqmti, " to", iqmtf
  call printline(unitout, "+")

  do iqmt = iqmti, iqmtf

    vqmt(:) = input%xs%qpointset%qpoint(:, iqmt)

    if(mpiglobal%rank == 0) then
      write(unitout, *)
      call printline(unitout, "+")
      write(unitout, '("Info(",a,"):", a)') trim(thisname), &
        & " Calculating screened Coulomb interaction matrix W"
      write(unitout, '("Info(",a,"):", a, i3)') trim(thisname), &
        & " Momentum tranfer list index: iqmt=", iqmt
      write(unitout, '("Info(",a,"):", a, 3f8.3)') trim(thisname), &
        & " Momentum tranfer: vqmtl=", vqmt(1:3)
      call printline(unitout, "+")
      write(unitout,*)
    end if

    select case(trim(casestring))

      case("RR","rr")

        ! RR block
        if(mpiglobal%rank == 0) then
          write(unitout, '("Info(",a,"):&
            & Calculating RR block of W")') trim(thisname) 
          write(unitout,*)
        end if
        ! b_scrcoulint(iqmt, fra=.false., fti=.false.)
        call b_scrcoulint(iqmt, .false., fti)
        call barrier(mpiglobal, callername=trim(thisname))

      case("RA","ra")

        ! RA block
        if(fcoup) then 
          if(mpiglobal%rank == 0) then
            call printline(unitout, "-")
            write(unitout, '("Info(",a,"):&
              & Calculating RA block of W")') trim(thisname) 
            if(fti) then
              write(unitout, '("Info(",a,"):&
                & Using time inverted anti-resonant basis")') trim(thisname) 
            end if
            write(unitout,*)
          end if
          call b_scrcoulint(iqmt, .true., fti)
          call barrier(mpiglobal, callername=trim(thisname))
        end if

      case("both","BOTH")

        ! RR block
        if(mpiglobal%rank == 0) then
          write(unitout, '("Info(",a,"):&
            & Calculating RR block of W")') trim(thisname) 
          write(unitout,*)
        end if
        ! b_scrcoulint(iqmt, fra=.false., fti=.false.)
        call b_scrcoulint(iqmt, .false., fti)
        call barrier(mpiglobal, callername=trim(thisname))

        ! RA block
        if(fcoup) then 
          if(mpiglobal%rank == 0) then
            call printline(unitout, "-")
            write(unitout, '("Info(",a,"):&
              & Calculating RA block of W")') trim(thisname) 
            if(fti) then
              write(unitout, '("Info(",a,"):&
                & Using time inverted anti-resonant basis")') trim(thisname) 
            end if
            write(unitout,*)
          end if
          call b_scrcoulint(iqmt, .true., fti)
          call barrier(mpiglobal, callername=trim(thisname))
        end if

      case default

        write(*,'("Error(",a,"): Unrecongnized casesting:", a)')&
          & trim(thisname), trim(casestring)
        call terminate

    end select

    if(mpiglobal%rank == 0) then
      call printline(unitout, "+")
      write(unitout, '("Info(",a,"): Screened coulomb interaction&
        & finished for iqmt=",i4)') trim(thisname), iqmt
      call printline(unitout, "+")
    end if

  end do

end subroutine b_scrcoulintlauncher
!EOC

