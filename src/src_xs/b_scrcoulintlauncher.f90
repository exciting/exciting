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
!   Launches the calculation of the direct term of the Bethe-Salpeter Hamiltonian
!   for each specified momentum transfer vector.
!
! !REVISION HISTORY:
!   Created. 2016 (Aurich)
!EOP
!BOC      

  implicit none

  logical :: fcoup
  integer(4) :: iqmt, iqmti, iqmtf, nqmt
  real(8) :: vqmt(3)
  character(256) :: casestring
  character(*), parameter :: thisname = "b_scrcoulintlauncher"

  ! Calculate RR or RA blocks
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

  ! Which momentum transfer Q points to consider 
  nqmt = size(input%xs%qpointset%qpoint, 2)
  iqmti = 1
  iqmtf = nqmt
  if(input%xs%bse%iqmtrange(1) /= -1) then 
    iqmti=input%xs%bse%iqmtrange(1)
    iqmtf=input%xs%bse%iqmtrange(2)
  end if
  if(iqmtf > nqmt .or. iqmti < -1 .or. iqmti > iqmtf) then 
    write(unitout, '("Error(",a,"):", a)') trim(thisname),&
      & " iqmtrange incompatible with qpointset list"
    call terminate
  end if

  ! Info output
  call printline(unitout, "+")
  write(unitout, '("Info(",a,"):", a)') trim(thisname),&
    & " Setting up screened interaction matrix."
  write(unitout, '("Info(",a,"):", a, i3, a, i3)') trim(thisname),&
    & " Using momentum transfer vectors from list : ", iqmti, " to", iqmtf
  call printline(unitout, "+")

  ! Loop over (subset of) Q points
  do iqmt = iqmti, iqmtf

    ! Get full Q vector for info out
    vqmt(:) = input%xs%qpointset%qpoint(:, iqmt)

    ! Info output
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
        call b_scrcoulint(iqmt, .false.)
        call barrier(mpiglobal, callername=trim(thisname))

      case("RA","ra")

        ! RA block
        if(fcoup) then 
          if(mpiglobal%rank == 0) then
            call printline(unitout, "-")
            write(unitout, '("Info(",a,"):&
              & Calculating RA block of W")') trim(thisname) 
            write(unitout,*)
          end if
          call b_scrcoulint(iqmt, .true.)
          call barrier(mpiglobal, callername=trim(thisname))
        end if

      case("both","BOTH")

        ! RR block
        if(mpiglobal%rank == 0) then
          write(unitout, '("Info(",a,"):&
            & Calculating RR block of W")') trim(thisname) 
          write(unitout,*)
        end if
        call b_scrcoulint(iqmt, .false.)
        call barrier(mpiglobal, callername=trim(thisname))

        ! RA block
        if(fcoup) then 
          if(mpiglobal%rank == 0) then
            call printline(unitout, "-")
            write(unitout, '("Info(",a,"):&
              & Calculating RA block of W")') trim(thisname) 
            write(unitout,*)
          end if
          call b_scrcoulint(iqmt, .true.)
          call barrier(mpiglobal, callername=trim(thisname))
        end if

      case default

        write(*,'("Error(",a,"): Unrecongnized casesting:", a)')&
          & trim(thisname), trim(casestring)
        call terminate

    end select

    ! Info out
    if(mpiglobal%rank == 0) then
      call printline(unitout, "+")
      write(unitout, '("Info(",a,"): Screened coulomb interaction&
        & finished for iqmt=",i4)') trim(thisname), iqmt
      call printline(unitout, "+")
    end if

  end do

end subroutine b_scrcoulintlauncher
!EOC

