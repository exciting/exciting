! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: exccoulintlauncher
! !INTERFACE:
subroutine exccoulintlauncher
! !USES:
  use modmpi
  use modxs, only: unitout
  use modinput, only: input
  use modbse
! !DESCRIPTION:
!   Launches the calculation of the exchange term of the Bethe-Salpeter Hamiltonian
!   for the specified momentum transfer vectors $\vec{Q}_\text{mt}$.
!
! !REVISION HISTORY:
!   Created. 2016 (Aurich)
!EOP
!BOC      

  implicit none

  logical :: fcoup, fchibarq
  integer(4) :: iqmt, iqmti, iqmtf, nqmt
  real(8) :: vqmt(3)
  character(256) :: casestring
  character(*), parameter :: thisname = "exccoulintlauncher"

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

  ! For BSE with coupling at finite Q, always use the full Chi
  ! and not \bar{Chi}. For TDA warn if \bar{Chi} is not used.
  fchibarq = input%xs%bse%chibarq
  ! If it is not the default of true
  if(.not. fchibarq) then 
    if(.not. fcoup) then
      write(unitout, '("Waring(",a,"):", a)') trim(thisname),&
        & " TDA and full Chi produce bad results for finite Q, use \bar{Chi}!&
        &   set input%xs%bse%chibarq = .true."
    end if
  ! Otherwise determin value with fcoup
  else
    ! Use full Chi if ra coupling is used
    if(fcoup) fchibarq = .false.
    ! Use truncated Coulomb potential (\bar{Chi}) for TDA
    if(.not. fcoup) fchibarq = .true.
    ! Overwrite input 
    input%xs%bse%chibarq = fchibarq
  end if

  ! Which momentum transfer Q points to consider 
  nqmt = size(input%xs%qpointset%qpoint, 2)
  iqmti = 1
  iqmtf = nqmt
  !   or use selected range
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
    & " Setting up exchange interaction matrix V."
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
        & "Calculating exchange interaction matrix V"
      write(unitout, '("Info(",a,"):", a, i3)') trim(thisname), &
        & " Momentum transfer list index: iqmt=", iqmt
      write(unitout, '("Info(",a,"):", a, 3f8.3)') trim(thisname), &
        & " Momentum transfer: vqmtl=", vqmt(1:3)
      call printline(unitout, "+")
      write(unitout,*)
    end if

    select case(trim(casestring))

      case("RR","rr")

        ! RR block
        if(mpiglobal%rank == 0) then
          write(unitout, '("Info(",a,"):&
            & Calculating RR block of V")') trim(thisname)
          write(unitout,*)
        end if
        call exccoulint(iqmt)
        call barrier(mpiglobal, callername=trim(thisname))

      case("RA","ra")

        ! RA block
        if(fcoup) then 
          call printline(unitout, "-")
          write(unitout, '("Info(",a,"):&
            & RR = RA^{tr} no further calculation needed.")') trim(thisname)
        end if
        
      case("both","BOTH")

        ! RR block
        if(mpiglobal%rank == 0) then
          write(unitout, '("Info(",a,"):&
            & Calculating RR block of V")') trim(thisname)
          write(unitout,*)
        end if
        call exccoulint(iqmt)
        call barrier(mpiglobal, callername=trim(thisname))

        ! RA block
        if(fcoup) then 
          call printline(unitout, "-")
          write(unitout, '("Info(",a,"):&
            & RR = RA^{tr} no further calculation needed.")') trim(thisname)
        end if

      case default

        write(*,'("Error(",a,"): Unrecongnized casesting:", a)') trim(thisname),&
          & trim(casestring)
        call terminate

    end select

    ! Info out
    call printline(unitout, "+")
    write(unitout, '("Info(",a,"): Exchange interaction&
      & finished for iqmt=",i4)') trim(thisname),  iqmt
    call printline(unitout, "+")

  end do

end subroutine exccoulintlauncher
!EOC

