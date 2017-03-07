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

  logical :: fcoup, fti
  integer(4) :: iqmt, iqmti, iqmtf
  character(256) :: casestring
  real(8) :: vqmt(3)

  ! Calculate RR, RA or RR and RA blocks
  casestring = input%xs%bse%blocks

  ! Also calculate coupling blocks
  if(input%xs%bse%coupling == .true.) then 
    fcoup = .true.
  else
    fcoup = .false.
    if(trim(casestring) /= "rr") then 
      write(*,*) "Ignoring input%xs%bse%blocks, since no RA coupling enabled"
      casestring="rr"
    end if
  end if

  ! Use time inverted anti-resonant basis
  if(input%xs%bse%ti == .true.) then 
    fti = .true.
  else
    fti = .false.
  end if

  ! Which Q points to consider 
  iqmti = 1
  iqmtf = size(input%xs%qpointset%qpoint, 2)
  if(input%xs%bse%iqmt /= -1) then 
    iqmti=input%xs%bse%iqmt
    iqmtf=input%xs%bse%iqmt
  end if

  do iqmt = iqmti, iqmtf

    vqmt(:) = input%xs%qpointset%qpoint(:, iqmt)

    if(mpiglobal%rank == 0) then
      write(unitout, *)
      call printline(unitout, "+")
      write(unitout, '("Info(b_exccoulintlauncher):", a)')&
        & "Calculating exchange interaction matrix V"
      write(unitout, '("Info(b_exccoulintlauncher):", a, i3)')&
        & " Momentum tranfer list index: iqmt=", iqmt
      write(unitout, '("Info(b_exccoulintlauncher):", a, 3f8.3)')&
        & " Momentum tranfer: vqmtl=", vqmt(1:3)
      call printline(unitout, "+")
      write(unitout,*)
    end if

    select case(trim(casestring))

      case("RR","rr")

        ! RR block
        if(mpiglobal%rank == 0) then
          write(unitout, '("Info(b_exccoulintlauncher):&
            & Calculating RR block of V")')
          write(unitout,*)
        end if
        call b_exccoulint(iqmt, .false., fti)
        call barrier(mpiglobal)

      case("RA","ra")

        ! RA block
        if(fcoup) then 
          if(mpiglobal%rank == 0) then
            call printline(unitout, "-")
            write(unitout, '("Info(b_exccoulintlauncher):&
              & Calculating RA block of V")')
          end if
          if(fti) then
            call printline(unitout, "-")
            write(unitout, '("Info(b_exccoulintlauncher):&
              & RR = RA^{ti} no further calculation needed.")')
          else
            call b_exccoulint(iqmt, .true., fti)
          end if
          call barrier(mpiglobal)
        end if
        
      case("both","BOTH")

        ! RR block
        if(mpiglobal%rank == 0) then
          write(unitout, '("Info(b_exccoulintlauncher):&
            & Calculating RR block of V")')
          write(unitout,*)
        end if
        call b_exccoulint(iqmt, .false., fti)
        call barrier(mpiglobal)

        ! RA block
        if(fcoup) then 
          if(mpiglobal%rank == 0) then
            call printline(unitout, "-")
            write(unitout, '("Info(b_exccoulintlauncher):&
              & Calculating RA block of V")')
          end if
          if(fti) then
            call printline(unitout, "-")
            write(unitout, '("Info(b_exccoulintlauncher):&
              & RR = RA^{ti} no further calculation needed.")')
          else
            call b_exccoulint(iqmt, .true., fti)
          end if
          call barrier(mpiglobal)
        end if

      case default

        write(*,*) "Error(b_exccoulintlauncher): Unrecongnized casesting:", trim(casestring)
        call terminate

    end select

    if(mpiglobal%rank == 0) then
      call printline(unitout, "+")
      write(unitout, '("Info(b_exccoulintlauncher): Exchange interaction&
        & finished for iqmt=",i4)') iqmt
      call printline(unitout, "+")
    end if

  end do

end subroutine b_exccoulintlauncher
!EOC

