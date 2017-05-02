! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_bselauncher
! !INTERFACE:
subroutine b_bselauncher
! !USES:
  use modmpi
  use modxs, only: unitout
  use modinput, only: input
! !DESCRIPTION:
!   Launches the construction and solving of the Bethe-Salpeter Hamiltonian.
!
! !REVISION HISTORY:
!   Created. 2016 (Aurich)
!EOP
!BOC      

  implicit none

  integer(4) :: iqmt, iqmti, iqmtf, nqmt, iq1, iq2
  real(8) :: vqmt(3)
  character(*), parameter :: thisname = "b_bselauncher"

  ! Q-point list entries
  !   Use all
  iqmti = 1
  iqmtf = size(input%xs%qpointset%qpoint, 2)
  !   or use only one
  if(input%xs%bse%iqmtrange(1)/= -1) then 
    iqmti=input%xs%bse%iqmtrange(1)
    iqmtf=input%xs%bse%iqmtrange(2)
  end if
  nqmt = iqmtf-iqmti+1

  call printline(unitout, "+")
  write(unitout, '("Info(",a,"):", a)') trim(thisname),&
    & " Setting up and diagonalizing BSE Hamiltonian."
  write(unitout, '("Info(",a,"):", a, i3, a, i3)') trim(thisname),&
    & " Using momentum transfer vectors from list : ", iqmti, " to", iqmtf
  call printline(unitout, "+")

  if(.not. input%xs%bse%distribute) then 
    write(unitout, '("Info(",a,"):", a, i3, a)') trim(thisname),&
      & " Distributing qmt-points over ", mpiglobal%procs, " processes."
    call printline(unitout, "+")
    iq1 = firstofset(mpiglobal%rank, nqmt, mpiglobal%procs)
    iq2 = lastofset(mpiglobal%rank, nqmt, mpiglobal%procs)
    write(*,*) "Rank=", rank, " iq1=", iq1, " iq2=", iq2
  else
    write(unitout, '("Info(",a,"):", a)') trim(thisname),&
      & " Distributing BSE matrix, not qmt-points"
    call printline(unitout, "+")
    iq1 = 1
    iq2 = nqmt
  end if

  do iqmt = iqmti+iq1-1, iqmti+iq2-1

    write(*,*) "Rank=", rank, " iqmt=", iqmt

    vqmt(:) = input%xs%qpointset%qpoint(:, iqmt)

    call printline(unitout, "-")
    write(unitout, '("Info(",a,"):", a, i3)') trim(thisname),&
      & " Momentum tranfer list index: iqmt=", iqmt
    write(unitout, '("Info(",a,"):", a, 3f8.3)') trim(thisname),&
      & " Momentum tranfer: vqmtl=", vqmt(1:3)
    call printline(unitout, "-")

    call b_bse(iqmt)

    call printline(unitout, "-")
    write(unitout, '("Info(",a,"): Spectrum finished for iqmt=", i3)')&
      &trim(thisname), iqmt
    call printline(unitout, "-")

    if(mpiglobal%rank == 0) then
      write(6, '(a,"BSE(q) Progress:", f10.3)', advance="no")&
        & achar( 13), 100.0d0*dble(iqmt-iqmti+1)/dble(iq2-iq1+1)
      flush(6)
    end if

  end do

  if(mpiglobal%rank == 0) then
    write(6, *)
  end if

  if(iq2<0) then
    write(unitout, '("Info(",a,"): Rank= ", i3, " is idle.")')&
      & trim(thisname), mpiglobal%rank
  end if

  call barrier(callername=thisname)

end subroutine b_bselauncher
!EOC

