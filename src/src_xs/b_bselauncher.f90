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

  integer(4) :: iqmt, iqmti, iqmtf
  real(8) :: vqmt(3)

  ! Testing 
  iqmti = 1
  iqmtf = size(input%xs%qpointset%qpoint, 2)
  if(input%xs%bse%iqmt /= -1) then 
    iqmti=input%xs%bse%iqmt
    iqmtf=input%xs%bse%iqmt
  end if

  if(mpiglobal%rank == 0) then
    !write(unitout, *)
    call printline(unitout, "+")
    write(unitout, '("Info(b_bselauncher):", a)')&
      & " Setting up and diagonalizing BSE Hamiltonian."
    write(unitout, '("Info(b_bselauncher):", a, i3, a, i3)')&
      & " Using momentum transfer vectors from list : ", iqmti, " to", iqmtf
    call printline(unitout, "+")
    !write(unitout,*)
  end if

  do iqmt = iqmti, iqmtf

    vqmt(:) = input%xs%qpointset%qpoint(:, iqmt)

    if(mpiglobal%rank == 0) then
      !write(unitout, *)
      call printline(unitout, "-")
      write(unitout, '("Info(b_bselauncher):", a, i3)')&
        & " Momentum tranfer list index: iqmt=", iqmt
      write(unitout, '("Info(b_bselauncher):", a, 3f8.3)')&
        & " Momentum tranfer: vqmtl=", vqmt(1:3)
      call printline(unitout, "-")
      !write(unitout,*)
    end if

    call b_bse(iqmt)

    if(mpiglobal%rank == 0) then
      call printline(unitout, "-")
      write(unitout, '("Info(b_bselauncher): Spectrum finished for iqmt=", i3)') iqmt
      call printline(unitout, "-")
    end if

  end do

end subroutine b_bselauncher
!EOC

