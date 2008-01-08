
! Copyright (C) 2002-2006 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeangmom
! !INTERFACE:
subroutine writeangmom(un)
! !USES:
  use modmain
  use modxs
! !INPUT/OUTPUT PARAMETERS:
! !DESCRIPTION:
!   Outputs information about the angular momentum cutoffs to the file
!   {\tt INFO_XS.OUT}.
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: un
  ! local variables
  character(*), parameter :: thisnam='writeangmom'
  if (calledxs.eq.1) then
     ! write angular momenta
     write(un,'(a)') 'Info('//thisnam//'): angular momenta:'
     write(un,'(a,2i8)') ' density and potential        : ',&
          lmaxvr, (lmaxvr+1)**2
     write(un,'(a,2i8)') ' APW functions                : ',&
          lmaxapw, (lmaxapw+1)**2
     write(un,'(a,2i8)') ' local orbitals               : ',&
          lolmax, (lolmax+1)**2
     write(un,'(a,2i8)') ' inner part of muffin-tin     : ',&
          lmaxinr, (lmaxinr+1)**2
     write(un,'(a,2i8)') ' Hamiltonian. matr. el.       : ',&
          lmaxmat, (lmaxmat+1)**2
     write(un,'(a,2i8)') ' PW matr. el.                 : ',&
          lmaxemat, (lmaxemat+1)**2
     write(un,'(a,2i8)') ' APW functions (PW matr. el.) : ',&
          lmaxapwtd, (lmaxapwtd+1)**2
     write(un,'(a,2i8)') ' overall                      : ',&
          lmaxmax, (lmaxmax+1)**2
  end if
end subroutine writeangmom
!EOC

