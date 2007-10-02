
! Copyright (C) 2002-2006 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeinfotd
! !INTERFACE:
subroutine writeinfotd
! !USES:
  use modmain
  use modtddft
! !INPUT/OUTPUT PARAMETERS:
!   fnum : unit specifier for INFOTD.OUT file (in,integer)
! !DESCRIPTION:
!   Outputs basic information about the run to the file {\tt INFOTD.OUT}.
!   Does not close the file afterwards.
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  character(*), parameter :: thisnam='writeinfotd'
  
  ! write angular momenta
  write(unitout,'(a,2i8)') 'Info('//thisnam//'): maximum angular momentum &
       &for Hamiltonian matrix elements     :', lmaxmat, lmmaxmat
  write(unitout,'(a,2i8)') 'Info('//thisnam//'): maximum angular momentum &
       &for the inner part of the muffin-tin:', lmaxinr, lmmaxinr
  write(unitout,'(a,2i8)') 'Info('//thisnam//'): maximum angular momentum &
       &for density and potential           :', lmaxvr, lmmaxvr
  write(unitout,'(a,2i8)') 'Info('//thisnam//'): maximum angular momentum &
       &for APW functions                   :', lmaxapw, lmmaxapw
  write(unitout,'(a,2i8)') 'Info('//thisnam//'): maximum angular momentum &
       &for APW functions (matr. el.)       :', lmaxapwtd, lmmaxapwtd
  write(unitout,'(a,2i8)') 'Info('//thisnam//'): maximum angular momentum &
       &for local orbitals                  :', lolmax, lolmmax
  write(unitout,'(a,2i8)') 'Info('//thisnam//'): maximum angular momentum &
       &for matrix elements of exponential  :', lmaxemat, lmmaxemat
  write(unitout,'(a,2i8)') 'Info('//thisnam//'): overall maximum angular  &
       &momentum                            :', lmaxmax, lmmaxmax
  ! information about Gamma point
  if (tq1gamma) then
     write(unitout,'(a)') 'Info('//thisnam//'): first q-point is the &
          &Gamma point'
  end if

end subroutine writeinfotd
!EOC

