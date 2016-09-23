!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine symsci0 (flag, scrnh0, scrnih0, scrnisym)
      Implicit None
  ! arguments
      Integer, Intent (In) :: flag
      Complex (8), Intent (In) :: scrnh0 (3)
      Complex (8), Intent (In) :: scrnih0 (3)
      Complex (8), Intent (Out) :: scrnisym
      Select Case (flag)
      Case (0)
     ! Peter's original choice in his BSE implementation:
     ! average the screening diagonal tensor components and take
     ! inverse of this value
         scrnisym = 1.d0 / ((sum(scrnh0))/3.d0)
      Case (1)
     ! Treatment found in the BSE implementation of R. Laskowski, also
     ! mentioned in [M. Hybertsen, PRB 34, 5390 (1986), p5411]: average
     ! the diagonal tensor components of the inverse of the screening matrix
         scrnisym = sum (scrnih0) / 3.d0
      Case Default
         Write (*,*)
         Write (*, '("Error(symsci0): not a valid choice for symmetrizi&
        &ng")')
         Write (*, '(" the screened Coulomb interaction for q=0:", i8)') flag
         Write (*,*)
         Call terminate
      End Select
End Subroutine symsci0
