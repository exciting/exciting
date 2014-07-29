! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created 1994, Claudia Draxl
! adapted for exciting 2014, Stefan Kontur
!
!
subroutine POTX6(maxp)
!
!  assumes coefficients for a polynomial of degree 6 and computes values at knots xpot(i)
!
   use raman_inter
   use raman_coeff
   use modinput
   implicit none
   integer :: i,maxp
   real(8) :: h1,x
!
   h1 = (input%properties%raman%xmax - input%properties%raman%xmin)/dble(maxp - 1)
   do i = 1,maxp
      xpot(i) = input%properties%raman%xmin + dble(i - 1)*h1
      x = xpot(i)
      pot(i) = A0 + A1*x + A2*x**2 + A3*x**3 + A4*x**4 + A5*x**5 + A6*x**6
!        write(77,*) x,pot(i)
   enddo
   return
!
end subroutine potx6
!
