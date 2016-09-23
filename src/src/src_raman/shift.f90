! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created 1994, Claudia Draxl
! adapted for exciting 2014, Stefan Kontur
!
!
subroutine SHIFT(zs)
!
! shift of potential polynomial of degree 6 by zs
! the altered coefficients are returned
!
   use raman_coeff
   implicit none
   real(8) :: zs
!
   a0 =       a0       +       a1*zs    +       a2*zs**2 + &
  &           a3*zs**3 +       a4*zs**4 +       a5*zs**5 + &
  &           a6*zs**6
   a1 =       a1       +  2.d0*a2*zs    +  3.d0*a3*zs**2 + &
  &      4.d0*a4*zs**3 +  5.d0*a5*zs**4 +  6.d0*a6*zs**5
   a2 =       a2       +  3.d0*a3*zs    +  6.d0*a4*zs**2 + &
  &     10.d0*a5*zs**3 + 15.d0*a6*zs**4
   a3 =       a3       +  4.d0*a4*zs    + 10.d0*a5*zs**2 + &
  &     20.d0*a6*zs**3
   a4 =       a4       +  5.d0*a5*zs    + 15.d0*a6*zs**2
   a5 =       a5       +  6.d0*a6*zs
   a6 = a6
!
   return
end subroutine shift
!
