! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created 1994, Claudia Draxl 
! adapted for exciting 2014, Stefan Kontur
!
!
subroutine EPSSHIFT(zs)
! 
! shift of Maclaurin series by zs up to n=6
! derivatives are updated by shifted ones
use raman_coeff
implicit none
real(8), intent(in) :: zs
integer :: oct1, oct2
do oct1 = 1, 3
   do oct2 = 1, 3
      if (oct1 .ne. oct2 .and. .not. offdiag) cycle
!
      deq (oct1, oct2) = deq (oct1, oct2)                 + &
                     &   d2eq(oct1, oct2)* zs             + &
                     &   d3eq(oct1, oct2)* zs**2 /   2.d0 + &
                     &   d4eq(oct1, oct2)* zs**3 /   6.d0 + &
                     &   d5eq(oct1, oct2)* zs**4 /  24.d0 + &
                     &   d6eq(oct1, oct2)* zs**5 / 120.d0
      d2eq(oct1, oct2) = d2eq(oct1, oct2)                 + &
                     &   d3eq(oct1, oct2)* zs             + &
                     &   d4eq(oct1, oct2)* zs**2 /   2.d0 + &
                     &   d5eq(oct1, oct2)* zs**3 /   6.d0 + &
                     &   d6eq(oct1, oct2)* zs**4 /  24.d0
      d3eq(oct1, oct2) = d3eq(oct1, oct2)                 + &
                     &   d4eq(oct1, oct2)* zs             + &
                     &   d5eq(oct1, oct2)* zs**2 /   2.d0 + &
                     &   d6eq(oct1, oct2)* zs**3 /   6.d0
      d4eq(oct1, oct2) = d4eq(oct1, oct2)                 + &
                     &   d5eq(oct1, oct2)* zs             + &
                     &   d6eq(oct1, oct2)* zs**2 /   2.d0
      d5eq(oct1, oct2) = d5eq(oct1, oct2)                 + &
                     &   d6eq(oct1, oct2)* zs
      d6eq(oct1, oct2) = d6eq(oct1, oct2)
   enddo
enddo
!
!deqx  =  deqx + d2eqx*zs + d3eqx*zs**2/2.d0 + d4eqx*zs**3/6.0d0 + d5eqx*zs**4/24.D0 + d6eqx*zs**5/120.d0
!deqy  =  deqy + d2eqy*zs + d3eqy*zs**2/2.d0 + d4eqy*zs**3/6.0d0 + d5eqy*zs**4/24.d0 + d6eqy*zs**5/120.d0
!deqz  =  deqz + d2eqz*zs + d3eqz*zs**2/2.d0 + d4eqz*zs**3/6.0d0 + d5eqz*zs**4/24.d9 + d6eqz*zs**5/120.d0
!d2eqx =         d2eqx    + d3eqx*zs         + d4eqx*zs**2/2.0d0 + d5eqx*zs**3/ 6.d0 + d6eqx*zs**4/ 24.d0
!d2eqy =         d2eqy    + d3eqy*zs         + d4eqy*zs**2/2.0d0 + d5eqy*zs**3/ 6.d0 + d6eqy*zs**4/ 24.d0
!d2eqz =         d2eqz    + d3eqz*zs         + d4eqz*zs**2/2.0d0 + d5eqz*zs**3/ 6.d0 + d6eqz*zs**4/ 24.d0
!d3eqx =                    d3eqx            + d4eqx*zs          + d5eqx*zs**2/ 2.d0 + d6eqx*zs**3/  6.d0
!d3eqy =                    d3eqy            + d4eqy*zs          + d5eqy*zs**2/ 2.d0 + d6eqy*zs**3/  6.d0
!d3eqz =                    d3eqz            + d4eqz*zs          + d5eqz*zs**2/ 2.d0 + d6eqz*zs**3/  6.d0
!d4eqx =                                       d4eqx             + d5eqx*zs          + d6eqx*zs**2/  2.d0
!d4eqy =                                       d4eqy             + d5eqy*zs          + d6eqy*zs**2/  2.d0
!d4eqz =                                       d4eqz             + d5eqz*zs          + d6eqz*zs**2/  2.d0
!d5eqx =                                                           d5eqx             + d6eqx*zs
!d5eqy =                                                           d5eqy             + d6eqy*zs
!d5eqz =                                                           d5eqz             + d6eqz*zs
!d6eqx =                                                                               d6eqx
!d6eqy =                                                                               d6eqy
!d6eqz =                                                                               d6eqz
!
return
end
!
!
