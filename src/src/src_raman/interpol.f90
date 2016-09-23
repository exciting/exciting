! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
!
subroutine INTERPOL
!
   use modinput
   use raman_inter, only : xpot,pot
   use raman_inter, only : b0,b1,b2,b3,b4
   implicit none
   integer :: i
   real(8) :: vx0,vx14,vx12,vx34,vx1
!
!    4th order interpolation for given knots of the potential
!    V = b0 + b1 xi + b2 xi^2 + b3 xi^3 + b4 xi^4
!    in the interval (0,1) knots are at 0, 1/4, 1/2, 3/4 and 1 (=0 of the next interval)
!    V(0) = b0
!    V(1/4) = b0 + 1/4 b1 + 1/16 b2 + 1/64 b3 + 1/256 b4
!    V(1/2) = b0 + 1/2 b1 + 1/4 b2 + 1/8 b3 + 1/16 b4
!    V(3/4) = b0 + 3/4 b1 + 9/16 b2 * 27/64 b3 + 81/256 b4
!    V(1) = b0 + b1 + b2 + b3 + b4
!    from where one gets:
!
      do i = 1, input%properties%raman%ninter
         b0(i) =  pot(4*i-3)               ! b0 = V(0)
         b1(i) = -pot(4*i-3)* 25.d0/3.d0 & ! b1 = -25/3 V(0) + 16 V(1/4) - 12 V(1/2) + 16/3 V(3/4) - V(1)
     &           +pot(4*i-2)* 16.d0 &
     &           -pot(4*i-1)* 12.d0 &
     &           +pot(4*i)  * 16.d0/3.d0 &
     &           -pot(4*i+1)
         b2(i) =  pot(4*i-3)* 70.d0/3.d0 & ! b2 = 70/3 V(0) - 208/3 V(1/4) + 76 V(1/2) - 112/3 V(3/4) + 22/3 V(1)
     &           -pot(4*i-2)*208.d0/3.d0 &
     &           +pot(4*i-1)* 76.d0 &
     &           -pot(4*i)  *112.d0/3.d0 &
     &           +pot(4*i+1)* 22.d0/3.d0
         b3(i) = -pot(4*i-3)* 80.d0/3.d0 & ! b3 = -80/3 V(0) + 96 V(1/4) - 128 V(1/2) + 224/3 V(3/4) - 16 V(1)
     &           +pot(4*i-2)* 96.d0 &
     &           -pot(4*i-1)*128.d0 &
     &           +pot(4*i)  *224.d0/3.d0 &
     &           -pot(4*i+1)* 16.d0
         b4(i) =  pot(4*i-3)* 32.d0/3.d0 & ! b4 = 32/3 V(0) - 128/3 V(1/4) + 64 V(1/2) - 128/3 V(3/4) + 32/3 V(1)
     &           -pot(4*i-2)*128.d0/3.d0 &
     &           +pot(4*i-1)* 64.d0 &
     &           -pot(4*i)  *128.d0/3.d0 &
     &           +pot(4*i+1)* 32.d0/3.d0
!
!     coefficients for [0,1]
!  
         vx0 = b0(i)
         vx14 = b0(i) + 1.d0/4.d0*b1(i) + 1.d0/16.d0*b2(i) + &
     &                  1.d0/64.d0*b3(i) + 1.d0/256.d0*b4(i)
         vx12 = b0(i) + 1.d0/2.d0*b1(i) + 1.d0/4.d0*b2(i) + &
     &                  1.d0/8.d0*b3(i) + 1.d0/16.d0*b4(i)
         vx34 = b0(i) + 3.d0/4.d0*b1(i) + 9.d0/16.d0*b2(i) + &
     &                  27.d0/64.d0*b3(i) + 81.d0/256.d0*b4(i)
         if (i.eq.1) then
            write(77,99) xpot(4*i-3),pot(4*i-3),vx0
         else
            vx1 = b0(i-1) + b1(i-1) + b2(i-1) + b3(i-1) + b4(i-1)
            write(77,100) xpot(4*i-3),pot(4*i-3),vx0,vx1
         endif
         write(77,99) xpot(4*i-2),pot(4*i-2),vx14
         write(77,99) xpot(4*i-1),pot(4*i-1),vx12
         write(77,99) xpot(4*i),pot(4*i),vx34
      enddo
      return
  99  format(3f14.8)
 100  format(4f14.8)
end subroutine interpol
!
