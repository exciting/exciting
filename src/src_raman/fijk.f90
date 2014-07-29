! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created 1994, Claudia Draxl
! adapted for exciting 2014, Stefan Kontur
!
! -----------------------------------------------------------------------------
SUBROUTINE FIJK
! -----------------------------------------------------------------------------
      use raman_fij
      implicit none
      integer :: ip(4,4),idp(2,4)
      integer :: i,j,k
!
!     calculation of following integrals:
!
!     FI(i,j,k)=Int(P(i)*P(j)*x**k;0,1)
!     DFI(i,j)=Int(d2P(i)/dx**2 * P(j); 0,1)
!
!
!     Coefficients of cubic polynomials
!
      ip(1,1) =  1            ! P1(xi) = ip(1,1)xi^0 + ip(2,1)xi^1 + ip(3,1)xi^2 + ip(4,1)xi^3
      ip(2,1) =  0            ! P2(xi) = ip(1,2)xi^0 + ip(2,2)xi^1 + ip(3,2)xi^2 + ip(4,2)x1^3
      ip(3,1) = -3            ! etc.
      ip(4,1) =  2
      ip(1,2) =  0
      ip(2,2) =  1
      ip(3,2) = -2
      ip(4,2) =  1
      ip(1,3) =  0
      ip(2,3) =  0
      ip(3,3) =  3
      ip(4,3) = -2
      ip(1,4) =  0
      ip(2,4) =  0
      ip(3,4) = -1
      ip(4,4) =  1
!
      idp(1,1) =  -6          ! P1''(xi) = idp(1,1)xi^0 + idp(2,1)xi^1
      idp(2,1) =  12          ! P2''(xi) = idp(1,2)xi^0 + idp(2,2)xi^1
      idp(1,2) =  -4          ! P3''(xi) = idp(1,3)xi^0 + idp(2,3)xi^1
      idp(2,2) =   6          ! etc.
      idp(1,3) =   6
      idp(2,3) = -12
      idp(1,4) =  -2
      idp(2,4) =   6
!
!  Precompute the two integrals
!
      do i = 1,4
         do j = 1,4
            do k = 0,6
               FI(i,j,k) = dble(ip(1,i)*ip(1,j)) / dble(k+1) + &
   &                       dble(ip(1,i)*ip(2,j) + ip(2,i)*ip(1,j)) / dble(k+2) + &
   &                       dble(ip(1,i)*ip(3,j) + ip(2,i)*ip(2,j) + ip(3,i)*ip(1,j)) / dble(k+3) + &
   &                       dble(ip(1,i)*ip(4,j) + ip(2,i)*ip(3,j) + ip(3,i)*ip(2,j) + ip(4,i)*ip(1,j)) / dble(k+4) + &
   &                                         dble(ip(2,i)*ip(4,j) + ip(3,i)*ip(3,j) + ip(4,i)*ip(2,j)) / dble(k+5) + &
   &                                                           dble(ip(3,i)*ip(4,j) + ip(4,i)*ip(3,j)) / dble(k+6) + &
   &                                                                             dble(ip(4,i)*ip(4,j)) / dble(k+7)
            enddo
            DFI(i,j) = dble(idp(1,i)*ip(1,j)) + &
   &                   dble(idp(1,i)*ip(2,j) + idp(2,i)*ip(1,j)) / 2.d0 + &
   &                   dble(idp(1,i)*ip(3,j) + idp(2,i)*ip(2,j)) / 3.d0 + &
   &                   dble(idp(1,i)*ip(4,j) + idp(2,i)*ip(3,j)) / 4.d0 + &
   &                                      dble(idp(2,i)*ip(4,j)) / 5.d0
         enddo
      enddo
!
      return
      end
!
