! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created 1994, Claudia Draxl
! adapted for exciting 1994, Stefan Kontur
!
! -----------------------------------------------------------------------------
SUBROUTINE MATRIX_TR(iidim)
! -----------------------------------------------------------------------------
!
!     determines matrix elements 
!
      use raman_trmat
      implicit none
      integer :: iidim,n,kk,kj,j
!
      n = iidim/2
!     DO KK = 1,iidim                                 ! set zero matrices T and R
!        DO KJ = 1,iidim
!           T(KK,KJ)=0.d0
!           R(KK,KJ)=0.d0
!        enddo
!     enddo
      T = 0.d0
      R = 0.d0
!
!  compute values of sparse matrices T and R.
!  only four elements arond diagonal have to be
!  computed, due to interval limits for the basis
!  functions u_i
!
      do j = 1,3
         T(1,j) =         S(1,2,j+1)
         R(1,j) =         O(2,j+1)
      enddo
      T(2,1) =            S(1,3,2)
      T(3,1) =            S(1,4,2)
      R(2,1) =            O(3,2)
      R(3,1) =            O(4,2)
!
      do j = 2,2*n-4,2
         T(j,j) =         S(j/2,3,3) + S(j/2+1,1,1)   
         T(j+1,j+1) =     S(j/2,4,4) + S(j/2+1,2,2)   
         T(j,j+1) =       S(j/2,3,4) + S(j/2+1,1,2)   
         T(j+1,j) =       S(j/2,4,3) + S(j/2+1,2,1)
!
         R(j,j) =         O(3,3) + O(1,1)
         R(j+1,j+1) =     O(4,4) + O(2,2)
         R(j,j+1) =       O(3,4) + O(1,2)
         R(j+1,j) =       O(4,3) + O(2,1)
!   
!
         T(j,j+2) =       S(j/2+1,1,3)   
         T(j,j+3) =       S(j/2+1,1,4)
         T(j+1,j+2) =     S(j/2+1,2,3)
         T(j+1,j+3) =     S(j/2+1,2,4)
         T(j+2,j) =       S(j/2+1,3,1)
         T(j+2,j+1) =     S(j/2+1,3,2)
         T(j+3,j) =       S(j/2+1,4,1)
         T(j+3,j+1) =     S(j/2+1,4,2)
! 
         R(j,j+2) =       O(1,3)
         R(j,j+3) =       O(1,4)
         R(j+1,j+2) =     O(2,3)
         R(j+1,j+3) =     O(2,4)
         R(j+2,j) =       O(3,1)
         R(j+2,j+1) =     O(3,2)
         R(j+3,j) =       O(4,1)
         R(j+3,j+1) =     O(4,2)
      enddo
!
      T(2*n-2,2*n-2) = S(n-1,3,3) + S(n,1,1)
      T(2*n-1,2*n-1) = S(n-1,4,4) + S(n,2,2)
      T(2*n-2,2*n-1) = S(n-1,3,4) + S(n,1,2)
      T(2*n-1,2*n-2) = S(n-1,4,3) + S(n,2,1)
      T(2*n-2,2*n) =   S(n,1,4)
      T(2*n-1,2*n) =   S(n,2,4)
      T(2*n,2*n-2) =   S(n,4,1)
      T(2*n,2*n-1) =   S(n,4,2)
      T(2*n,2*n) =     S(n,4,4)
!
      R(2*n-2,2*n-2) = O(3,3) + O(1,1)
      R(2*n-1,2*n-1) = O(4,4) + O(2,2)
      R(2*n-2,2*n-1) = O(3,4) + O(1,2)
      R(2*n-1,2*n-2) = O(4,3) + O(2,1)
      R(2*n-2,2*n) =   O(1,4)
      R(2*n-1,2*n) =   O(2,4)
      R(2*n,2*n-2) =   O(4,1)
      R(2*n,2*n-1) =   O(4,2)
      R(2*n,2*n) =     O(4,4)
!

Contains
!
      real(8) function s(i,j,k)
!
!   determines matrix elements S(j,k) for interval i
!        
      use raman_fij       
      use raman_inter
      implicit none
      integer :: i,j,k,l  
      double precision :: a(0:4)
      double precision :: Shelp
!        
      a(0) = b0(i)                  ! b0(i)...b4(i) are coefficients from interpolation of potential for interval i
      a(1) = b1(i)
      a(2) = b2(i)
      a(3) = b3(i)
      a(4) = b4(i)
!  S(j,k) = -h_bar/2M 1/h^2 DFI(k,j) + Sum_(l=0)^4 a(k) FI(k,j,l)
      Shelp = -sfact * DFI(k,j) / h**2
      do l = 0, 4
         Shelp = Shelp + a(l)*FI(k,j,l)
      enddo
      S = Shelp
      return
      end function s
!
!
      real(8) function O(j,k)
!
!   determines overlap matrix elements O(j,k)
!
      use raman_inter, only: h
      use raman_fij
      implicit none
      integer :: j,k
!
      O = FI(k,j,0)              ! O(j,k) = FI(k,j,0)
!
      return
      end function o
!
!
end subroutine matrix_tr
!
