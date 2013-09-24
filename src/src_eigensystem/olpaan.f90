!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine olpaan (overlap, is, ia, ngp, apwalm)
      Use modmain
      Use modinput
      Use modfvsystem
      Implicit None
! arguments
      Type (hermitianmatrix), Intent (Inout) :: overlap
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Complex (8) :: x (ngp), y (ngp)
!
! local variables
      Integer :: ias, l, m, lm, io, naa3
      real(8) :: ta, tb
      Complex (8) :: viens
      Complex (8), allocatable :: zm3(:,:)
! external functions
      Complex (8) zdotu
      External zdotu
      call timesec(ta)
      viens=dcmplx(1d0,0d0)
      ias = idxas (ia, is)
      allocate(zm3(haaijSize,ngp))
      zm3=dcmplx(0d0,0d0)
      naa3=0
      Do l = 0, input%groundstate%lmaxapw
         Do m = - l, l
            lm = idxlm (l, m)
            Do io = 1, apword (l, is)
              naa3=naa3+1
              zm3(naa3,:)=apwalm(1:ngp, io, lm, ias)
            End Do
         End Do
      End Do
call timesec(tb)
!write(*,*) tb-ta
ta=tb
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                 ngp, &          ! M ... rows of op( A ) = rows of C
                 ngp, &           ! N ... cols of op( B ) = cols of C
                 haaijSize, &          ! K ... cols of op( A ) = rows of op( B )
                 viens, &          ! alpha
                 zm3, &           ! A
                 haaijSize,&           ! LDA ... leading dimension of A
                 zm3, &           ! B
                 haaijSize, &          ! LDB ... leading dimension of B
                 viens, &          ! beta
                 overlap%za, &  ! C
                 overlap%rank &      ! LDC ... leading dimension of C
                )


call timesec(tb)
    deallocate(zm3)


!      Do l = 0, input%groundstate%lmaxmat
!         Do m = - l, l
!            lm = idxlm (l, m)
!            Do io = 1, apword (l, is)
!               x = conjg (apwalm(1:ngp, io, lm, ias))
!               y = conjg (apwalm(1:ngp, io, lm, ias))
!               Call Hermitianmatrix_rank2update (overlap, ngp, zhalf, &
!              & x, y)
!            End Do
!         End Do
!      End Do
      Return
End Subroutine
