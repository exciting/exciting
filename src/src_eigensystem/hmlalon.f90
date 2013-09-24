!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine hmlalon (hamilton, is, ia, ngp, apwalm)
      Use modmain
      Use modinput
      Use modfvsystem
      Implicit None
! arguments
      Type (hermitianMatrix), Intent (Inout) :: hamilton
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
!
!
! local variables
      Integer :: ias, io, ilo, i, j, k, naa3, nalo1, if1, if3, is1
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3, inonz, ireset3, maxnlo
      Complex (8) zsum, zt1, viens
      Complex (8), allocatable :: zm3(:,:),integrals(:,:)

      ias = idxas (ia, is)
      viens= dcmplx(1d0,0d0)
      allocate(zm3(haaijSize,ngp))
      zm3=dcmplx(0d0,0d0)
      naa3=0
      Do l3 = 0, input%groundstate%lmaxapw
         Do m3 = - l3, l3
            lm3 = idxlm (l3, m3)
            Do io = 1, apword (l3, is)
              naa3=naa3+1
              zm3(naa3,:)=apwalm(1:ngp, io, lm3, ias)
            End Do
         End Do
      End Do

      maxnlo=size(haloij,1)
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                 nalo1, &          ! M ... rows of op( A ) = rows of C
                 ngp, &           ! N ... cols of op( B ) = cols of C
                 haaijSize, &          ! K ... cols of op( A ) = rows of op( B )
                 viens, &          ! alpha
                 haloij(:,:,ias), &           ! A
                 maxnlo,&           ! LDA ... leading dimension of A
                 zm3, &           ! B
                 haaijSize, &          ! LDB ... leading dimension of B
                 viens, &          ! beta
                 hamilton%za(ngp+idxlo (lm1, 1, ias),1), &  ! C
                 hamilton%rank &      ! LDC ... leading dimension of C
                )

      l3 = lorbl (1, is)
      lm3 = idxlm (l3, -l3)
      do i=idxlo (lm3,1,ias)+ngp, nalo1+ngp+idxlo (lm3,1,ias)-1
         hamilton%za(1:ngp,i)=conjg(hamilton%za(i,1:ngp))
      enddo
 deallocate(zm3)
      Return
End Subroutine
