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
      Integer :: ias, io, ilo, i, j, k, naa3, nalo1, if1, if3
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3, inonz, ireset3
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

     ilo=nlorb (is)
     l1 = lorbl (ilo, is)
     lm1 = idxlm (l1, l1)
     l3 = lorbl (1, is)
     lm3 = idxlm (l3, -l3)
     nalo1=idxlo (lm1, ilo, ias)- idxlo (lm3, 1, ias)+1
     allocate(integrals(nalo1,haaijSize))
     integrals=dcmplx(0d0,0d0)


     if1=0 
      Do ilo = 1, nlorb (is)
         l1 = lorbl (ilo, is)
         inonz=gntnonzlindex(l1)
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            i = ngp + idxlo (lm1, ilo, ias)
            if1=if1+1
            if3=0
            Do l3 = 0, input%groundstate%lmaxmat
               Do m3 = - l3, l3
                  lm3 = idxlm (l3, m3)
                  ireset3=inonz
                  Do io = 1, apword (l3, is)
                     if3=if3+1
                     zsum = 0.d0
                     do while ((gntnonzlm3(inonz).eq.lm3).and.(gntnonzlm1(inonz).eq.lm1))
                       zsum=zsum+gntnonz(inonz)*hloa(gntnonzlm2(inonz),io, l3, ilo, ias)
                       inonz=inonz+1
                     enddo
                     if (io.ne.apword(l3,is)) inonz=ireset3
                  End Do
               End Do
            End Do
         End Do
      End Do

     l1 = lorbl (1, is)
     lm1 = idxlm (l1, -l1)

      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                 nalo1, &          ! M ... rows of op( A ) = rows of C
                 ngp, &           ! N ... cols of op( B ) = cols of C
                 haaijSize, &          ! K ... cols of op( A ) = rows of op( B )
                 viens, &          ! alpha
                 integrals, &           ! A
                 nalo1,&           ! LDA ... leading dimension of A
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
 deallocate(integrals,zm3)
      Return
End Subroutine
