!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: hmlaa
! !INTERFACE:
!
!
Subroutine hmlaan (hamilton, is, ia, ngp, apwalm)
! !USES:
      Use modinput
      Use modmain
      Use modfvsystem
! !INPUT/OUTPUT PARAMETERS:
!   tapp   : .true. if the Hamiltonian is to be applied to the input vector,
!            .false. if the full matrix is to be calculated (in,logical)
!   is     : species number (in,integer)
!   ia     : atom number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   v      : input vector to which H is applied if tapp is .true., otherwise
!            not referenced (in,complex(nmatmax))
!   h      : H applied to v if tapp is .true., otherwise it is the Hamiltonian
!            matrix in packed form (inout,complex(npmatmax))
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the Hamiltonian matrix.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
!
! arguments
      Type (HermitianMatrix), Intent (Inout) :: hamilton
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Complex (8) :: x (ngp), y (ngp)
!
! local variables
      Integer :: ias, io1, io2, naa1, naa3, inonz,ireset1,ireset3
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
      Real (8) :: t1
      Complex (8) zt1, zsum, viens
      Real(8) :: alpha,a2,energyref,ta,tb
      Parameter (alpha=1d0 / 137.03599911d0)
      complex(8),allocatable::zm1(:,:),zm2(:,:),middle(:,:),zm3(:,:)
      Complex (8) zdotu
      external zdotu
    
! automatic arrays
      Complex (8) zv (ngp)
! external functions
      Real (8) :: polynom
      Complex (8) zdotc
      External polynom, zdotc

      call timesec(ta)
      ias = idxas (ia, is)
      allocate(zm3(haaijSize,ngp))
      allocate(zm2(haaijSize,ngp))
      zm3=dcmplx(0d0,0d0)
      naa3=0
      Do l3 = 0, input%groundstate%lmaxapw
         Do m3 = - l3, l3
            lm3 = idxlm (l3, m3)
            Do io2 = 1, apword (l3, is)
              naa3=naa3+1
              zm3(naa3,:)=apwalm(1:ngp, io2, lm3, ias)
            End Do
         End Do
      End Do
call timesec(tb)
!write(*,*) tb-ta
ta=tb

      zm2=zzero
      viens=dcmplx(1d0,0)
      call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                 haaijSize, &          ! M ... rows of op( A ) = rows of C
                 ngp, &           ! N ... cols of op( B ) = cols of C
                 haaijSize, &          ! K ... cols of op( A ) = rows of op( B )
                 viens, &          ! alpha
                 haaij(:,:,ias), &        ! A
                 haaijSize,&           ! LDA ... leading dimension of A
                 zm3, &           ! B
                 haaijSize, &          ! LDB ... leading dimension of B
                 viens, &          ! beta
                 zm2, &  ! C
                 haaijSize &      ! LDC ... leading dimension of C
                )
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
                 zm2, &           ! B
                 haaijSize, &          ! LDB ... leading dimension of B
                 viens, &          ! beta
                 hamilton%za, &  ! C
                 hamilton%rank &      ! LDC ... leading dimension of C
                )

       
call timesec(tb)
    deallocate(zm2,zm3)
      Return
End Subroutine
!EOC
