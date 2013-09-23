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
      Complex (8) zt1, zsum, integrals(lmmaxapw),viens
! It should have been integrals(lmmaxvr), but the hard-coded condition lmmaxvr<=lmmaxapw permits it
      Real(8) :: alpha,a2,energyref,ta,tb
      Parameter (alpha=1d0 / 137.03599911d0)
      complex(8),allocatable::zm1(:,:),zm2(:,:),middle(:,:),zm3(:,:)
      logical :: test1
      Complex (8) zdotu
      external zdotu
    
! automatic arrays
      Complex (8) zv (ngp)
! external functions
      Real (8) :: polynom
      Complex (8) zdotc
      External polynom, zdotc

      call timesec(ta)
      energyref=input%groundstate%energyref
      if (input%groundstate%ValenceRelativity.eq."scalar") then
         a2=0.5d0*alpha**2
      else
         a2=0d0
      endif
 ias = idxas (ia, is)

      naa1=0
      naa3=0
      Do l1 = 0, input%groundstate%lmaxmat
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do io1 = 1, apword (l1, is)
              naa1=naa1+1
            End Do
         End Do
      End Do
      Do l3 = 0, input%groundstate%lmaxapw
         Do m3 = - l3, l3
            lm3 = idxlm (l3, m3)
            Do io2 = 1, apword (l3, is)
              naa3=naa3+1
            End Do
         End Do
      End Do
      allocate(middle(naa1,naa3))
      allocate(zm1(ngp,naa1))
      allocate(zm3(naa3,ngp))
      allocate(zm2(naa1,ngp))
      naa1=0
      naa3=0
      Do l1 = 0, input%groundstate%lmaxmat
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do io1 = 1, apword (l1, is)
              naa1=naa1+1
              zm1(:,naa1)=conjg(apwalm(1:ngp, io1, lm1, ias))
            End Do
         End Do
      End Do
      Do l3 = 0, input%groundstate%lmaxapw
         Do m3 = - l3, l3
            lm3 = idxlm (l3, m3)
            Do io2 = 1, apword (l3, is)
              naa3=naa3+1
              zm3(naa3,:)=apwalm(1:ngp, io2, lm3, ias)
            End Do
         End Do
      End Do
      middle=zzero
      
      inonz=1
      naa1=0
call timesec(tb)
!write(*,*) tb-ta
ta=tb
      Do l1 = 0, input%groundstate%lmaxmat
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            ireset1=inonz
            Do io1 = 1, apword (l1, is)
!               zv (:) = 0.d0
               naa1=naa1+1
               naa3=0
               Do l3 = 0, input%groundstate%lmaxapw
                  Do m3 = - l3, l3
                     lm3 = idxlm (l3, m3)
                        ireset3=inonz
                        Do io2 = 1, apword (l3, is)
                          naa3=naa3+1
                          zsum = 0.d0
                          do while ((gntnonzlm3(inonz).eq.lm3).and.(gntnonzlm1(inonz).eq.lm1))
                            zsum=zsum+gntnonz(inonz)*haa(gntnonzlm2(inonz), io2, l3, io1, l1, ias)
                            inonz=inonz+1
                          enddo
!                          integrals(1:lmmaxvr)=dcmplx(haa (1:lmmaxvr, io2, l3, io1, l1, ias),0d0)
!                          zsum=zdotu(lmmaxvr,gntryy (1:lmmaxvr, lm3, lm1),1,integrals(1:lmmaxvr),1)
                          middle(naa1,naa3)=zsum
                          if (io2.ne.apword (l3, is)) inonz=ireset3
                        End Do
                  End Do
               End Do
              if (io1.ne.apword (l1, is)) inonz=ireset1
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
                 naa1, &          ! M ... rows of op( A ) = rows of C
                 ngp, &           ! N ... cols of op( B ) = cols of C
                 naa3, &          ! K ... cols of op( A ) = rows of op( B )
                 viens, &          ! alpha
                 middle, &        ! A
                 naa1,&           ! LDA ... leading dimension of A
                 zm3, &           ! B
                 naa3, &          ! LDB ... leading dimension of B
                 viens, &          ! beta
                 zm2, &  ! C
                 naa1 &      ! LDC ... leading dimension of C
                )
call timesec(tb)
!write(*,*) tb-ta
ta=tb 
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                 ngp, &          ! M ... rows of op( A ) = rows of C
                 ngp, &           ! N ... cols of op( B ) = cols of C
                 naa1, &          ! K ... cols of op( A ) = rows of op( B )
                 viens, &          ! alpha
                 zm1, &           ! A
                 ngp,&           ! LDA ... leading dimension of A
                 zm2, &           ! B
                 naa1, &          ! LDB ... leading dimension of B
                 viens, &          ! beta
                 hamilton%za, &  ! C
                 hamilton%rank &      ! LDC ... leading dimension of C
                )

       
call timesec(tb)
!write(*,*) tb-ta
!stop
     if (.not.input%groundstate%SymmetricKineticEnergy) then
! kinetic surface contribution
      t1 = 0.25d0 * rmt (is) ** 2
      Do l1 = 0, input%groundstate%lmaxmat
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do io1 = 1, apword (l1, is)
               Do io2 = 1, apword (l1, is)
                  zt1 = t1 * apwfr (nrmt(is), 1, io1, l1, ias) * apwdfr (io2, l1, ias) *1d0/(1d0+(energyref-veffmt(1, nrmt(is), ias)*y00)*a2)
                  x = conjg (apwalm(1:ngp, io1, lm1, ias))
                  y = conjg (apwalm(1:ngp, io2, lm1, ias))
                  Call Hermitianmatrix_rank2update (hamilton, ngp, zt1, &
                 & x, y)
               End Do
            End Do
         End Do
      End Do
    endif
if (test1) then
    deallocate(zm1,zm2,zm3,middle)
else
    deallocate(zm1,zm2)
endif
!    write(*,*) hamilton%za(2,2)
      Return
End Subroutine
!EOC
