!
!
!
! Copyright (C) 2002-2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine hamiltonsetup (system, ngp, apwalm, igpig, vgpc)
      Use modfvsystem
      Use modinput
      Use mod_eigensystem
      Use mod_atoms
      Use mod_timing
      Use mod_muffin_tin
      Use mod_APW_LO
      Use mod_gkvector
!
      Implicit None
      Type (evsystem) :: system
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, natmtot)
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Integer :: n
      Character (256) :: prefix
!local variables
      Integer, Save :: ikc
      Real (8), Save :: cputot
      Real (8) :: cpuaa, cpualo, cpulolo, cpui, cpu00, cpu01,ts0,ts1
      Integer :: i, is, ia, maxaa
      Complex (8) v(1),viens
      Real (8) :: cpu0, cpu1
      Real (8) :: threshold
      Complex (8), allocatable :: apwi(:,:),zm(:,:),apwi2(:,:)
      integer if1,if3,l3,m3,lm3,io1,io2,ias,maxnlo,ilo,j1,j2,j3,lm1,lm2,j,io,l,ilo1,ilo2,l1

!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!
!
!
      Call timesec (cpu0)
! set the matrices to zero
!
! muffin-tin contributions
      maxaa=mt_hscf%maxaa
      maxnlo=mt_hscf%maxnlo

      allocate(apwi(maxaa,ngp))
      allocate(apwi2(ngp,maxaa))
      allocate(zm(maxaa,ngp))
      
      Do is = 1, nspecies
        Do ia = 1, natoms (is)
!--Hamiltonian--
! APW-APW part 
          Call timesec (ts0)
          ias = idxas (ia, is)
          apwi=dcmplx(0d0,0d0)
          apwi2=dcmplx(0d0,0d0)
          if3=0
          Do l3 = 0, input%groundstate%lmaxmat
            Do m3 = - l3, l3
            lm3 = idxlm (l3, m3)
              Do io2 = 1, apword (l3, is)
                if3=if3+1
                apwi(if3,:)=apwalm(1:ngp, io2, lm3, ias)
 !               apwi2(:,if3)=conjg(apwalm(1:ngp, io2, lm3, ias))
              End Do
            End Do
          End Do
          zm=zzero
          viens=dcmplx(1d0,0)

          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      maxaa, &          ! M ... rows of op( A ) = rows of C
                      ngp, &           ! N ... cols of op( B ) = cols of C
                      maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      viens, &          ! alpha
                      mt_hscf%main%aa(1,1,ias), &
                      maxaa,&           ! LDA ... leading dimension of A
                      apwi, &           ! B
                      maxaa, &          ! LDB ... leading dimension of B
                      viens, &          ! beta
                      zm, &  ! C
                      maxaa &      ! LDC ... leading dimension of C
                      )
if (.true.) then
          call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      ngp, &          ! M ... rows of op( A ) = rows of C
                      ngp, &           ! N ... cols of op( B ) = cols of C
                      maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      viens, &          ! alpha
                      apwi, &           ! A
                      maxaa,&           ! LDA ... leading dimension of A
                      zm, &           ! B
                      maxaa, &          ! LDB ... leading dimension of B
                      viens, &          ! beta
                      system%hamilton%za, &  ! C
                      system%hamilton%rank &      ! LDC ... leading dimension of C
                     )
else
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      ngp, &          ! M ... rows of op( A ) = rows of C
                      ngp, &           ! N ... cols of op( B ) = cols of C
                      maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      viens, &          ! alpha
                      apwi2, &           ! A
                      ngp,&           ! LDA ... leading dimension of A
                      zm, &           ! B
                      maxaa, &          ! LDB ... leading dimension of B
                      viens, &          ! beta
                      system%hamilton%za, &  ! C
                      system%hamilton%rank &      ! LDC ... leading dimension of C
                     )
endif          
          Call timesec (ts1)
          time_hmlaan=ts1-ts0+time_hmlaan

!What if it is, say, LAPW calculation without any local orbitals?
        if (nlorb(is).ne.0) then 
! APW-LO part
          Call timesec (ts0)
          maxnlo=mt_hscf%maxnlo
          l3 = lorbl (1, is)
          lm3 = idxlm (l3, -l3)

          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                     mt_hscf%losize(is), &          ! M ... rows of op( A ) = rows of C
                     ngp, &           ! N ... cols of op( B ) = cols of C
                     maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                     viens, &          ! alpha
!                     haloij(:,:,ias), &           ! A                     
                     mt_hscf%main%loa(1,1,ias), &           ! A
                     maxnlo,&           ! LDA ... leading dimension of A
                     apwi, &           ! B
                     maxaa, &          ! LDB ... leading dimension of B
                     viens, &          ! beta
                     system%hamilton%za(ngp+idxlo (lm3, 1, ias),1), &  ! C
                     system%hamilton%rank &      ! LDC ... leading dimension of C
                     )

          do i=idxlo (lm3,1,ias)+ngp, mt_hscf%losize(is)+ngp+idxlo (lm3,1,ias)-1
            system%hamilton%za(1:ngp,i)=conjg(system%hamilton%za(i,1:ngp))
          enddo
          Call timesec (ts1)
          time_hmlalon=ts1-ts0+time_hmlalon
! LO-LO part
          Call timesec (ts0)
          l1=lorbl (1, is)
          lm1=idxlm (l1,-l1)
          j1= ngp + idxlo (lm1, 1, ias)

          l3=lorbl (nlorb(is),is)
          lm3=idxlm (l3,l3)
          j3= ngp + idxlo (lm3,nlorb(is),ias)
          system%hamilton%za(j1:j3,j1:j3)=system%hamilton%za(j1:j3,j1:j3)+ mt_hscf%main%lolo(1:1+j3-j1,1:1+j3-j1,ias) !hloloij(1:1+j3-j1,1:1+j3-j1,ias)!

          Call timesec (ts1)
          time_hmllolon=ts1-ts0+time_hmllolon
        endif

        End Do
      End Do

      deallocate(apwi,apwi2,zm)

! interstitial contributions
       Call timesec (ts0)
       Call hmlistln (system%hamilton, ngp, igpig, vgpc)
       Call timesec (ts1)
       time_hmlistln=ts1-ts0+time_hmlistln


#ifdef DEBUGHO
      Write (*,*) "apwalm", apwalm
      prefix = "H"
      Call HermitianMatrixToFiles (system%hamilton, prefix)
      prefix = "O"
      Call HermitianMatrixToFiles (system%overlap, prefix)
      Write (*,*) "wrote"
      Stop
#endif
!
      Call timesec (cpu1)
      timemat = timemat + cpu1 - cpu0

End Subroutine
