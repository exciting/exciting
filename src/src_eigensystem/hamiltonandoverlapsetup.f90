
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine hamiltonandoverlapsetup (system, ngp, apwalm, igpig, vgpc)
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
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Integer :: n
      Character (256) :: prefix
!local variables
      Integer, Save :: ikc
      Real (8), Save :: cputot
      Real (8) :: cpuaa, cpualo, cpulolo, cpui, cpu00, cpu01,ts0,ts1
      Integer :: i, is, ia
      Complex (8) v(1),viens,zzero
      Real (8) :: cpu0, cpu1
      Real (8) :: threshold
      Complex (8), allocatable :: apwi(:,:),zm(:,:)
      integer if1,if3,l3,m3,lm3,io2,ias,maxnlo,ilo,j1,j2,lm1,lm2,j,io,l
!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!
!
!
      Call timesec (cpu0)
! set the matrices to zero
!
! muffin-tin contributions

      allocate(apwi(haaijSize,ngp))
      allocate(zm(haaijSize,ngp))
      zzero=dcmplx(0d0,0d0)
      Do is = 1, nspecies
        Do ia = 1, natoms (is)
!--Hamiltonian--
! APW-APW part 
          Call timesec (ts0)
          ias = idxas (ia, is)
          apwi=dcmplx(0d0,0d0)
          if3=0
          Do l3 = 0, input%groundstate%lmaxapw
            Do m3 = - l3, l3
            lm3 = idxlm (l3, m3)
              Do io2 = 1, apword (l3, is)
                if3=if3+1
                apwi(if3,:)=apwalm(1:ngp, io2, lm3, ias)
              End Do
            End Do
          End Do
          zm=zzero
          viens=dcmplx(1d0,0)

          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      haaijSize, &          ! M ... rows of op( A ) = rows of C
                      ngp, &           ! N ... cols of op( B ) = cols of C
                      haaijSize, &          ! K ... cols of op( A ) = rows of op( B )
                      viens, &          ! alpha
                      haaij(:,:,ias), &        ! A
                      haaijSize,&           ! LDA ... leading dimension of A
                      apwi, &           ! B
                      haaijSize, &          ! LDB ... leading dimension of B
                      viens, &          ! beta
                      zm, &  ! C
                      haaijSize &      ! LDC ... leading dimension of C
                      )
          call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      ngp, &          ! M ... rows of op( A ) = rows of C
                      ngp, &           ! N ... cols of op( B ) = cols of C
                      haaijSize, &          ! K ... cols of op( A ) = rows of op( B )
                      viens, &          ! alpha
                      apwi, &           ! A
                      haaijSize,&           ! LDA ... leading dimension of A
                      zm, &           ! B
                      haaijSize, &          ! LDB ... leading dimension of B
                      viens, &          ! beta
                      system%hamilton%za, &  ! C
                      system%hamilton%rank &      ! LDC ... leading dimension of C
                     )
          Call timesec (ts1)
          time_hmlaan=ts1-ts0+time_hmlaan

! APW-LO part
          Call timesec (ts0)
          maxnlo=size(haloij,1)
          l3 = lorbl (1, is)
          lm3 = idxlm (l3, -l3)
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                     haloijSize(is), &          ! M ... rows of op( A ) = rows of C
                     ngp, &           ! N ... cols of op( B ) = cols of C
                     haaijSize, &          ! K ... cols of op( A ) = rows of op( B )
                     viens, &          ! alpha
                     haloij(:,:,ias), &           ! A
                     maxnlo,&           ! LDA ... leading dimension of A
                     apwi, &           ! B
                     haaijSize, &          ! LDB ... leading dimension of B
                     viens, &          ! beta
                     system%hamilton%za(ngp+idxlo (lm3, 1, ias),1), &  ! C
                     system%hamilton%rank &      ! LDC ... leading dimension of C
                     )
          do i=idxlo (lm3,1,ias)+ngp, haloijSize(is)+ngp+idxlo (lm3,1,ias)-1
            system%hamilton%za(1:ngp,i)=conjg(system%hamilton%za(i,1:ngp))
          enddo
          Call timesec (ts1)
          time_hmlalon=ts1-ts0+time_hmlalon
! LO-LO part
          Call timesec (ts0)
          Call hmllolon (system%hamilton, is, ia, ngp)
          Call timesec (ts1)
          time_hmllolon=ts1-ts0+time_hmllolon
!--Overlap--
! APW-APW part
          Call timesec (ts0)
          call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                     ngp, &          ! M ... rows of op( A ) = rows of C
                     ngp, &           ! N ... cols of op( B ) = cols of C
                     haaijSize, &          ! K ... cols of op( A ) = rows of op( B )
                     viens, &          ! alpha
                     apwi, &           ! A
                     haaijSize,&           ! LDA ... leading dimension of A
                     apwi, &           ! B
                     haaijSize, &          ! LDB ... leading dimension of B
                     viens, &          ! beta
                     system%overlap%za, &  ! C
                     system%overlap%rank &      ! LDC ... leading dimension of C
                     )

          Call timesec (ts1)
          time_olpaan=ts1-ts0+time_olpaan
! APW-LO part
          Call timesec (ts0)
          Do ilo = 1, nlorb (is)
            l = lorbl (ilo, is)
            lm1 = idxlm (l,-l)
            lm2 = idxlm (l, l)
            j1 = ngp + idxlo (lm1, ilo, ias)
            j2 = ngp + idxlo (lm2, ilo, ias)
            Do io = 1, apword (l, is)
              system%overlap%za(1:ngp,j1:j2)=system%overlap%za(1:ngp,j1:j2)+conjg(apwalm(:, io, lm1:lm2, ias) * oalo (io, ilo, ias))
            End Do
            do j=j1,j2
              system%overlap%za(j,1:ngp)=conjg(system%overlap%za(1:ngp,j))
            End Do
          End Do
          Call timesec (ts1)
          time_olpalon=ts1-ts0+time_olpalon
! LO-LO part
          Call timesec (ts0)
          Call olplolon (system%overlap, is, ia, ngp)
          Call timesec (ts1)
          time_olplolon=ts1-ts0+time_olplolon


         if (input%groundstate%ValenceRelativity.eq.'lkh') then
             Call timesec (ts0)
            Call hml1aan (system%h1, is, ia, ngp, apwalm)
             Call timesec (ts1)
             time_hmlaan=ts1-ts0+time_hmlaan
             Call timesec (ts0)
            Call hml1alon (system%h1, is, ia, ngp, apwalm)
             Call timesec (ts1)
             time_hmlalon=ts1-ts0+time_hmlalon
             Call timesec (ts0)
            Call hml1lolon (system%h1, is, ia, ngp)
             Call timesec (ts1)
             time_hmllolon=ts1-ts0+time_hmllolon
         endif
         End Do
      End Do
      deallocate(apwi,zm)

! interstitial contributions
       Call timesec (ts0)
      Call hmlistln (system%hamilton, ngp, igpig, vgpc)
       Call timesec (ts1)
       time_hmlistln=ts1-ts0+time_hmlistln
       Call timesec (ts0)
      Call olpistln (system%overlap, ngp, igpig)
       Call timesec (ts1)
       time_olpistln=ts1-ts0+time_olpistln
         if (input%groundstate%ValenceRelativity.eq.'lkh') then
            Call timesec (ts0)
            Call hml1istln (system%h1, ngp, igpig, vgpc)
            Call timesec (ts1)
            time_hmlistln=ts1-ts0+time_hmlistln
         endif
      threshold = 1e-16
!call HermitianMatrixTruncate(system%hamilton,threshold)
!call HermitianMatrixTruncate(system%overlap,threshold)
!
!
!
       If ( .Not. ispacked(system%hamilton)) Then
          Call hamiltonoverlapocopy_UL (system)
       End If
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
