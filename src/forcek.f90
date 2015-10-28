!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: forcek
! !INTERFACE:
!
!
Subroutine forcek (ik, ffacg)
! !USES:
      Use modinput
      Use modmain
      Use mod_eigensystem
      Use modfvsystem
!      use mod_constants
! !DESCRIPTION:
!   Computes the {\bf k}-dependent contribution to the incomplete basis set
!   (IBS) force. See the calling routine {\tt force} for a full description.
!
! !REVISION HISTORY:
!   Created June 2006 (JKD)
!   Updated for spin-spiral case, May 2007 (Francesco Cricchio and JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Real (8), Intent (In) :: ffacg (ngvec, nspecies)
! local variables
      Integer :: np, ispn, jspn
      Integer :: is, ia, ias, ist, jst
      Integer :: i, j, k, l, iv (3), ig
      Real (8) :: sum, t1, ta,tb
      Complex (8) zt1, zt2
      Complex (8) v (1),viens
      integer if3,l3,m3,lm3,io1,io2,maxnlo,lm1,lm2,io,j1,j2,ilo
! allocatable arrays
      Type (evsystem) :: system
      Integer, Allocatable :: ijgij(:,:)
      Real (8), Allocatable :: dpij(:,:)
      Real (8), Allocatable :: evalfv (:, :)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :, :),evecfv2(:,:,:)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: dlhij(:,:)
      Complex (8), Allocatable :: dloij(:,:)
      Complex (8), Allocatable :: ffv (:, :)
      Complex (8), Allocatable :: y (:)
      Complex (8), Allocatable :: apwi(:,:),zm(:,:),apwi2(:,:)
! external functions

      Complex (8) zdotc
      External zdotc

 !write(*,*)
      np = npmat (1, ik)
      If (isspinspiral()) np = Max (np, npmat(2, ik))
! allocate local arrays
      Allocate (evalfv(nstfv, nspnfv))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv(nmatmax, nstfv, nspnfv))
      Allocate (evecfv2(nmatmax, nstfv, nspnfv))
      Allocate (evecsv(nstsv, nstsv))
      Allocate (ffv(nstfv, nstfv))
      Allocate (y(nstfv))
! get the eigenvalues/vectors and occupancies from file
      Call getevalfv (vkl(:, ik), evalfv)
      Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
      Call getevecsv (vkl(:, ik), evecsv)
      Call getoccsv (vkl(:, ik), occsv(:, ik))

      do ispn=1, nspnfv
        do i=1, nstfv
          evecfv2(:,i,ispn)=evecfv(:,i,ispn)*evalfv(i,ispn)
        enddo
      enddo  



!      call timesec(ta)
! begin loop over first-variational spin components
      Do ispn = 1, nspnfv
! find the matching coefficients
         Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, &
        & ik), sfacgk(:, :, ispn, ik), apwalm)

       allocate(ijgij(ngk (ispn, ik),ngk (ispn, ik)))
       allocate(dpij(ngk (ispn, ik),ngk (ispn, ik)))
         Do j = 1, ngk (ispn, ik)
            Do i = 1, ngk (ispn, ik)
               iv (:) = ivg (:, igkig(i, ispn, ik)) - ivg (:, igkig(j, ispn, ik))
               iv (:) = modulo (iv(:)-intgv(:, 1), ngrid(:)) + intgv (:, 1)
               ijgij (i,j) = ivgig (iv(1), iv(2), iv(3))
               dpij (i,j) = 0.5d0 * dot_product (vgkc(:, i, ispn, ik), vgkc(:, j, ispn, ik))
            End Do
         End Do


         allocate(dlhij(nmat (ispn, ik),nmat (ispn, ik)))
         allocate(dloij(nmat (ispn, ik),nmat (ispn, ik)))
         call newsystem(system,.false.,nmat (ispn, ik))
         allocate(apwi(haaijSize,ngk(ispn, ik)))
         allocate(apwi2(ngk(ispn, ik),haaijSize))
!       call timesec(tb)
!       write(*,*) 'init',tb-ta
! loop over species and atoms
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               system%hamilton%za=dcmplx(0d0,0d0)
               system%overlap%za=dcmplx(0d0,0d0)
               ias = idxas (ia, is)


!        call timesec(ta)
          apwi=dcmplx(0d0,0d0)
          apwi2=dcmplx(0d0,0d0)
          if3=0
          Do l3 = 0, input%groundstate%lmaxmat
            Do m3 = - l3, l3
            lm3 = idxlm (l3, m3)
              Do io2 = 1, apword (l3, is)
                if3=if3+1
                apwi(if3,:)=apwalm(1:ngk(ispn, ik), io2, lm3, ias)
                apwi2(:,if3)=apwalm(1:ngk(ispn, ik), io2, lm3, ias)
              End Do
            End Do
          End Do
          allocate(zm(haaijSize,ngk(ispn, ik)))
          zm=zzero
          viens=dcmplx(1d0,0d0)
!APW-APW
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      haaijSize, &          ! M ... rows of op( A ) = rows of C
                      ngk(ispn, ik), &           ! N ... cols of op( B ) = cols of C
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
                      ngk(ispn, ik), &          ! M ... rows of op( A ) = rows of C
                      ngk(ispn, ik), &           ! N ... cols of op( B ) = cols of C
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
          deallocate(zm)
!APW-LO
        if (nlorb(is).ne.0) then
! APW-LO part
          maxnlo=size(haloij,1)
          l3 = lorbl (1, is)
          lm3 = idxlm (l3, -l3)
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                     haloijSize(is), &          ! M ... rows of op( A ) = rows of C
                     ngk (ispn, ik), &           ! N ... cols of op( B ) = cols of C
                     haaijSize, &          ! K ... cols of op( A ) = rows of op( B )
                     viens, &          ! alpha
                     haloij(:,:,ias), &           ! A
                     maxnlo,&           ! LDA ... leading dimension of A
                     apwi, &           ! B
                     haaijSize, &          ! LDB ... leading dimension of B
                     viens, &          ! beta
                     system%hamilton%za(ngk (ispn, ik)+idxlo (lm3, 1, ias),1), &  ! C
                     system%hamilton%rank &      ! LDC ... leading dimension of C
                     )
          do i=idxlo (lm3,1,ias)+ngk (ispn, ik), haloijSize(is)+ngk (ispn, ik)+idxlo (lm3,1,ias)-1
            system%hamilton%za(1:ngk (ispn, ik),i)=conjg(system%hamilton%za(i,1:ngk (ispn, ik)))
          enddo
        endif

!APW-APW
         if (input%groundstate%ValenceRelativity.ne.'iora*') then
!if (.false.) then
            call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                       'N', &           ! TRANSB = 'N'  op( B ) = B.
                       ngk (ispn, ik), &          ! M ... rows of op( A ) = rows of C
                       ngk (ispn, ik), &           ! N ... cols of op( B ) = cols of C
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
          else
            allocate(zm(ngk (ispn, ik),haaijSize))
            zm=zzero
            if3=0
            Do l3 = 0, input%groundstate%lmaxmat
              Do m3 = - l3, l3
              lm3 = idxlm (l3, m3)

                Do io2 = 1, apword (l3, is)
                  Do io1 = 1, apword (l3, is)
                    zm(:,if3+io2)=zm(:,if3+io2)+h1aa(io1,io2,l3,ias)*conjg(apwi2(:,if3+io1))
                  enddo
                  zm(:,if3+io2)=zm(:,if3+io2)+conjg(apwi2(:,if3+io2))
                End Do
                if3=if3+apword (l3, is)
              End Do
            End Do
            call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                       'C', &           ! TRANSB = 'N'  op( B ) = B.
                        ngk (ispn, ik), &          ! M ... rows of op( A ) = rows of C
                        ngk (ispn, ik), &           ! N ... cols of op( B ) = cols of C
                        haaijSize, &          ! K ... cols of op( A ) = rows of op( B )
                        viens, &          ! alpha
                        apwi, &           ! A
                        haaijSize,&           ! LDA ... leading dimension of A
                        zm, &           ! B
                        ngk (ispn, ik), &          ! LDB ... leading dimension of B
                        viens, &          ! beta
                        system%overlap%za, &  ! C
                        system%overlap%rank &      ! LDC ... leading dimension of C
                       )
             deallocate(zm)
          endif
!APW-LO
        if (nlorb(is).ne.0) then
! APW-LO part
          Do ilo = 1, nlorb (is)
            l = lorbl (ilo, is)
            lm1 = idxlm (l,-l)
            lm2 = idxlm (l, l)
            j1 = ngk (ispn, ik) + idxlo (lm1, ilo, ias)
            j2 = ngk (ispn, ik) + idxlo (lm2, ilo, ias)
            Do io = 1, apword (l, is)
              system%overlap%za(1:ngk (ispn, ik),j1:j2)=system%overlap%za(1:ngk (ispn, ik),j1:j2)+conjg(apwalm(1:ngk (ispn, ik), io, lm1:lm2, ias) * (oalo (io, ilo, ias) +h1loa(io, ilo, ias)))
            End Do
            do j=j1,j2
              system%overlap%za(j,1:ngk (ispn, ik))=conjg(system%overlap%za(1:ngk (ispn, ik),j))
            End Do
          End Do
        endif


! loop over Cartesian directions
               Do l = 1, 3

! APW-APW contribution
!        call timesec(ta)
                  Do j = 1, ngk (ispn, ik)
                     Do i = 1, ngk (ispn, ik)
                        ig = ijgij (i,j)
                        t1 = vgc (l, ig)
                        zt1 = - ffacg (ig, is) * conjg (sfacg(ig, ias))
                        dlhij (i,j) = (dpij(i,j)*zt1+system%hamilton%za(i,j)) * t1*dcmplx(0d0,1d0)
                        dloij (i,j) = -(zt1+system%overlap%za(i,j)) * t1*dcmplx(0d0,1d0)
                     End Do
                  End Do
! APW-local-orbital contribution
                  Do j = ngk (ispn, ik) + 1, nmat (ispn, ik)
                     Do i = 1, ngk (ispn, ik)
                        t1 = vgkc (l, i, ispn, ik)
                        dlhij (i,j) = system%hamilton%za(i,j) * t1*dcmplx(0d0,1d0)
                        dlhij (j,i) = conjg(dlhij (i,j))
                        dloij (i,j) = -system%overlap%za(i,j) * t1*dcmplx(0d0,1d0)
                        dloij (j,i) = conjg(dloij (i,j))
                     End Do
                  End Do
                  dlhij(ngk(ispn,ik)+1:nmat(ispn,ik),ngk(ispn,ik)+1:nmat(ispn,ik))=0d0
                  dloij(ngk(ispn,ik)+1:nmat(ispn,ik),ngk(ispn,ik)+1:nmat(ispn,ik))=0d0


!       call timesec(tb)
!       write(*,*) 'matrix2',tb-ta

!call timesec(ta)
! compute the force matrix elements in the first-variational basis
            allocate(zm(nmat(ispn, ik),nstfv))
            call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                       'N', &           ! TRANSB = 'N'  op( B ) = B.
                       nmat (ispn, ik), &          ! M ... rows of op( A ) = rows of C
                       nstfv, &          ! N ... cols of op( B ) = cols of C
                       nmat (ispn, ik), &          ! K ... cols of op( A ) = rows of op( B )
                       viens, &          ! alpha
                       dlhij, &           ! A
                       nmat (ispn, ik),&           ! LDA ... leading dimension of A
                       evecfv(:,:,ispn), &           ! B
                       nmatmax, &          ! LDB ... leading dimension of B
                       zzero, &          ! beta
                       zm, &  ! C
                       nmat(ispn, ik) &      ! LDC ... leading dimension of C
                       )
            call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                       'N', &           ! TRANSB = 'N'  op( B ) = B.
                       nmat (ispn, ik), &          ! M ... rows of op( A ) = rows of C
                       nstfv, &          ! N ... cols of op( B ) = cols of C
                       nmat (ispn, ik), &          ! K ... cols of op( A ) = rows of op( B )
                       viens, &          ! alpha
                       dloij, &           ! A
                       nmat (ispn, ik),&           ! LDA ... leading dimension of A
                       evecfv2(:,:,ispn), &           ! B
                       nmatmax, &          ! LDB ... leading dimension of B
                       viens, &          ! beta
                       zm, &  ! C
                       nmat(ispn, ik) &      ! LDC ... leading dimension of C
                       )
            call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                       'N', &           ! TRANSB = 'N'  op( B ) = B.
                       nstfv, &          ! M ... rows of op( A ) = rows of C
                       nstfv, &          ! N ... cols of op( B ) = cols of C
                       nmat (ispn, ik), &          ! K ... cols of op( A ) = rows of op( B )
                       viens, &          ! alpha 
                       evecfv(:,:,ispn), &           ! A
                       nmatmax,&           ! LDA ... leading dimension of A
                       zm, &           ! B
                       nmat(ispn, ik), &          ! LDB ... leading dimension of B
                       zzero, &          ! beta
                       ffv, &  ! C
                       nstfv &      ! LDC ... leading dimension of C
                       )

            deallocate(zm)
                   
!       call timesec(tb)
!       write(*,*) 'fv',tb-ta
!call timesec(ta)
! compute the force using the second-variational coefficients if required
                  sum = 0.d0
                  If (input%groundstate%tevecsv) Then
                     If (isspinspiral()) Then
! spin-spiral case
                        Do j = 1, nstsv
                           t1 = occsv (j, ik)
                           i = (ispn-1) * nstfv + 1
                           Call zgemv ('N', nstfv, nstfv, zone, ffv, &
                          & nstfv, evecsv(i, j), 1, zzero, y, 1)
                           zt1 = zdotc (nstfv, evecsv(i, j), 1, y, 1)
                           sum = sum + t1 * dble (zt1)
                        End Do
                     Else
! normal spin-polarised case
                        Do j = 1, nstsv
                           t1 = occsv (j, ik)
                           Do jspn = 1, nspinor
                              i = (jspn-1) * nstfv + 1
                              Call zgemv ('N', nstfv, nstfv, zone, ffv, &
                             & nstfv, evecsv(i, j), 1, zzero, y, 1)
                              zt1 = zdotc (nstfv, evecsv(i, j), 1, y, &
                             & 1)
                              sum = sum + t1 * dble (zt1)
                           End Do
                        End Do
                     End If
                  Else
! spin-unpolarised case
                     Do j = 1, nstsv
                        sum = sum + occsv (j, ik) * dble (ffv(j, j))
!                       write(*,*) ffv(j,j),occsv(j, ik)
                     End Do
!                    write(*,*) sum

                  End If
!$OMP CRITICAL
                  forceibs (l, ias) = forceibs (l, ias) + wkpt (ik) * &
                 & sum
!$OMP END CRITICAL
!       call timesec(tb)
!       write(*,*) 'sv',tb-ta

! end loop over Cartesian components
               End Do
!               write(*,*)
! end loop over atoms and species
            End Do
         End Do
! end loop over first-variational spins
      End Do
!       stop
      Deallocate (evalfv, apwalm, evecfv, evecsv, evecfv2)
      Deallocate (ffv, y)
      call deleteystem(system)
      deallocate(apwi)
      deallocate(ijgij,dpij,dlhij,dloij)

      Return
End Subroutine
!EOC
