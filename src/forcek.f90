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
      Integer, Allocatable :: ijg (:),ijgij(:,:)
      Real (8), Allocatable :: dp (:),dpij(:,:)
      Real (8), Allocatable :: evalfv (:, :)
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :, :)
      Complex (8), Allocatable :: evecsv (:, :)
!      Complex (8), Allocatable :: h (:)
!      Complex (8), Allocatable :: o (:)
      Complex (8), Allocatable :: dlh (:),dlhij(:,:)
      Complex (8), Allocatable :: dlo (:),dloij(:,:)
      Complex (8), Allocatable :: vh (:)
      Complex (8), Allocatable :: vo (:)
      Complex (8), Allocatable :: ffv (:, :)
      Complex (8), Allocatable :: y (:)
      Complex (8), Allocatable :: apwi(:,:),zm(:,:),apwi2(:,:)
! external functions

      Complex (8) zdotc
      External zdotc
      np = npmat (1, ik)
      If (isspinspiral()) np = Max (np, npmat(2, ik))
! allocate local arrays
      Allocate (ijg(np))
      Allocate (dp(np))
      Allocate (evalfv(nstfv, nspnfv))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfv(nmatmax, nstfv, nspnfv))
      Allocate (evecsv(nstsv, nstsv))
!     Allocate (h(np), o(np))
      Allocate (dlh(np), dlo(np))
      Allocate (vh(nmatmax))
      Allocate (vo(nmatmax))
      Allocate (ffv(nstfv, nstfv))
      Allocate (y(nstfv))
! get the eigenvalues/vectors and occupancies from file
      Call getevalfv (vkl(:, ik), evalfv)
      Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
      Call getevecsv (vkl(:, ik), evecsv)
      Call getoccsv (vkl(:, ik), occsv(:, ik))

      call timesec(ta)
! begin loop over first-variational spin components
      Do ispn = 1, nspnfv
! find the matching coefficients
         Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, &
        & ik), sfacgk(:, :, ispn, ik), apwalm)
         Do j = 1, ngk (ispn, ik)
            k = ((j-1)*j) / 2
            Do i = 1, j
               k = k + 1
               iv (:) = ivg (:, igkig(i, ispn, ik)) - ivg (:, igkig(j, &
              & ispn, ik))
               iv (:) = modulo (iv(:)-intgv(:, 1), ngrid(:)) + intgv &
              & (:, 1)
               ijg (k) = ivgig (iv(1), iv(2), iv(3))
               dp (k) = 0.5d0 * dot_product (vgkc(:, i, ispn, ik), &
              & vgkc(:, j, ispn, ik))
            End Do
         End Do


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
         allocate(zm(haaijSize,ngk(ispn, ik)))
       call timesec(tb)
       write(*,*) 'init',tb-ta
! loop over species and atoms
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               system%hamilton%za=dcmplx(0d0,0d0)
               system%overlap%za=dcmplx(0d0,0d0)
               ias = idxas (ia, is)


        call timesec(ta)
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
          zm=zzero
          viens=dcmplx(1d0,0)
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
          if (input%groundstate%ValenceRelativity.ne.'lkh') then
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
            deallocate(zm)
            allocate(zm(ngk (ispn, ik),haaijSize))
            zm=zzero
            if3=0
            Do l3 = 0, input%groundstate%lmaxmat
              Do m3 = - l3, l3
              lm3 = idxlm (l3, m3)

                Do io2 = 1, apword (l3, is)
                  Do io1 = 1, apword (l3, is)
!                    zm(if3+io2,:)=zm(if3+io2,:)+h1aa(io1,io2,l3,ias)*conjg(apwi2(:,if3+io1))
                    zm(:,if3+io2)=zm(:,if3+io2)+h1aa(io1,io2,l3,ias)*apwi2(:,if3+io1)
                  enddo
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
             allocate(zm(haaijSize,ngk (ispn, ik)))
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
              system%overlap%za(1:ngk (ispn, ik),j1:j2)=system%overlap%za(1:ngk (ispn, ik),j1:j2)+conjg(apwalm(:, io, lm1:lm2, ias) * oalo (io, ilo, ias))
            End Do
            do j=j1,j2
              system%overlap%za(j,1:ngk (ispn, ik))=conjg(system%overlap%za(1:ngk (ispn, ik),j))
            End Do
          End Do
        endif





! loop over Cartesian directions
               Do l = 1, 3

! APW-APW contribution
        call timesec(ta)
                  Do j = 1, ngk (ispn, ik)
                     k = ((j-1)*j) / 2
                     Do i = 1, j
                        k = k + 1
                        ig = ijg (k)
                        t1 = vgc (l, ig)
                        zt1 = - ffacg (ig, is) * conjg (sfacg(ig, ias))
                        
                        dlh (k) = (dp(k)*zt1+system%hamilton%za(i,j)) * t1
                        dlo (k) = (zt1+system%overlap%za(i,j)) * t1
                     End Do
                  End Do
! APW-local-orbital contribution
                  Do j = ngk (ispn, ik) + 1, nmat (ispn, ik)
                     k = ((j-1)*j) / 2
                     Do i = 1, ngk (ispn, ik)
                        k = k + 1
                        t1 = vgkc (l, i, ispn, ik)
                        dlh (k) = system%hamilton%za(i,j) * t1
                        dlo (k) = system%overlap%za(i,j) * t1
                     End Do
                     Do i = ngk (ispn, ik) + 1, j
                        k = k + 1
                        dlh (k) = 0.d0
                        dlo (k) = 0.d0
                     End Do
                  End Do
       call timesec(tb)
       write(*,*) 'matrix2',tb-ta

call timesec(ta)
! multiply by i
                  Do k = 1, npmat (ispn, ik)
                     dlh (k) = cmplx (-aimag(dlh(k)), dble(dlh(k)), 8)
                     dlo (k) = cmplx (-aimag(dlo(k)), dble(dlo(k)), 8)
                  End Do
! compute the force matrix elements in the first-variational basis
                  Do jst = 1, nstfv
                     Call zhpmv ('U', nmat(ispn, ik), zone, dlh, &
                    & evecfv(:, jst, ispn), 1, zzero, vh, 1)
                     Call zhpmv ('U', nmat(ispn, ik), zone, dlo, &
                    & evecfv(:, jst, ispn), 1, zzero, vo, 1)
                     t1 = evalfv (jst, ispn)
                     Do ist = 1, nstfv
                        zt1 = zdotc (nmat(ispn, ik), evecfv(:, ist, &
                       & ispn), 1, vh, 1)
                        zt2 = zdotc (nmat(ispn, ik), evecfv(:, ist, &
                       & ispn), 1, vo, 1)
                        ffv (ist, jst) = zt1 - t1 * zt2
                     End Do
                  End Do
       call timesec(tb)
       write(*,*) 'fv',tb-ta
call timesec(ta)
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
                     End Do
                  End If
!$OMP CRITICAL
                  forceibs (l, ias) = forceibs (l, ias) + wkpt (ik) * &
                 & sum
!$OMP END CRITICAL
       call timesec(tb)
       write(*,*) 'sv',tb-ta

! end loop over Cartesian components
               End Do
! end loop over atoms and species
            End Do
         End Do
! end loop over first-variational spins
      End Do
      Deallocate (ijg, dp, evalfv, apwalm, evecfv, evecsv)
      Deallocate (dlh, dlo, vh, vo, ffv, y)
      call deleteystem(system)
      deallocate(apwi)
      deallocate(ijgij,dpij,dlhij,dloij)

      Return
End Subroutine
!EOC
