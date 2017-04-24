!
!
!
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl
! F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine seceqnsv (ik, apwalm, evalfv, evecfv, evecsv)
      Use modinput, only: input, isspinorb
      Use mod_constants, only: y00, zone, zzero, zi
      Use mod_LDA_LU, only: ldapu, llu, vmatlu, lmmaxlu
      Use mod_Gvector, only: ngrtot, cfunir, igfft, ngrid
      Use mod_Gkvector, only: ngkmax, ngk, igkig
      Use mod_atoms, only: natmtot, nspecies, natoms, idxas, spr
      Use mod_muffin_tin, only: lmmaxvr, nrcmtmax, lmmaxapw, nrmtmax,&
                              & nrmt, nrcmt, idxlm, rcmt
      Use mod_potential_and_density, only: bxcmt, veffmt, bxcir, ex_coef
      Use mod_SHT, only: rbshtvr, zbshtvr, zfshtvr
      Use mod_eigensystem, only: nmatmax
      Use mod_spin, only: ncmag, nspinor, ndmag
      Use mod_eigenvalue_occupancy, only: nstfv, nstsv, evalsv
      Use mod_APW_LO, only: apwordmax
      Use mod_hybrids, only: ihyb, bxnl
      Use mod_timing, only: timesv
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Real (8), Intent (In) :: evalfv (nstfv)
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv)
      Complex (8), Intent (Out) :: evecsv (nstsv, nstsv)
! local variables
      Integer :: ispn, jspn, ia, is, ias
      Integer :: ist, jst, i, j, k, l, lm, nm
      Integer :: ir, irc, igk, ifg
      Integer :: nsc, lwork, info
! fine structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
! electron g factor
      Real (8), Parameter :: ge = 2.0023193043718d0
      Real (8), Parameter :: ga4 = ge * alpha / 4.d0
      Real (8), Parameter :: a24 = alpha ** 2 / 4.d0
      Real (8) :: rm, t1
      Real (8) :: ts0, ts1, ta,tb,tc,td
! automatic arrays
      Complex (8) zlflm (lmmaxvr, 3)
! allocatable arrays
      Real (8), Allocatable :: bmt (:, :, :)
      Real (8), Allocatable :: bir (:, :)
      Real (8), Allocatable :: vr (:)
      Real (8), Allocatable :: drv (:)
      Real (8), Allocatable :: cf (:, :)
      Real (8), Allocatable :: sor (:)
      Real (8), Allocatable :: rwork (:)
      Complex (8), Allocatable :: wfmt1 (:, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :)
      Complex (8), Allocatable :: zfft1 (:)
      Complex (8), Allocatable :: zfft2 (:)
      Complex (8), Allocatable :: zv (:, :)
      Complex (8), Allocatable :: work (:)
!     Complex (8), Allocatable :: wfmt3 (:, :)
      Complex (8) :: wfmt3 (lmmaxvr, nrcmtmax),wfmt4 (lmmaxvr, nrcmtmax)
! external functions
      Complex (8) zdotc, zfmtinp
      External zdotc, zfmtinp
! spin-unpolarised case
      If (( .Not. associated(input%groundstate%spin)) .And. (ldapu .Eq. &
     & 0)) Then
         Do i = 1, nstsv
            evalsv (i, ik) = evalfv (i)
         End Do
         evecsv (:, :) = 0.d0
         Do i = 1, nstsv
            evecsv (i, i) = 1.d0
         End Do
         Return
      End If
! number of spin combinations after application of Hamiltonian
      If (associated(input%groundstate%spin)) Then
         If ((ncmag) .Or. (isspinorb())) Then
            nsc = 3
         Else
            nsc = 2
         End If
      Else
         nsc = 1
      End If
      Call timesec (ts0)
      call timesec(ta)
      Allocate (bmt(lmmaxvr, nrcmtmax, 3))
      Allocate (bir(ngrtot, 3))
      Allocate (vr(nrmtmax))
      Allocate (drv(nrmtmax))
      Allocate (cf(3, nrmtmax))
      Allocate (sor(nrcmtmax))
      Allocate (rwork(3*nstsv))
      Allocate (wfmt1(lmmaxvr, nrcmtmax, nstfv))
      Allocate (wfmt2(lmmaxvr, nrcmtmax, nsc))
!      Allocate (zfft1(ngrtot))
!      Allocate (zfft2(ngrtot))
!      Allocate (zv(ngkmax, nsc))
      lwork = 2 * nstsv
      Allocate (work(lwork))
! zero the second-variational Hamiltonian (stored in the eigenvector array)
      evecsv (:, :) = 0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
      Do is = 1, nspecies
!        allocate(wfmt3(lmmaxvr, nrcmt(is)))
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            If (associated(input%groundstate%spin)) Then
! exchange-correlation magnetic field in spherical coordinates
               call timesec(tc)
               If (ncmag) Then
! non-collinear
                  irc = 0
                  Do ir = 1, nrmt (is), input%groundstate%lradstep
                     irc = irc + 1
                     Do i = 1, 3
                        Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, &
                       & rbshtvr, lmmaxvr, bxcmt(:, ir, ias, i), 1, &
                       & 0.d0, bmt(:, irc, i), 1)
                     End Do
                  End Do
               Else
! collinear
                  irc = 0
                  Do ir = 1, nrmt (is), input%groundstate%lradstep
                     irc = irc + 1
                     bmt (:, irc, 1:2) = 0.d0
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                    & lmmaxvr, bxcmt(:, ir, ias, 1), 1, 0.d0, bmt(:, &
                    & irc, 3), 1)
                  End Do
               End If
               call timesec(td)
!               write(*,*) td-tc
               call timesec(tc)
! external muffin-tin magnetic field
               Do irc = 1, nrcmt (is)
                  Do i = 1, 3
                     bmt (:, irc, i) = bmt (:, irc, i) + ga4 * &
                    & (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(i)+input%groundstate%spin%bfieldc(i))
                  End Do
               End Do
               call timesec(td)
!               write(*,*) td-tc
! spin-orbit radial function
               If (isspinorb()) Then
! radial derivative of the spherical part of the potential
                  vr (1:nrmt(is)) = veffmt (1, 1:nrmt(is), ias) * y00
                  Call fderiv (1, nrmt(is), spr(:, is), vr, drv, cf)
! spin-orbit coupling prefactor
                  irc = 0
                  Do ir = 1, nrmt (is), input%groundstate%lradstep
                     irc = irc + 1
                     rm = 1.d0 - 0.5d0 * (alpha**2) * vr (ir)
                     sor (irc) = a24 * drv (ir) / (spr(ir, is)*rm**2)
                  End Do
               End If
            End If
! compute the first-variational wavefunctions
            call timesec(tc)
            Do ist = 1, nstfv
               Call wavefmt (input%groundstate%lradstep, &
              & input%groundstate%lmaxvr, is, ia, ngk(1, ik), apwalm, &
              & evecfv(:, ist), lmmaxvr, wfmt1(:, :, ist))
            End Do
            call timesec(td)
!            write(*,*) td-tc
            call timesec(tc)
! begin loop over states
!xOMP PARALLEL DEFAULT(SHARED) PRIVATE(wfmt3,wfmt2,i,j,wfmt4)
!xOMP DO
            Do jst = 1, nstfv
               If (associated(input%groundstate%spin)) Then
                  wfmt4(:, :)=wfmt1(:, :, jst)
                  Call zgemm ('N', 'N', lmmaxvr, &
                                & nrcmt(is), lmmaxvr, zone, zbshtvr, &
                                & lmmaxvr, wfmt4(:,:), lmmaxvr, zzero, &
                                & wfmt3(:, :), lmmaxvr)
                  wfmt4(:, :)=wfmt3(:, :)*bmt (:, :, 3)
                  Call zgemm ('N', 'N', lmmaxvr, &
                                & nrcmt(is), lmmaxvr, zone, zfshtvr, &
                                & lmmaxvr, wfmt4(:,:), lmmaxvr, zzero, &
                                & wfmt2(:, :,1), lmmaxvr)
                  wfmt2 (:, :, 2) = - wfmt2 (:, :, 1)
                  If (nsc .Eq. 3) Then
                     wfmt4(:, :)=wfmt3(:, :) * cmplx (bmt(:, :, &
                          & 1),-bmt(:, :, 2), 8)
                     Call zgemm ('N', 'N', lmmaxvr, &
                          & nrcmt(is), lmmaxvr, zone, zfshtvr, &
                          & lmmaxvr, wfmt4(:,:), lmmaxvr, zzero, &
                          & wfmt2(:, :,3), lmmaxvr)
                  End If

                  If (isspinorb()) Then
                     Do irc = 1, nrcmt (is)
                        Call lopzflm (input%groundstate%lmaxvr, &
                             & wfmt1(:, irc, jst), lmmaxvr, zlflm)
                        t1 = sor (irc)
                        Do lm = 1, lmmaxvr
                           wfmt2 (lm, irc, 1) = wfmt2 (lm, irc, 1) + t1 &
                                & * zlflm (lm, 3)
                           wfmt2 (lm, irc, 2) = wfmt2 (lm, irc, 2) - t1 &
                                & * zlflm (lm, 3)
                           wfmt2 (lm, irc, 3) = wfmt2 (lm, irc, 3) + t1 &
                                & * (zlflm(lm, 1)-zi*zlflm(lm, 2))
                        End Do
                     End Do
                  End If
               Else
                  wfmt2 (:, :, :) = 0.d0
               End If
! apply LDA+U potential if required
               If ((ldapu .Ne. 0) .And. (llu(is) .Ge. 0)) Then
                  l = llu (is)
                  nm = 2 * l + 1
                  lm = idxlm (l,-l)
                  Do k = 1, nsc
                     If (k .Eq. 1) Then
                        ispn = 1
                        jspn = 1
                     Else If (k .Eq. 2) Then
                        ispn = 2
                        jspn = 2
                     Else
                        ispn = 1
                        jspn = 2
                     End If
                     Call zgemm ('N', 'N', nm, nrcmt(is), nm, zone, &
                    & vmatlu(lm, lm, ispn, jspn, ias), lmmaxlu, &
                    & wfmt1(lm, 1, jst), lmmaxvr, zone, wfmt2(lm, 1, &
                    & k), lmmaxvr)
                  End Do
               End If
! second-variational Hamiltonian matrix
               Do ist = 1, nstfv
                  Do k = 1, nsc
                     If (k .Eq. 1) Then
                        i = ist
                        j = jst
                     Else If (k .Eq. 2) Then
                        i = ist + nstfv
                        j = jst + nstfv
                     Else
                        i = ist
                        j = jst + nstfv
                     End If
!                     If (i .Le. j) Then
                        evecsv (i, j) = evecsv (i, j) + zfmtinp &
                       & (.True., input%groundstate%lmaxmat, nrcmt(is), &
                       & rcmt(:, is), lmmaxvr, wfmt1(:, :, ist), &
                       & wfmt2(:, :, k))
!                     End If
                  End Do
               End Do
            End Do
!xOMP END DO
!xOMP END PARALLEL
            call timesec(td)
!            write(*,*) td-tc
! end loops over atoms and species
         End Do
!         deallocate(wfmt3)
      End Do
      call timesec(tb)
!      write(*,*) 'sv / MT part',tb-ta
!---------------------------!
!     interstitial part     !
!---------------------------!
      call timesec(ta)
      If (associated(input%groundstate%spin)) Then
         If (ncmag) Then
! non-collinear
            Do ir = 1, ngrtot
               bir (ir, :) = &
               &  (bxcir(ir,:)+ga4*input%groundstate%spin%bfieldc(:))*cfunir(ir)
            End Do
         Else
! collinear
            Do ir = 1, ngrtot
               bir (ir, 1:2) = 0.d0
               bir (ir, 3) = (bxcir(ir, &
              & 1)+ga4*input%groundstate%spin%bfieldc(3)) * cfunir (ir)
            End Do
         End If
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) SHARED(nstfv,ngk,igfft,igkig,ngrid,ik,evecfv,evecsv,nsc,bir,ngkmax,ngrtot) PRIVATE(zfft1,zfft2,jst,ist,igk,ifg,zv,k,i,j)
#endif
      Allocate (zfft1(ngrtot))
      Allocate (zfft2(ngrtot))
      Allocate (zv(ngkmax, nsc))
#ifdef USEOMP
!$OMP DO
#endif 

         Do jst = 1, nstfv
            zfft1 (:) = 0.d0
            Do igk = 1, ngk (1, ik)
               ifg = igfft (igkig(igk, 1, ik))
               zfft1 (ifg) = evecfv (igk, jst)
            End Do
! Fourier transform wavefunction to real-space
            Call zfftifc (3, ngrid, 1, zfft1)
! multiply with magnetic field and transform to G-space
            zfft2 (:) = zfft1 (:) * bir (:, 3)
            Call zfftifc (3, ngrid,-1, zfft2)
            Do igk = 1, ngk (1, ik)
               ifg = igfft (igkig(igk, 1, ik))
               zv (igk, 1) = zfft2 (ifg)
               zv (igk, 2) = - zfft2 (ifg)
            End Do
            If (nsc .Eq. 3) Then
               zfft2 (:) = zfft1 (:) * cmplx (bir(:, 1),-bir(:, 2), 8)
               Call zfftifc (3, ngrid,-1, zfft2)
               Do igk = 1, ngk (1, ik)
                  ifg = igfft (igkig(igk, 1, ik))
                  zv (igk, 3) = zfft2 (ifg)
               End Do
            End If
! add to Hamiltonian matrix
            Do ist = 1, nstfv
               Do k = 1, nsc
                  If (k .Eq. 1) Then
                     i = ist
                     j = jst
                  Else If (k .Eq. 2) Then
                     i = ist + nstfv
                     j = jst + nstfv
                  Else
                     i = ist
                     j = jst + nstfv
                  End If
                  If (i .Le. j) Then
                     evecsv (i, j) = evecsv (i, j) + zdotc (ngk(1, ik), &
                    & evecfv(:, ist), 1, zv(:, k), 1)
                  End If
               End Do
            End Do
         End Do
#ifdef USEOMP
!$OMP END DO
#endif
      deallocate (zfft1,zfft2,zv)
#ifdef USEOMP
!$OMP END PARALLEL
#endif 
      End If

!-----------------------------------------
! add the diagonal first-variational part
!-----------------------------------------
      i = 0
      Do ispn = 1, nspinor
         Do ist = 1, nstfv
            i = i + 1
            evecsv (i, i) = evecsv (i, i) + evalfv (ist)
         End Do
      End Do
      
!---------------------------------------------------------      
! add the second-variational \Sigma_x matrix elements
!---------------------------------------------------------
      if (associated(input%groundstate%Hybrid)) then
         if (input%groundstate%Hybrid%exchangetypenumber == 1) then
            ! Update Hamiltonian
            if (ihyb>0) evecsv(:,:) = &
            &  evecsv(:,:) + ex_coef*bxnl(:,:,ik)
         end if
      end if      
      
! diagonalise second-variational Hamiltonian
      call timesec(ta)
      If (ndmag .Eq. 1) Then
! collinear: block diagonalise H
!         do i=1,nstfv
!           write(*,*) dble(evecsv(i,i)),dimag(evecsv(i,i))
!         enddo
!       do i=1,nstfv
!        do j=1,nstfv
!          write(*,*) dble(evecsv(j,i)),dimag(evecsv(j,i))
!        enddo
!       enddo
!stop
!         write(*,*) dble(sum(evecsv(1:nstfv,1:nstfv))),dimag(sum(evecsv(1:nstfv,1:nstfv)))
!         write(*,*) dble(sum(evecsv(nstfv+1:2*nstfv,nstfv+1:2*nstfv))),dimag(sum(evecsv(nstfv+1:2*nstfv,nstfv+1:2*nstfv)))
!         write(*,*) 'sv'
!         stop
         Call zheev ('V', 'U', nstfv, evecsv, nstsv, evalsv(:, ik), &
        & work, lwork, rwork, info)
         If (info .Ne. 0) Go To 20
         i = nstfv + 1
         Call zheev ('V', 'U', nstfv, evecsv(i, i), nstsv, evalsv(i, &
        & ik), work, lwork, rwork, info)
         If (info .Ne. 0) Go To 20
         Do i = 1, nstfv
            Do j = 1, nstfv
               evecsv (i, j+nstfv) = 0.d0
               evecsv (i+nstfv, j) = 0.d0
            End Do
         End Do
      Else
! non-collinear or spin-unpolarised: full diagonalisation
         Call zheev ('V', 'U', nstsv, evecsv, nstsv, evalsv(:, ik), &
        & work, lwork, rwork, info)
         If (info .Ne. 0) Go To 20
      End If
      Deallocate (bmt, bir, vr, drv, cf, sor, rwork)
      Deallocate (wfmt1, wfmt2, work)
      Call timesec (ts1)
      call timesec(tb)
!      write(*,*) 'sv / diagonalization', tb-ta

      timesv = timesv + ts1 - ts0       
!$OMP CRITICAL
!!      timesv = timesv + ts1 - ts0
!$OMP END CRITICAL
      Return
20    Continue
      Write (*,*)
      Write (*, '("Error(seceqnsv):& 
     & diagonalisation of the second-variational Hamiltonian failed")')
      Write (*, '(" for k-point ", I8)') ik
      Write (*, '(" ZHEEV returned INFO = ", I8)') info
      Write (*,*)
      Stop
End Subroutine
