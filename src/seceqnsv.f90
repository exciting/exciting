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
      Use modmain
      Use modinput
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
      Real (8) :: ts0, ts1
! automatic arrays
      Complex (8) zftp1 (lmmaxvr), zftp2 (lmmaxvr)
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
      Allocate (bmt(lmmaxvr, nrcmtmax, 3))
      Allocate (bir(ngrtot, 3))
      Allocate (vr(nrmtmax))
      Allocate (drv(nrmtmax))
      Allocate (cf(3, nrmtmax))
      Allocate (sor(nrcmtmax))
      Allocate (rwork(3*nstsv))
      Allocate (wfmt1(lmmaxvr, nrcmtmax, nstfv))
      Allocate (wfmt2(lmmaxvr, nrcmtmax, nsc))
      Allocate (zfft1(ngrtot))
      Allocate (zfft2(ngrtot))
      Allocate (zv(ngkmax, nsc))
      lwork = 2 * nstsv
      Allocate (work(lwork))
! zero the second-variational Hamiltonian (stored in the eigenvector array)
      evecsv (:, :) = 0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            If (associated(input%groundstate%spin)) Then
! exchange-correlation magnetic field in spherical coordinates
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
! external muffin-tin magnetic field
               Do irc = 1, nrcmt (is)
                  Do i = 1, 3
                     bmt (:, irc, i) = bmt (:, irc, i) + ga4 * &
                    & (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(i)+input%groundstate%spin%bfieldc(i))
                  End Do
               End Do
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
            Do ist = 1, nstfv
               Call wavefmt (input%groundstate%lradstep, &
              & input%groundstate%lmaxvr, is, ia, ngk(1, ik), apwalm, &
              & evecfv(:, ist), lmmaxvr, wfmt1(:, :, ist))
            End Do
! begin loop over states
            Do jst = 1, nstfv
               If (associated(input%groundstate%spin)) Then
                  Do irc = 1, nrcmt (is)
! convert wavefunction to spherical coordinates
                     Call zgemv ('N', lmmaxvr, lmmaxvr, zone, zbshtvr, &
                    & lmmaxvr, wfmt1(:, irc, jst), 1, zzero, zftp1, 1)
! apply effective magnetic field and convert to spherical harmonics
                     zftp2 (:) = zftp1 (:) * bmt (:, irc, 3)
                     Call zgemv ('N', lmmaxvr, lmmaxvr, zone, zfshtvr, &
                    & lmmaxvr, zftp2, 1, zzero, wfmt2(:, irc, 1), 1)
                     wfmt2 (:, irc, 2) = - wfmt2 (:, irc, 1)
                     If (nsc .Eq. 3) Then
                        zftp2 (:) = zftp1 (:) * cmplx (bmt(:, irc, &
                       & 1),-bmt(:, irc, 2), 8)
                        Call zgemv ('N', lmmaxvr, lmmaxvr, zone, &
                       & zfshtvr, lmmaxvr, zftp2, 1, zzero, wfmt2(:, &
                       & irc, 3), 1)
                     End If
! apply spin-orbit coupling if required
                     If (isspinorb()) Then
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
                     End If
                  End Do
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
                     If (i .Le. j) Then
                        evecsv (i, j) = evecsv (i, j) + zfmtinp &
                       & (.True., input%groundstate%lmaxmat, nrcmt(is), &
                       & rcmt(:, is), lmmaxvr, wfmt1(:, :, ist), &
                       & wfmt2(:, :, k))
                     End If
                  End Do
               End Do
            End Do
! end loops over atoms and species
         End Do
      End Do
!---------------------------!
!     interstitial part     !
!---------------------------!
      If (associated(input%groundstate%spin)) Then
         If (ncmag) Then
! non-collinear
            Do ir = 1, ngrtot
               bir (ir, :) = (bxcir(ir, &
              & :)+ga4*input%groundstate%spin%bfieldc(:)) * cfunir (ir)
            End Do
         Else
! collinear
            Do ir = 1, ngrtot
               bir (ir, 1:2) = 0.d0
               bir (ir, 3) = (bxcir(ir, &
              & 1)+ga4*input%groundstate%spin%bfieldc(3)) * cfunir (ir)
            End Do
         End If
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
      End If
! add the diagonal first-variational part
      i = 0
      Do ispn = 1, nspinor
         Do ist = 1, nstfv
            i = i + 1
            evecsv (i, i) = evecsv (i, i) + evalfv (ist)
         End Do
      End Do
! diagonalise second-variational Hamiltonian
      If (ndmag .Eq. 1) Then
! collinear: block diagonalise H
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
      Deallocate (wfmt1, wfmt2, zfft1, zfft2, zv, work)
      Call timesec (ts1)
!$OMP CRITICAL
      timesv = timesv + ts1 - ts0
!$OMP END CRITICAL
      Return
20    Continue
      Write (*,*)
      Write (*, '("Error(seceqnsv): diagonalisation of the second-varia&
     &tional Hamiltonian failed")')
      Write (*, '(" for k-point ", I8)') ik
      Write (*, '(" ZHEEV returned INFO = ", I8)') info
      Write (*,*)
      Stop
End Subroutine
