!
!
!
! Copyright (C) 2006 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine seceqnss (ik, apwalm, evalfv, evecfv, evecsv)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot, nspnfv)
      Real (8), Intent (In) :: evalfv (nstfv, nspnfv)
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv, nspnfv)
      Complex (8), Intent (Out) :: evecsv (nstsv, nstsv)
! local variables
      Integer :: ispn, jspn, is, ia, ias
      Integer :: ist, jst, i, j, k, l, lm, nm
      Integer :: ir, irc, igk, ifg
      Integer :: lwork, info
! fine structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
! electron g factor
      Real (8), Parameter :: ge = 2.0023193043718d0
      Real (8), Parameter :: ga4 = ge * alpha / 4.d0
      Real (8) :: ts0, ts1
! automatic arrays
      Complex (8) zftp1 (lmmaxvr, nspnfv), zftp2 (lmmaxvr)
! allocatable arrays
      Real (8), Allocatable :: bmt (:, :, :)
      Real (8), Allocatable :: bir (:, :)
      Real (8), Allocatable :: rwork (:)
      Complex (8), Allocatable :: wfmt1 (:, :, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :)
      Complex (8), Allocatable :: zfft1 (:, :)
      Complex (8), Allocatable :: zfft2 (:)
      Complex (8), Allocatable :: zv (:, :)
      Complex (8), Allocatable :: work (:)
! external functions
      Complex (8) zdotc
      Complex (8) zfmtinp
      External zdotc, zfmtinp
      If ( .Not. associated(input%groundstate%spin)) Then
         Write (*,*)
         Write (*, '("Error(seceqnss): spin-unpolarised calculation")')
         Write (*,*)
         Stop
      End If
      Call timesec (ts0)
      Allocate (bmt(lmmaxvr, nrcmtmax, 3))
      Allocate (bir(ngrtot, 3))
      Allocate (rwork(3*nstsv))
      Allocate (wfmt1(lmmaxvr, nrcmtmax, nstfv, nspnfv))
      Allocate (wfmt2(lmmaxvr, nrcmtmax, 3))
      Allocate (zfft1(ngrtot, nspnfv))
      Allocate (zfft2(ngrtot))
      Allocate (zv(ngkmax, 3))
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
! exchange-correlation magnetic field in spherical coordinates
            irc = 0
            Do ir = 1, nrmt (is), input%groundstate%lradstep
               irc = irc + 1
               Do i = 1, 3
                  Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                 & lmmaxvr, bxcmt(:, ir, ias, i), 1, 0.d0, bmt(:, irc, &
                 & i), 1)
               End Do
            End Do
! external muffin-tin magnetic field
            Do irc = 1, nrcmt (is)
               Do i = 1, 3
                  bmt (:, irc, i) = bmt (:, irc, i) + ga4 * &
                 & (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(i)+input%groundstate%spin%bfieldc(i))
               End Do
            End Do
! compute the first-variational wavefunctions
            Do ispn = 1, nspnfv
               Do ist = 1, nstfv
                  Call wavefmt (input%groundstate%lradstep, &
                 & input%groundstate%lmaxvr, is, ia, ngk(ispn, ik), &
                 & apwalm(:, :, :, :, ispn), evecfv(:, ist, ispn), &
                 & lmmaxvr, wfmt1(:, :, ist, ispn))
               End Do
            End Do
            Do jst = 1, nstfv
               Do irc = 1, nrcmt (is)
! convert wavefunctions to spherical coordinates
                  Do ispn = 1, nspnfv
                     Call zgemv ('N', lmmaxvr, lmmaxvr, zone, zbshtvr, &
                    & lmmaxvr, wfmt1(:, irc, jst, ispn), 1, zzero, &
                    & zftp1(:, ispn), 1)
                  End Do
! apply effective magnetic field and convert to spherical harmonics
                  zftp2 (:) = zftp1 (:, 1) * bmt (:, irc, 3)
                  Call zgemv ('N', lmmaxvr, lmmaxvr, zone, zfshtvr, &
                 & lmmaxvr, zftp2, 1, zzero, wfmt2(:, irc, 1), 1)
                  zftp2 (:) = - zftp1 (:, 2) * bmt (:, irc, 3)
                  Call zgemv ('N', lmmaxvr, lmmaxvr, zone, zfshtvr, &
                 & lmmaxvr, zftp2, 1, zzero, wfmt2(:, irc, 2), 1)
                  zftp2 (:) = zftp1 (:, 2) * cmplx (bmt(:, irc, &
                 & 1),-bmt(:, irc, 2), 8)
                  Call zgemv ('N', lmmaxvr, lmmaxvr, zone, zfshtvr, &
                 & lmmaxvr, zftp2, 1, zzero, wfmt2(:, irc, 3), 1)
               End Do
! apply LDA+U potential if required
               If ((ldapu .Ne. 0) .And. (llu(is) .Ge. 0)) Then
                  l = llu (is)
                  nm = 2 * l + 1
                  lm = idxlm (l,-l)
                  Do k = 1, 3
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
                    & wfmt1(lm, 1, jst, jspn), lmmaxvr, zone, wfmt2(lm, &
                    & 1, k), lmmaxvr)
                  End Do
               End If
! second-variational Hamiltonian matrix
               Do ist = 1, nstfv
                  Do k = 1, 3
                     If (k .Eq. 1) Then
                        ispn = 1
                        i = ist
                        j = jst
                     Else If (k .Eq. 2) Then
                        ispn = 2
                        i = ist + nstfv
                        j = jst + nstfv
                     Else
                        ispn = 1
                        i = ist
                        j = jst + nstfv
                     End If
                     If (i .Le. j) Then
                        evecsv (i, j) = evecsv (i, j) + zfmtinp &
                       & (.True., input%groundstate%lmaxmat, nrcmt(is), &
                       & rcmt(:, is), lmmaxvr, wfmt1(:, :, ist, ispn), &
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
      Do ir = 1, ngrtot
         bir (ir, :) = (bxcir(ir, &
        & :)+ga4*input%groundstate%spin%bfieldc(:)) * cfunir (ir)
      End Do
      Do jst = 1, nstfv
         Do ispn = 1, nspnfv
            zfft1 (:, ispn) = 0.d0
            Do igk = 1, ngk (ispn, ik)
               ifg = igfft (igkig(igk, ispn, ik))
               zfft1 (ifg, ispn) = evecfv (igk, jst, ispn)
            End Do
! Fourier transform wavefunction to real-space
            Call zfftifc (3, ngrid, 1, zfft1(:, ispn))
         End Do
! multiply with magnetic field and transform to G-space
         Do k = 1, 3
            If (k .Eq. 1) Then
               ispn = 1
               zfft2 (:) = zfft1 (:, 1) * bir (:, 3)
            Else If (k .Eq. 2) Then
               ispn = 2
               zfft2 (:) = - zfft1 (:, 2) * bir (:, 3)
            Else
               ispn = 1
               zfft2 (:) = zfft1 (:, 2) * cmplx (bir(:, 1),-bir(:, 2), &
              & 8)
            End If
            Call zfftifc (3, ngrid,-1, zfft2)
            Do igk = 1, ngk (ispn, ik)
               ifg = igfft (igkig(igk, ispn, ik))
               zv (igk, k) = zfft2 (ifg)
            End Do
         End Do
         Do ist = 1, nstfv
            Do k = 1, 3
               If (k .Eq. 1) Then
                  ispn = 1
                  i = ist
                  j = jst
               Else If (k .Eq. 2) Then
                  ispn = 2
                  i = ist + nstfv
                  j = jst + nstfv
               Else
                  ispn = 1
                  i = ist
                  j = jst + nstfv
               End If
               If (i .Le. j) Then
                  evecsv (i, j) = evecsv (i, j) + zdotc (ngk(ispn, ik), &
                 & evecfv(:, ist, ispn), 1, zv(:, k), 1)
               End If
            End Do
         End Do
      End Do
! add the diagonal first-variational part
      i = 0
      Do ispn = 1, nspinor
         Do ist = 1, nstfv
            i = i + 1
            evecsv (i, i) = evecsv (i, i) + evalfv (ist, ispn)
         End Do
      End Do
! diagonalise the Hamiltonian
      Call zheev ('V', 'U', nstsv, evecsv, nstsv, evalsv(:, ik), work, &
     & lwork, rwork, info)
      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(seceqnss): diagonalisation of the second-va&
        &riational Hamiltonian failed")')
         Write (*, '(" for k-point ", I8)') ik
         Write (*, '(" ZHEEV returned INFO = ", I8)') info
         Write (*,*)
         Stop
      End If
      Deallocate (bmt, bir, rwork)
      Deallocate (wfmt1, wfmt2, zfft1, zfft2, work)
      Call timesec (ts1)
      timesv = timesv + ts1 - ts0
!$OMP CRITICAL
!!      timesv = timesv + ts1 - ts0
!$OMP END CRITICAL
      Return
End Subroutine
