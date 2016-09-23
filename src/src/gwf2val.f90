!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine gwf2val (ik, evecfv, evecsv, gw2fmt, gw2fir)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv)
      Complex (8), Intent (In) :: evecsv (nstsv, nstsv)
      Real (8), Intent (Inout) :: gw2fmt (lmmaxvr, nrmtmax, natmtot)
      Real (8), Intent (Inout) :: gw2fir (ngrtot)
! local variables
      Integer :: ispn, ist, is, ia, ias
      Integer :: ir, itp, igk, ifg, i, j, n
      Real (8) :: wo, t1
      Complex (8) zt1
! automatic arrays
      Logical :: done (nstfv)
      Complex (8) zftp (lmmaxvr)
! allocatable arrays
      Complex (8), Allocatable :: apwalm (:, :, :, :)
      Complex (8), Allocatable :: wfmt1 (:, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :)
      Complex (8), Allocatable :: gzfmt (:, :, :)
      Complex (8), Allocatable :: zfft1 (:, :)
      Complex (8), Allocatable :: zfft2 (:)
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (wfmt1(lmmaxvr, nrmtmax, nstfv))
      Allocate (wfmt2(lmmaxvr, nrmtmax, nspinor))
      Allocate (gzfmt(lmmaxvr, nrmtmax, 3))
      Allocate (zfft1(ngrtot, nspinor))
      Allocate (zfft2(ngrtot))
! find the matching coefficients
      Call match (ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, 1, ik), &
     & sfacgk(:, :, 1, ik), apwalm)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
      Do is = 1, nspecies
         n = lmmaxvr * nrmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            done (:) = .False.
            Do j = 1, nstsv
               wo = wkpt (ik) * occsv (j, ik)
               If (Abs(wo) .Gt. input%groundstate%epsocc) Then
                  If (input%groundstate%tevecsv) Then
! generate spinor wavefunction from second-variational eigenvectors
                     wfmt2 (:, :, :) = 0.d0
                     i = 0
                     Do ispn = 1, nspinor
                        Do ist = 1, nstfv
                           i = i + 1
                           zt1 = evecsv (i, j)
                           If (Abs(dble(zt1))+Abs(aimag(zt1)) .Gt. &
                          & input%groundstate%epsocc) Then
                              If ( .Not. done(ist)) Then
                                 Call wavefmt (1, &
                                & input%groundstate%lmaxvr, is, ia, &
                                & ngk(1, ik), apwalm, evecfv(:, ist), &
                                & lmmaxvr, wfmt1(:, :, ist))
                                 done (ist) = .True.
                              End If
! add to spinor wavefunction
                              Call zaxpy (n, zt1, wfmt1(:, :, ist), 1, &
                             & wfmt2(:, :, ispn), 1)
                           End If
                        End Do
                     End Do
                  Else
! spin-unpolarised wavefunction
                     Call wavefmt (1, input%groundstate%lmaxvr, is, ia, &
                    & ngk(1, ik), apwalm, evecfv(:, j), lmmaxvr, wfmt2)
                  End If
! compute the gradient of the wavefunction
                  Do ispn = 1, nspinor
                     Call gradzfmt (input%groundstate%lmaxvr, nrmt(is), &
                    & spr(:, is), lmmaxvr, nrmtmax, wfmt2(:, :, ispn), &
                    & gzfmt)
! convert gradient from spherical harmonics to spherical coordinates
                     Do i = 1, 3
                        Do ir = 1, nrmt (is)
                           Call zgemv ('N', lmmaxvr, lmmaxvr, zone, &
                          & zbshtvr, lmmaxvr, gzfmt(:, ir, i), 1, &
                          & zzero, zftp, 1)
                           Do itp = 1, lmmaxvr
                              t1 = wo * &
                             & (dble(zftp(itp))**2+aimag(zftp(itp))**2)
                              gw2fmt (itp, ir, ias) = gw2fmt (itp, ir, &
                             & ias) + t1
                           End Do
                        End Do
                     End Do
                  End Do
               End If
            End Do
         End Do
      End Do
!---------------------------!
!     interstitial part     !
!---------------------------!
      Do j = 1, nstsv
         wo = wkpt (ik) * occsv (j, ik)
         If (Abs(wo) .Gt. input%groundstate%epsocc) Then
            t1 = wo / omega
            zfft1 (:, :) = 0.d0
            If (input%groundstate%tevecsv) Then
! generate spinor wavefunction from second-variational eigenvectors
               i = 0
               Do ispn = 1, nspinor
                  Do ist = 1, nstfv
                     i = i + 1
                     zt1 = evecsv (i, j)
                     If (Abs(dble(zt1))+Abs(aimag(zt1)) .Gt. &
                    & input%groundstate%epsocc) Then
                        Do igk = 1, ngk (1, ik)
                           ifg = igfft (igkig(igk, 1, ik))
                           zfft1 (ifg, ispn) = zfft1 (ifg, ispn) + zt1 &
                          & * evecfv (igk, ist)
                        End Do
                     End If
                  End Do
               End Do
            Else
! spin-unpolarised wavefunction
               Do igk = 1, ngk (1, ik)
                  ifg = igfft (igkig(igk, 1, ik))
                  zfft1 (ifg, 1) = evecfv (igk, j)
               End Do
            End If
! compute gradient of wavefunction
            Do ispn = 1, nspinor
               Do i = 1, 3
                  zfft2 (:) = 0.d0
                  Do igk = 1, ngk (1, ik)
                     ifg = igfft (igkig(igk, 1, ik))
                     zfft2 (ifg) = zi * vgkc (i, igk, 1, ik) * zfft1 &
                    & (ifg, ispn)
                  End Do
! Fourier transform gradient to real-space
                  Call zfftifc (3, ngrid, 1, zfft2)
                  Do ir = 1, ngrtot
                     gw2fir (ir) = gw2fir (ir) + t1 * &
                    & (dble(zfft2(ir))**2+aimag(zfft2(ir))**2)
                  End Do
               End Do
            End Do
         End If
      End Do
      Deallocate (apwalm, wfmt1, wfmt2, gzfmt, zfft1, zfft2)
      Return
End Subroutine
