!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genbmatk (bmt, bir, wfmt, wfir, bmat)
! calculates the magnetic field matrix elements
      Use modinput
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (In) :: bmt (lmmaxvr, nrcmtmax, natmtot, ndmag)
      Real (8), Intent (In) :: bir (ngrtot, ndmag)
      Complex (8), Intent (In) :: wfmt (lmmaxvr, nrcmtmax, natmtot, &
     & nspinor, nstsv)
      Complex (8), Intent (In) :: wfir (ngrtot, nspinor, nstsv)
      Complex (8), Intent (Out) :: bmat (nstsv, nstsv)
! local variables
      Integer :: is, ia, ias, nrc, ir, irc
      Integer :: ist, jst, idm, ispn
      Real (8) :: t1
      Complex (8) zt1
! automatic arrays
      Complex (8) zflm (lmmaxvr)
! allocatable arrays
      Real (8), Allocatable :: rvfir (:, :)
      Complex (8), Allocatable :: zfmt (:, :, :)
      Complex (8), Allocatable :: zfir (:, :)
! external functions
      Complex (8) zfmtinp, zdotc
      External zfmtinp, zdotc
      Allocate (rvfir(ngrtot, ndmag))
      Allocate (zfmt(lmmaxvr, nrcmtmax, nspinor))
      Allocate (zfir(ngrtot, nspinor))
! zero the matrix elements
      bmat (:, :) = 0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
      Do jst = 1, nstsv
         Do is = 1, nspecies
            nrc = nrcmt (is)
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
! apply magnetic field to spinor wavefunction
               Do irc = 1, nrc
                  zfmt (:, irc, 1) = bmt (:, irc, ias, ndmag) * wfmt &
                 & (:, irc, ias, 1, jst)
                  zfmt (:, irc, 2) = - bmt (:, irc, ias, ndmag) * wfmt &
                 & (:, irc, ias, 2, jst)
                  If (ncmag) Then
                     zflm (:) = cmplx (bmt(:, irc, ias, 1), bmt(:, irc, &
                    & ias, 2), 8)
                     zfmt (:, irc, 1) = zfmt (:, irc, 1) + conjg &
                    & (zflm(:)) * wfmt (:, irc, ias, 2, jst)
                     zfmt (:, irc, 2) = zfmt (:, irc, 2) + zflm (:) * &
                    & wfmt (:, irc, ias, 1, jst)
                  End If
               End Do
               Do ist = 1, jst
! compute inner product (functions are in spherical coordinates)
                  Do ispn = 1, nspinor
                     zt1 = zfmtinp (.False., input%groundstate%lmaxvr, &
                    & nrc, rcmt(:, is), lmmaxvr, wfmt(:, :, ias, ispn, &
                    & ist), zfmt(:, :, ispn))
                     bmat (ist, jst) = bmat (ist, jst) + zt1
                  End Do
               End Do
            End Do
         End Do
      End Do
!---------------------------!
!     interstitial part     !
!---------------------------!
      Do idm = 1, ndmag
         rvfir (:, idm) = bir (:, idm) * cfunir (:)
      End Do
      t1 = omega / dble (ngrtot)
      Do jst = 1, nstsv
! apply magnetic field to spinor wavefunction
         Do ir = 1, ngrtot
            zfir (ir, 1) = rvfir (ir, ndmag) * wfir (ir, 1, jst)
            zfir (ir, 2) = - rvfir (ir, ndmag) * wfir (ir, 2, jst)
         End Do
         If (ncmag) Then
            Do ir = 1, ngrtot
               zt1 = cmplx (rvfir(ir, 1), rvfir(ir, 2), 8)
               zfir (ir, 1) = zfir (ir, 1) + conjg (zt1) * wfir (ir, 2, &
              & jst)
               zfir (ir, 2) = zfir (ir, 2) + zt1 * wfir (ir, 1, jst)
            End Do
         End If
         Do ist = 1, jst
            Do ispn = 1, nspinor
               zt1 = zdotc (ngrtot, wfir(:, ispn, ist), 1, zfir(:, &
              & ispn), 1)
               bmat (ist, jst) = bmat (ist, jst) + t1 * zt1
            End Do
         End Do
      End Do
! lower triangular part
      Do ist = 1, nstsv
         Do jst = 1, ist - 1
            bmat (ist, jst) = conjg (bmat(jst, ist))
         End Do
      End Do
      Deallocate (rvfir, zfmt, zfir)
      Return
End Subroutine
