!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genvmatk (vmt, vir, wfmt, wfir, vmat)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Real (8), Intent (In) :: vmt (lmmaxvr, nrcmtmax, natmtot)
      Real (8), Intent (In) :: vir (ngrtot)
      Complex (8), Intent (In) :: wfmt (lmmaxvr, nrcmtmax, natmtot, &
     & nspinor, nstsv)
      Complex (8), Intent (In) :: wfir (ngrtot, nspinor, nstsv)
      Complex (8), Intent (Out) :: vmat (nstsv, nstsv)
! local variables
      Integer :: is, ia, ias, nrc, irc
      Integer :: ispn, ist, jst
      Real (8) :: t1
      Complex (8) zt1
! allocatable arrays
      Real (8), Allocatable :: rfir (:)
      Complex (8), Allocatable :: zfmt (:, :)
      Complex (8), Allocatable :: zfir (:)
! external functions
      Complex (8) zfmtinp, zdotc
      External zfmtinp, zdotc
! allocate local arrays
      Allocate (rfir(ngrtot))
      Allocate (zfmt(lmmaxvr, nrcmtmax))
      Allocate (zfir(ngrtot))
! zero the matrix elements
      vmat (:, :) = 0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
      Do jst = 1, nstsv
         Do is = 1, nspecies
            nrc = nrcmt (is)
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ispn = 1, nspinor
! apply potential to wavefunction
                  Do irc = 1, nrc
                     zfmt (:, irc) = vmt (:, irc, ias) * wfmt (:, irc, &
                    & ias, ispn, jst)
                  End Do
                  Do ist = 1, jst
! compute inner product (functions are in spherical coordinates)
                     zt1 = zfmtinp (.False., input%groundstate%lmaxvr, &
                    & nrc, rcmt(:, is), lmmaxvr, wfmt(:, :, ias, ispn, &
                    & ist), zfmt)
                     vmat (ist, jst) = vmat (ist, jst) + zt1
                  End Do
               End Do
            End Do
         End Do
      End Do
!---------------------------!
!     interstitial part     !
!---------------------------!
      rfir (:) = vir (:) * cfunir (:)
      t1 = omega / dble (ngrtot)
      Do jst = 1, nstsv
         Do ispn = 1, nspinor
! apply potential to wavefunction
            zfir (:) = rfir (:) * wfir (:, ispn, jst)
            Do ist = 1, jst
               zt1 = zdotc (ngrtot, wfir(:, ispn, ist), 1, zfir, 1)
               vmat (ist, jst) = vmat (ist, jst) + t1 * zt1
            End Do
         End Do
      End Do
! lower triangular part
      Do ist = 1, nstsv
         Do jst = 1, ist - 1
            vmat (ist, jst) = conjg (vmat(jst, ist))
         End Do
      End Do
      Deallocate (rfir, zfmt, zfir)
      Return
End Subroutine
