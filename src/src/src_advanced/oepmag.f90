!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine oepmag (tsh, wfmt1, wfmt2, wfir1, wfir2, zmagmt, zmagir)
      Use modmain
      Implicit None
! arguments
      Logical, Intent (In) :: tsh
      Complex (8), Intent (In) :: wfmt1 (lmmaxvr, nrcmtmax, natmtot, &
     & nspinor)
      Complex (8), Intent (In) :: wfmt2 (lmmaxvr, nrcmtmax, natmtot, &
     & nspinor)
      Complex (8), Intent (In) :: wfir1 (ngrtot, nspinor)
      Complex (8), Intent (In) :: wfir2 (ngrtot, nspinor)
      Complex (8), Intent (Out) :: zmagmt (lmmaxvr, nrcmtmax, natmtot, &
     & ndmag)
      Complex (8), Intent (Out) :: zmagir (ngrtot, ndmag)
! local variables
      Integer :: is, ia, ias, nrc, ir, idm
      Complex (8) zt1, zt2
! allocatable arrays
      Complex (8), Allocatable :: zvfmt (:, :, :)
      Allocate (zvfmt(lmmaxvr, nrcmtmax, ndmag))
! muffin-tin part
      Do is = 1, nspecies
         nrc = nrcmt (is)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Call oepmagmt (tsh, is, wfmt1(:, :, ias, 1), wfmt1(:, :, &
           & ias, 2), wfmt2(:, :, ias, 1), wfmt2(:, :, ias, 2), zvfmt)
            Do idm = 1, ndmag
               zmagmt (:, 1:nrc, ias, idm) = zvfmt (:, 1:nrc, idm)
            End Do
         End Do
      End Do
! interstitial part
      Do ir = 1, ngrtot
! calculate the z-component of mangetisation: up-up - dn-dn
         zmagir (ir, ndmag) = conjg (wfir1(ir, 1)) * wfir2 (ir, 1) - &
        & conjg (wfir1(ir, 2)) * wfir2 (ir, 2)
         If (ncmag) Then
! up-dn spin density
            zt1 = conjg (wfir1(ir, 1)) * wfir2 (ir, 2)
! dn-up spin density
            zt2 = conjg (wfir1(ir, 2)) * wfir2 (ir, 1)
! calculate the x-component: up-dn + dn-up
            zmagir (ir, 1) = zt1 + zt2
! calculate the y-component: i*(dn-up - up-dn)
            zmagir (ir, 2) = zi * (zt2-zt1)
         End If
      End Do
      Deallocate (zvfmt)
      Return
End Subroutine
