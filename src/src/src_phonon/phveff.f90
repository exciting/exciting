!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine phveff
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: is, ia, ias, jas, i
      Integer :: i1, i2, i3, ir
      Real (8) :: v1 (3), v2 (3)
! external functions
      Real (8) :: rfirvec
      External rfirvec
! muffin-tin density
      ias = 0
      jas = 0
      Do is = 1, nspecies
         Do ia = 1, natoms0 (is)
            ias = ias + 1
            Do i = 1, nphcell
               jas = jas + 1
               veffmt (:, :, jas) = veffmt0 (:, :, ias)
            End Do
         End Do
      End Do
! interstitial density
      ir = 0
      Do i3 = 0, ngrid (3) - 1
         v1 (3) = dble (i3) / dble (ngrid(3))
         Do i2 = 0, ngrid (2) - 1
            v1 (2) = dble (i2) / dble (ngrid(2))
            Do i1 = 0, ngrid (1) - 1
               v1 (1) = dble (i1) / dble (ngrid(1))
               ir = ir + 1
               Call r3mv (input%structure%crystal%basevect, v1, v2)
               veffir (ir) = rfirvec (ngrid0, ainv0, v2, veffir0)
            End Do
         End Do
      End Do
      Call genveffig
      Return
End Subroutine
