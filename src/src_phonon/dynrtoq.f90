!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine dynrtoq (vpl, dynr, dynp)
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (In) :: vpl (3)
      Complex (8), Intent (In) :: dynr (3*natmtot, 3*natmtot, &
     & ngridq(1)*ngridq(2)*ngridq(3))
      Complex (8), Intent (Out) :: dynp (3*natmtot, 3*natmtot)
! local variables
      Integer :: i1, i2, i3, ir, i, j
      Real (8) :: t1
      Complex (8) zt1
      dynp (:, :) = 0.d0
! loop over R-vectors
      ir = 0
      Do i3 = ngridq (3) / 2 - ngridq (3) + 1, ngridq (3) / 2
         Do i2 = ngridq (2) / 2 - ngridq (2) + 1, ngridq (2) / 2
            Do i1 = ngridq (1) / 2 - ngridq (1) + 1, ngridq (1) / 2
               ir = ir + 1
               t1 = - twopi * &
              & (vpl(1)*dble(i1)+vpl(2)*dble(i2)+vpl(3)*dble(i3))
               zt1 = cmplx (Cos(t1), Sin(t1), 8)
               Do i = 1, 3 * natmtot
                  Do j = 1, 3 * natmtot
                     dynp (i, j) = dynp (i, j) + zt1 * dynr (i, j, ir)
                  End Do
               End Do
            End Do
         End Do
      End Do
! symmetrise the dynamical matrix
      Call dynsym (vpl, dynp)
      Return
End Subroutine
