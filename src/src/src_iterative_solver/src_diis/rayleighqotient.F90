
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine rayleighqotient (n, m, evecfv, h, s, evalfv)
!
      Implicit None
      Integer, Intent (In) :: n, m
      Complex (8), Intent (In) :: h (n, m), s (n, m), evecfv (n, m)
      Real (8), Intent (Out) :: evalfv (m)
      Complex (8) zdotc
      External zdotc
      Complex (8) :: vwork (n)
      Complex (8) :: vhv, vsv
      Integer :: i
      Do i = 1, m
         vhv = zdotc (n, evecfv(1, i), 1, h(1, i), 1)
         vsv = zdotc (n, evecfv(1, i), 1, s(1, i), 1)
         evalfv (i) = vhv / vsv
!
      End Do
End Subroutine rayleighqotient
