
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine setupprojectedhamilton (n, m, h, o, nmatmax, evecfv, evalfv, &
& da, hprojected, oprojected)
      Implicit None
      Integer, Intent (In) :: n, m, nmatmax
      Complex (8), Intent (In) :: evecfv (nmatmax, m), h (n*(n+1)/2), o &
     & (n*(n+1)/2), da (n, m)
      Real (8), Intent (In) :: evalfv (m)
      Complex (8), Intent (Out) :: hprojected (2*m*(2*m+1)/2), &
     & oprojected (2*m*(2*m+1)/2)
      Complex (8) :: basis (n, 2*m), vec (n)
      Complex (8) :: zdotc
      External zdotc
      Integer :: i, j, pi
      Do i = 1, 2 * m
         If (i .Le. m) Then
            Call zcopy (n, evecfv(1, i), 1, basis(1, i), 1)
         Else
            Call zcopy (n, da(1, i-m), 1, basis(1, i), 1)
         End If
      End Do
      pi = 1
      Do i = 1, 2 * m
         Do j = 1, i
            vec = 0.0
            Call zhpmv ('U', n, (1.d0, 0.d0), h, basis(1, j), 1, (0.d0, &
           & 0.d0), vec, 1)
            hprojected (pi) = zdotc (n, basis(1, i), 1, vec, 1)
            oprojected (pi) = zdotc (n, basis(1, i), 1, basis(1, j), 1)
            pi = pi + 1
         End Do
      End Do
!
End Subroutine setupprojectedhamilton
