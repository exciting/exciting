
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine updateevecfv (n, m, da, nmatmax, evecfv, evalfv, evecp, &
& evalp)
      Implicit None
      Integer, Intent (In) :: n, m, nmatmax
      Complex (8), Intent (In) :: da (n, m), evecp (2*m, m)
      Complex (8), Intent (Inout) :: evecfv (nmatmax, m)
      Real (8), Intent (In) :: evalp (2*m)
      Real (8), Intent (Out) :: evalfv (m)
      Integer :: i, j
!
!local vars
      Complex (8) :: basis (n, 2*m)
      Do i = 1, 2 * m
         If (i .Le. m) Then
            Call zcopy (n, evecfv(1, i), 1, basis(1, i), 1)
         Else
            Call zcopy (n, da(1, i-m), 1, basis(1, i), 1)
         End If
      End Do
      evecfv (:, :) = 0.d0
      Do i = 1, m
         Do j = 1, 2 * m
            Call zaxpy (n, evecp(j, i), basis(1, j), 1, evecfv(1, i), &
           & 1)
		!evecfv(:,i)=evecfv(:,i)+basis(:,j)*evecp(j,i)
         End Do
      End Do
!
End Subroutine
