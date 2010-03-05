
! Copyright (C) 2002-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine hamiltonoverlapocopy_UL (system)
      Use modfvsystem
      Implicit None
      Type (evsystem) system
      Integer :: i, ohrank
      Complex (8), Allocatable :: tmp (:)
      Complex (8), Pointer :: h (:, :), o (:, :)
      If (ispacked(system%hamilton) .Eqv. .False.) Then
         Allocate (tmp(getrank(system%hamilton)))
         ohrank = getrank (system%hamilton)
         h => system%hamilton%za
         o => system%overlap%za
         Do i = 1, ohrank - 1
            Call zcopy (ohrank-i, h(i, i+1), ohrank, tmp, 1)
            tmp (1:ohrank-1) = conjg (tmp(1:ohrank-1))
            Call zcopy (ohrank-i, tmp, 1, h(i+1, i), 1)
         End Do
!
         Do i = 1, ohrank - 1
            Call zcopy (ohrank-i, o(i, i+1), ohrank, tmp, 1)
            tmp (1:ohrank-1) = conjg (tmp(1:ohrank-1))
            Call zcopy (ohrank-i, tmp, 1, o(i+1, i), 1)
         End Do
      End If
End Subroutine
