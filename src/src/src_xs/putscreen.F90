!
! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Subroutine putscreen (un, tq0, n, chi0, chi0h, chi0w)
      Use mod_constants, Only: krondelta
      Implicit None
  ! input parameters
      Logical, Intent (In) :: tq0
      Integer, Intent (In) :: un, n
      Complex (8), Intent (In) :: chi0 (n, n), chi0h (3, 3), chi0w (n, &
     & 2, 3)
  ! local variables
      Integer :: ig1, ig2, i, j
      Real (8) :: r1
      Do ig1 = 1, n
         Do ig2 = 1, n
            r1 = 0.d0
            If (ig1 .Eq. ig2) r1 = 1.d0
            If (tq0) Then
               If ((ig1 .Eq. 1) .And. (ig2 .Eq. 1)) Then
                  Write (un, '(2i8,3g18.10)') ((-i,-j, &
                 & dble(krondelta(i, j))-chi0h(i, j), &
                 & Abs(dble(krondelta(i, j))-chi0h(i, j)), j=1, 3), &
                 & i=1, 3)
               End If
               If ((ig1 .Eq. 1) .And. (ig2 .Ne. 1)) Then
                  Write (un, '(2i8,3g18.10)') (-i, ig2,-chi0w(ig2, 1, &
                 & i), Abs(-chi0w(ig2, 1, i)), i=1, 3)
               End If
               If ((ig1 .Ne. 1) .And. (ig2 .Eq. 1)) Then
                  Write (un, '(2i8,3g18.10)') (ig1,-j,-chi0w(ig1, 2, &
                 & j), Abs(-chi0w(ig1, 2, j)), j=1, 3)
               End If
               If ((ig1 .Ne. 1) .And. (ig2 .Ne. 1)) Then
                  Write (un, '(2i8,3g18.10)') ig1, ig2, r1 - chi0 (ig1, &
                 & ig2), Abs (r1-chi0(ig1, ig2))
               End If
            Else
               Write (un, '(2i8,3g18.10)') ig1, ig2, r1 - chi0 (ig1, &
              & ig2), Abs (r1-chi0(ig1, ig2))
            End If
         End Do
      End Do
End Subroutine
