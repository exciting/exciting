!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: getngqmax
! !INTERFACE:
!
!
Subroutine getngqmax
! !USES:
      Use modinput
      Use modmain
      use modmpi
      Use modxs
! !DESCRIPTION:
!   Determines the largest number of ${\bf G+q}$-vectors with length less than
!   {\tt gqmax} over all the ${\bf q}$-points and stores it in the global
!   variable {\tt ngqmax}. This variable is used for allocating arrays.
!   Based upon the routine {\tt getngkmax}.
!
! !REVISION HISTORY:
!   Created October 2006 (Sagmeister)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: iq, i, j, ig
      Real (8) :: v1 (3), v2 (3), t1, t2
      If (input%xs%gqmax .Lt. input%structure%epslat) Then
         intgqv (:, :) = 0
         ngqmax = 1
         ngridgq (:) = 1
         Return
      End If
      t1 = input%xs%gqmax ** 2
      intgqv (:, :) = 0
      ngqmax = 0
      Do iq = 1, nqpt
         v1 (:) = vqc (:, iq)
         i = 0
         Do ig = 1, ngvec
            v2 (:) = vgc (:, ig) + v1 (:)
            ! cutoff type for G-vectors
            if (tgqmaxg) then
              t2 = vgc(1,ig) ** 2 + vgc(2,ig) ** 2 + vgc(3,ig) ** 2
            else
              t2 = v2 (1) ** 2 + v2 (2) ** 2 + v2 (3) ** 2
            end if
            If (t2 .Lt. t1) Then
               i = i + 1
               Do j = 1, 3
                  intgqv (j, 1) = Min (intgqv(j, 1), ivg(j, ig))
                  intgqv (j, 2) = Max (intgqv(j, 2), ivg(j, ig))
               End Do
            End If
         End Do
         ngqmax = Max (ngqmax, i)
      End Do
      ngridgq (:) = intgqv (:, 2) - intgqv (:, 1) + 1
! debug output
      If (input%xs%dbglev .Gt. 1) Then
         Write (*, '(a)') 'Debug(getngqmax): intgqv:'
         Write (*, '(2i6)') intgqv (1, 1), intgqv (1, 2)
         Write (*, '(2i6)') intgqv (2, 1), intgqv (2, 2)
         Write (*, '(2i6)') intgqv (3, 1), intgqv (3, 2)
         Write (*, '(i8)') ngqmax
         Write (*,*)
      End If
      If (ngqmax .Lt. 1) Then
         Write (*,*)
         Write (*, '("Error(getngqmax): no G-vectors found - increase g&
        &qmax")')
         Write (*,*)
         Call terminate
      End If
End Subroutine getngqmax
!EOC
