!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: getngkmax
! !INTERFACE:
!
!
Subroutine getngkmax
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Determines the largest number of ${\bf G+k}$-vectors with length less than
!   {\tt gkmax} over all the $k$-points and stores it in the global variable
!   {\tt ngkmax}. This variable is used for allocating arrays.
!
! !REVISION HISTORY:
!   Created October 2004 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ispn, ik, i, ig
      Real (8) :: v1 (3), v2 (3), t1, t2
      t1 = gkmax ** 2
      ngkmax = 0
      Do ispn = 1, nspnfv
         Do ik = 1, nkpt
            If (isspinspiral()) Then
! spin-spiral case
               If (ispn .Eq. 1) Then
                  v1 (:) = vkc (:, ik) + 0.5d0 * vqcss (:)
               Else
                  v1 (:) = vkc (:, ik) - 0.5d0 * vqcss (:)
               End If
            Else
               v1 (:) = vkc (:, ik)
            End If
            i = 0
            Do ig = 1, ngvec
               v2 (:) = vgc (:, ig) + v1 (:)
               t2 = v2 (1) ** 2 + v2 (2) ** 2 + v2 (3) ** 2
               If (t2 .Lt. t1) i = i + 1
            End Do
            ngkmax = Max (ngkmax, i)
         End Do
      End Do
      Return
End Subroutine
!EOC
