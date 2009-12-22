!
!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genjlgq0r (gq0, jlgq0r)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Real (8), Intent (In) :: gq0
      Real (8), Intent (Out) :: jlgq0r (0:input%groundstate%lmaxvr, &
     & nrcmtmax, nspecies)
! local variables
      Integer :: is, irc
      Real (8) :: t1
      Do is = 1, nspecies
         Do irc = 1, nrcmt (is)
            t1 = gq0 * rcmt (irc, is)
            Call sbessel (input%groundstate%lmaxvr, t1, jlgq0r(:, irc, &
           & is))
         End Do
      End Do
      Return
End Subroutine
