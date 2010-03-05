
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine precondspectrumupdate (n, m, hamilton, overlap, P, w)
      Use modmain, Only: nstfv, nmatmax
      Use diisinterfaces
      Implicit None
      Integer, Intent (In) :: n, m
      Complex (8), Intent (In) :: hamilton (n, n), overlap (n, n)
      Complex (8), Intent (In) :: P (nmatmax, nmatmax)
      Real (8), Intent (Inout) :: w (nmatmax)
      Complex (8) :: h (n, m), s (n, m)
End Subroutine precondspectrumupdate
