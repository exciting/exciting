
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine genscclieff (iqr, nmax, n, scieff)
      Use modmain
      Use modxs
      Implicit None
  ! arguments
      Integer, Intent (In) :: iqr, n, nmax
      Complex (8), Intent (Out) :: scieff (nmax, nmax)
  ! local variables
      Logical :: tq0
      Complex (8), Allocatable :: scrn (:, :), scrnw (:, :, :), scrnh &
     & (:, :)
      Logical, External :: tqgamma
      Allocate (scrn(n, n), scrnw(n, 2, 3), scrnh(3, 3))
  ! read screening from file
      Call getscreen (iqr, n, scrnh, scrnw, scrn)
      tq0 = tqgamma (iqr)
      If (tq0) Then
     ! averaging using Lebedev-Laikov spherical grids
         Call angavsc0 (n, nmax, scrnh, scrnw, scrn, scieff)
      Else
     ! averaging using numerical method and extrapolation
         Call avscq (iqr, n, nmax, scrn, scieff)
      End If
End Subroutine genscclieff
