
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine avscq (iqr, n, nmax, scrn, scieff)
      Use modmain
      Use modinput
      Use modxs
      Use invert
      Implicit None
  ! arguments
      Integer, Intent (In) :: iqr, n, nmax
      Complex (8), Intent (In) :: scrn (n, n)
      Complex (8), Intent (Out) :: scieff (nmax, nmax)
  ! local variables
      Integer :: iqrnr, j1, j2, flg
      Real (8) :: clwt
  ! find reduced q-point in non-reduced set
      iqrnr = iqmap (ivqr(1, iqr), ivqr(2, iqr), ivqr(3, iqr))
  ! invert dielectric matrix
      Call zinvert_hermitian (input%xs%BSE%scrherm, scrn, scieff(:n, &
     & :n))
      Do j1 = 1, n
         Do j2 = 1, j1
            If ((input%xs%BSE%sciavqhd .And. (j1 .Eq. 1) .And. (j2 .Eq. &
           & 1)) .Or. (input%xs%BSE%sciavqwg .And. (j1 .Ne. 1) .And. &
           & (j2 .Eq. 1)) .Or. (input%xs%BSE%sciavqwg .And. (j1 .Eq. 1) &
           & .And. (j2 .Ne. 1)) .Or. (input%xs%BSE%sciavqbd .And. (j1 &
           & .Ne. 1) .And. (j2 .Ne. 1))) Then
           ! numerical averaging on grids with extrapolation to continuum
               flg = 2
            Else
           ! analytic expression, no averaging
               flg = 0
            End If
        ! generate the (averaged) symmetrized Coulomb potential
            Call genwiqggp (flg, iqrnr, j1, j2, clwt)
        ! multiply with averaged Coulomb potential
            scieff (j1, j2) = scieff (j1, j2) * clwt
        ! set upper triangle
            scieff (j2, j1) = conjg (scieff(j1, j2))
         End Do
      End Do
End Subroutine avscq
