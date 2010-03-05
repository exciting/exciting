
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine exactupdatevectors (n, iunconverged, hamilton, overlap, r, &
& rhizvalue, eigenvector, trialvecs)
!calculate update equation with linsolver
!solvefor dA:  dA=(H-e*S)\R
! dA 	Update step to zero residual
! H 	Hamilton
! S 	Overlap
! e 	Rhitz Value
! R 	Residual
! trialvecs=eigenvector+dA
      Use modfvsystem
      Use modmain, Only: zone
      Integer, Intent (In) :: n, iunconverged
      Type (HermiteanMatrix), Intent (In) :: hamilton, overlap
      Complex (8), Intent (In) :: r (n, iunconverged), eigenvector (n, &
     & iunconverged)
      Real (8), Intent (In) :: rhizvalue (iunconverged)
      Complex (8), Intent (Out) :: trialvecs (n, iunconverged)
      Complex (8) :: dA (n, iunconverged), sigma
      Type (HermiteanMatrix) :: HES
      Call zcopy (n*iunconverged, eigenvector(1, 1), 1, trialvecs(1, &
     & 1), 1)
      Call zcopy (n*iunconverged, r(1, 1), 1, dA(1, 1), 1)
      Call newmatrix (HES, .False., n)
      Do i = 1, iunconverged
         Call HermiteanMatrixcopy (hamilton, HES)
         sigma = dcmplx (-rhizvalue(i), 0.0)
         Call HermiteanMatrixAXPY (sigma, overlap, HES)
         Call HermiteanmatrixLU (HES)
         Call Hermiteanmatrixlinsolve (HES, dA(:, i))
         Call zaxpy (n,-zone, dA(1, i), 1, trialvecs(1, i), 1)
      End Do
      Call deletematrix (HES)
End Subroutine
