!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_dysonsym
      Implicit None
Contains
!
!BOP
! !ROUTINE: dysonsym
! !INTERFACE:
!
!
      Subroutine dysonsym (n, s0, k, s)
! !USES:
         Use invert
         Use modxs
! !INPUT/OUTPUT PARAMETERS:
!   n     : matrix size of local field effects (in,integer)
!   s0    : S0 matrix (in,complex(:,:))
!   k     : kernel matrix multiplied by S0 from both sides (in,complex(:,:))
!   s     : S (solution) matrix (in,complex(:,:))
! !DESCRIPTION:
!   Solve symmetric form of Dyson's equation
!     $$   S = S_0 + S_0 (1 + S_0^{-1} K S_0^{-1}) S  $$
!   for $S$ by inversion;
!     $$ S = S_0\left[ S_0(1-S0) - T\right]^{-1} S_0. $$
!   The inversion is carried out using the LAPACK routines {\tt zgetrf} and
!   {\tt zgetri}.
!
! !REVISION HISTORY:
!   Created October 2008 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Intent (In) :: n
         Complex (8), Intent (In) :: s0 (:, :), k (:, :)
         Complex (8), Intent (Out) :: s (:, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'dysonsym'
         Complex (8), Parameter :: zone = (1.d0, 0.d0), zzero = (0.d0, &
        & 0.d0)
         Complex (8), Allocatable :: mt (:, :), mt2 (:, :)
         Integer :: shs0 (2), shk (2), shs (2), nmin, nmax, j
!
         Complex (8), Allocatable :: solv (:), s0row (:), s0col (:), u &
        & (:, :), vh (:, :), work (:)
         Real (8), Allocatable :: singv (:), rwork (:)
         Integer, Allocatable :: ipiv (:)
         Integer :: info, lwork
         Real (8) :: eps
!
         Complex (8), External :: zdotu
!
!
    ! check matrix sizes
         shs0 = shape (s0)
         shk = shape (k)
         shs = shape (s)
         nmin = minval ( (/ shs0, shk, shs /))
         nmax = maxval ( (/ shs0, shk, shs /))
         If ((nmin .Ne. nmax) .Or. (nmin .Lt. n)) Then
            Write (*, '("Error(", a, "): inconsistent matrix sizes")') &
           & trim (thisnam)
            Write (*, '("  n :", i9)') n
            Write (*, '("  S0:", 2i9)') shs0
            Write (*, '("  K :", 2i9)') shk
            Write (*, '("  S :", 2i9)') shs
            Call terminate
         End If
!
    ! allocate
         Allocate (mt(n, n), mt2(n, n))
!
    ! calculate matrix 1-S0
         mt (:, :) = zzero
         Forall (j=1:n) mt (j, j) = 1.d0
         mt (:, :) = mt (:, :) - s0 (:, :)
!
    ! calculate S0(1-S0)
         Call zgemm ('n', 'n', n, n, n, zone, s0, n, mt, n, zzero, mt2, &
        & n)
!
    ! calculate X := S0(1-S0) - K
         mt2 (:, :) = mt2 (:, :) - k (:, :)
!
!    ! calculate [S0(1-S0) - K]^-1 =: Y = X^-1
!    call zinvert_lapack(mt2,mt)
!
!    ! calculate S0 Y
!    call zgemm('n','n', n, n, n, zone, s0, n, mt, n, zzero, mt2, n )
!
!    ! calculate solution S = S0 Y S0 = S0 [S0(1-S0) - K]^(-1) S0
!    call zgemm('n','n', n, n, n, zone, mt2, n, s0, n, zzero, s, n )
!
!
!
!!!goto 100
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
    ! F. Sottile, PhD thesis, p. 167 (Appendix E)
    ! solve linear system of equations instead of direct inversion
!
    ! first column of S0 is RHS of system of equations
         Allocate (solv(n), s0row(n), s0col(n))
         s0row (:) = s0 (1, :)
         s0col (:) = s0 (:, 1)
         solv (:) = s0col (:)
!
    !------------------------------------ solve linear system of equations
         Allocate (ipiv(n))
         Call zgetrf (n, n, mt2, n, ipiv, info)
         If (info .Ne. 0) Then
            Write (*,*)
            Write (*, '("Error(", a, "): zgetrf returned non-zero info : ", I8)') thisnam, info
            Write (*,*)
            Call terminate
         End If
         Call zgetrs ('n', n, 1, mt2, n, ipiv, solv, n, info)
         If (info .Ne. 0) Then
            Write (*,*)
            Write (*, '("Error(", a, "): zgetrs returned non-zero info : ", I8)') thisnam, info
            Write (*,*)
            Call terminate
         End If
         Deallocate (ipiv)
    !------------------------------------
!
    ! calculate S_00 = S0 * Solv
         s (1, 1) = zdotu (n, solv, 1, s0row, 1)
         Deallocate (solv, s0row, s0col)
!!!100 continue
         Go To 200
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
         lwork = 3 * n
         Allocate (singv(n), u(n, n), vh(n, n), work(lwork), &
        & rwork(5*n))
!
    ! try SVD for inversion of matrix X = S0(1 - S0) - K
         Call ZGESVD ('a', 'a', n, n, mt2, n, singv, u, n, vh, n, work, &
        & lwork, rwork, info)
         If (info .Ne. 0) Then
            Write (*,*)
            Write (*, '("Error(", a, "): zgesvd returned non-zero info : ", I8)') thisnam, info
            Write (*,*)
            Call terminate
         End If
!
   ! invert singular values above cutoff
         eps = 1.d-3
         Do j = 1, n
            If (singv(j) .Lt. eps) Then
               singv (j) = 0.d0
            Else
               singv (j) = 1.d0 / singv (j)
            End If
	! multiply singular values with U-matrix
            u (:, j) = u (:, j) * singv (j)
         End Do
    ! inverse of matrix X^+:
         Call zgemm ('n', 'n', n, n, n, zone, u, n, vh, n, zzero, s, n)
         s = conjg (transpose(s))
!
    ! left and right multiply with S0
    ! calculate S0 Y
         Call zgemm ('n', 'n', n, n, n, zone, s0, n, s, n, zzero, mt2, &
        & n)
!
    ! calculate solution S = S0 Y S0 = S0 [S0(1-S0) - K]^(-1) S0
         Call zgemm ('n', 'n', n, n, n, zone, mt2, n, s0, n, zzero, s, &
        & n)
!
!
         Deallocate (mt, mt2, singv, u, vh, work, rwork)
200      Continue
      End Subroutine dysonsym
!EOC
!
End Module m_dysonsym
