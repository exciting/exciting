!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module invert
  use modmpi
      Implicit None
Contains
!
!
      Subroutine zinvert_lapack (m, mi)
         Implicit None
    ! arguments
         Complex (8), Intent (In) :: m (:, :)
         Complex (8), Intent (Out) :: mi (:, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'zinvert_lapack'
         Complex (8), Allocatable :: zwork (:)
         Integer, Allocatable :: ipiv (:)
         Integer :: lwork, info, sh (2), n
         sh = shape (m)
         n = sh (1)
         Allocate (ipiv(n))
         lwork = 2 * n
         Allocate (zwork(lwork))
         mi (:, :) = m (:, :)
         Call zgetrf (n, n, mi, n, ipiv, info)
         If (info .Ne. 0) Then
            Write (*,*)
            Write (*, '("Error(", a, "): zgetrf returned non-zero info : ", I8)') thisnam, info
            Write (*,*)
            Call terminate
         End If
         Call zgetri (n, mi, n, ipiv, zwork, lwork, info)
         If (info .Ne. 0) Then
            Write (*,*)
            Write (*, '("Error(", a, "): zgetri returned non-zero info : ", I8)') thisnam, info
            Write (*,*)
            Call terminate
         End If
         Deallocate (ipiv, zwork)
      End Subroutine zinvert_lapack
!
!
      Subroutine zinvert_hermitian (flag, m, mi)
         Implicit None
    ! arguments
         Integer, Intent (In) :: flag
         Complex (8), Intent (In) :: m (:, :)
         Complex (8), Intent (Out) :: mi (:, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'zinvert_hermitian'
         Character (1) :: uplo
         Integer :: info, n, j, sh (2)
         Complex (8), Allocatable :: tm (:, :)
    ! we do not check if both arguments have same shapes and are square matrices
         sh = shape (m)
         n = sh (1)
         Allocate (tm(sh(1), sh(2)))
         tm (:, :) = m (:, :)
         info = 0
         Select Case (flag)
         Case (0)
       ! invert full matrix (matrix is allowed to be not strictly Hermitian)
            Call zinvert_lapack (tm, mi)
         Case (1)
       ! Hermitian average matrix
            tm = 0.5d0 * (tm+conjg(transpose(tm)))
            uplo = 'u'
         Case (2)
       ! assume Hermitian and use upper triangle for inversion
            uplo = 'u'
         Case (3)
       ! assume Hermitian and use lower triangle for inversion
            uplo = 'l'
         Case Default
            Write (*,*)
            Write (*, '("Error(", a, "): not a valid flag:", i6)') trim &
           & (thisnam), flag
            Write (*,*)
            Call terminate
         End Select
         Select Case (flag)
         Case (1, 2, 3)
       ! set up unity matrix for zposv
            mi (:, :) = (0.d0, 0.d0)
            Forall (j=1:n)
               mi (j, j) = (1.d0, 0.d0)
            End Forall
       ! invert using upper/lower triangle of matrix
            Call zposv (uplo, n, n, tm, n, mi, n, info)
         End Select
         If (info .Ne. 0) Then
            Write (*,*)
            Write (*, '("Error(", a, "): zposv returned non-zero info :&
           & ", I8)') trim (thisnam), info
            Write (*,*)
            Call terminate
         End If
         Deallocate (tm)
      End Subroutine zinvert_hermitian
!
End Module invert
