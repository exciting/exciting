
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Module for setting up the eigensystem
! it is designed in a way that all other subroutines
! dealing with setting up and solving the system can acsess the
! data transparently allowing to choose from different datatypes
! more easily
Module modfvsystem
      Implicit None
!
      Type HermiteanMatrix
!
         Integer :: rank
         Logical :: packed, ludecomposed
         Integer, Pointer :: ipiv (:)
         Complex (8), Pointer :: za (:, :), zap (:)
      End Type HermiteanMatrix
!
      Type evsystem
         Type (HermiteanMatrix) :: hamilton, overlap
      End Type evsystem
!
Contains
!
!
      Subroutine newmatrix (self, packed, rank)
         Type (HermiteanMatrix), Intent (Inout) :: self
         Logical, Intent (In) :: packed
         Integer, Intent (In) :: rank
         self%rank = rank
         self%packed = packed
         self%ludecomposed = .False.
         If (packed .Eqv. .True.) Then
            Allocate (self%zap(rank*(rank+1)/2))
            self%zap = 0.0
         Else
            Allocate (self%za(rank, rank))
            self%za = 0.0
         End If
      End Subroutine newmatrix
!
!
      Subroutine deletematrix (self)
         Type (HermiteanMatrix), Intent (Inout) :: self
         If (self%packed .Eqv. .True.) Then
            Deallocate (self%zap)
         Else
            Deallocate (self%za)
         End If
         If (self%ludecomposed) deallocate (self%ipiv)
      End Subroutine deletematrix
!
!
      Subroutine newsystem (self, packed, rank)
         Type (evsystem), Intent (Out) :: self
         Logical, Intent (In) :: packed
         Integer, Intent (In) :: rank
         Call newmatrix (self%hamilton, packed, rank)
         Call newmatrix (self%overlap, packed, rank)
      End Subroutine newsystem
!
!
      Subroutine deleteystem (self)
         Type (evsystem), Intent (Inout) :: self
         Call deletematrix (self%hamilton)
         Call deletematrix (self%overlap)
      End Subroutine deleteystem
!
!
      Subroutine Hermiteanmatrix_rank2update (self, n, alpha, x, y)
         Type (HermiteanMatrix), Intent (Inout) :: self
         Integer, Intent (In) :: n
         Complex (8), Intent (In) :: alpha, x (:), y (:)
!
         If (self%packed) Then
            Call ZHPR2 ('U', n, alpha, x, 1, y, 1, self%zap)
         Else
            Call ZHER2 ('U', n, alpha, x, 1, y, 1, self%za, self%rank)
         End If
      End Subroutine Hermiteanmatrix_rank2update
!
!
      Subroutine Hermiteanmatrix_indexedupdate (self, i, j, z)
         Type (HermiteanMatrix), Intent (Inout) :: self
         Integer :: i, j
         Complex (8) :: z
         Integer :: ipx
         If (self%packed .Eqv. .True.) Then
            ipx = ((i-1)*i) / 2 + j
            self%zap (ipx) = self%zap(ipx) + z
         Else
            If (j .Le. i) Then
               self%za (j, i) = self%za(j, i) + z
            Else
               Write (*,*) "warning lower part of hamilton updated"
            End If
         End If
         Return
      End Subroutine Hermiteanmatrix_indexedupdate
!
!
      Subroutine Hermiteanmatrixvector (self, alpha, vin, beta, vout)
         Implicit None
         Type (HermiteanMatrix), Intent (Inout) :: self
         Complex (8), Intent (In) :: alpha, beta
         Complex (8), Intent (Inout) :: vin (:)
         Complex (8), Intent (Inout) :: vout (:)
!
         If (self%packed .Eqv. .True.) Then
            Call zhpmv ("U", self%rank, alpha, self%zap, vin, 1, beta, &
           & vout, 1)
         Else
            Call zhemv ("U", self%rank, alpha, self%za, self%rank, vin, &
           & 1, beta, vout, 1)
         End If
      End Subroutine Hermiteanmatrixvector
!
!
      Function ispacked (self)
         Logical :: ispacked
         Type (HermiteanMatrix) :: self
         ispacked = self%packed
      End Function ispacked
!
!
      Function getrank (self)
         Integer :: getrank
         Type (HermiteanMatrix) :: self
         getrank = self%rank
      End Function getrank
!
!
      Subroutine HermiteanmatrixLU (self)
         Type (HermiteanMatrix) :: self
         Integer :: info
         If ( .Not. self%ludecomposed) allocate (self%ipiv(self%rank))
!
         If ( .Not. self%ludecomposed) Then
            If ( .Not. ispacked(self)) Then
               Call ZGETRF (self%rank, self%rank, self%za, self%rank, &
              & self%ipiv, info)
            Else
               Call ZHPTRF ('U', self%rank, self%zap, self%ipiv, info)
            End If
            If (info .Ne. 0) Then
               Write (*,*) "error in iterativearpacksecequn  Hermiteanm&
              &atrixLU ", info
               Stop
            End If
            self%ludecomposed = .True.
         End If
      End Subroutine HermiteanmatrixLU
!
!
      Subroutine Hermiteanmatrixlinsolve (self, b)
         Type (HermiteanMatrix) :: self
         Complex (8), Intent (Inout) :: b (:)
         Integer :: info
         If (self%ludecomposed) Then
            If ( .Not. ispacked(self)) Then
               Call ZGETRS ('N', self%rank, 1, self%za, self%rank, &
              & self%ipiv, b, self%rank, info)
            Else
               Call ZHPTRS ('U', self%rank, 1, self%zap, self%ipiv, b, &
              & self%rank, info)
            End If
            If (info .Ne. 0) Then
               Write (*,*) "error in iterativearpacksecequn Hermiteanma&
              &trixlinsolve ", info
               Stop
            End If
         End If
      End Subroutine Hermiteanmatrixlinsolve
!
!
      Subroutine HermiteanMatrixAXPY (alpha, x, y)
         Complex (8) :: alpha
         Type (HermiteanMatrix) :: x, y
         Integer :: mysize
         If (ispacked(x)) Then
            mysize = (x%rank*(x%rank+1)) / 2
            Call zaxpy (mysize, alpha, x%zap, 1, y%zap, 1)
         Else
            mysize = x%rank * (x%rank)
            Call zaxpy (mysize, alpha, x%za, 1, y%za, 1)
         End If
      End Subroutine HermiteanMatrixAXPY
!
!
      Subroutine HermiteanMatrixcopy (x, y)
         Complex (8) :: alpha
         Type (HermiteanMatrix) :: x, y
         Integer :: mysize
         If (ispacked(x)) Then
            mysize = (x%rank*(x%rank+1)) / 2
            Call zcopy (mysize, x%zap, 1, y%zap, 1)
         Else
            mysize = x%rank * (x%rank)
            Call zcopy (mysize, x%za, 1, y%za, 1)
         End If
      End Subroutine HermiteanMatrixcopy
!
!
      Subroutine HermiteanMatrixToFiles (self, prefix)
         Implicit None
         Type (HermiteanMatrix), Intent (In) :: self
         Character (256), Intent (In) :: prefix
         Character (256) :: filename
         If (ispacked(self)) Then
            filename = trim (prefix) // ".packed.real.OUT"
            Open (888, File=filename)
            Write (888,*) dble (self%zap)
         Else
            filename = trim (prefix) // ".real.OUT"
            Open (888, File=filename)
            Write (888,*) dble (self%za)
         End If
         Close (888)
!
         If (ispacked(self)) Then
            filename = trim (prefix) // ".packed.imag.OUT"
            Open (888, File=filename)
            Write (888,*) aimag (self%zap)
         Else
            filename = trim (prefix) // ".imag.OUT"
            Open (888, File=filename)
            Write (888,*) aimag (self%za)
         End If
         Close (888)
      End Subroutine HermiteanMatrixToFiles
!
!
      Subroutine HermiteanMatrixTruncate (self, threshold)
         Implicit None
         Type (HermiteanMatrix), Intent (Inout) :: self
         Real (8), Intent (In) :: threshold
         Integer :: n, i, j
         n = self%rank
         If (ispacked(self)) Then
            Do i = 1, n * (n+1) / 2
               If (Abs(dble(self%zap(i))) .Lt. threshold) self%zap(i) = &
              & self%zap(i) - dcmplx (dble(self%zap(i)), 0)
               If (Abs(aimag(self%zap(i))) .Lt. threshold) self%zap(i) &
              & = self%zap(i) - dcmplx (0, aimag(self%zap(i)))
            End Do
         Else
            Do j = 1, n
               Do i = 1, n
                  If (Abs(dble(self%za(i, j))) .Lt. threshold) &
                 & self%za(i, j) = self%za(i, j) - dcmplx &
                 & (dble(self%za(i, j)), 0)
                  If (Abs(aimag(self%za(i, j))) .Lt. threshold) &
                 & self%za(i, j) = self%za(i, j) - dcmplx (0, &
                 & aimag(self%za(i, j)))
               End Do
            End Do
         End If
      End Subroutine
!
!
      Subroutine HermiteanMatrixdiagonal (self, d)
         Implicit None
         Type (HermiteanMatrix), Intent (In) :: self
         Complex (8), Intent (Out) :: d (self%rank)
         Integer :: i
         If (ispacked(self)) Then
            Do i = 1, self%rank
               d (i) = self%zap((i*(i+1))/2)
            End Do
         Else
            Do i = 1, self%rank
               d (i) = self%za(i, i)
            End Do
         End If
      End Subroutine
!
End Module modfvsystem
