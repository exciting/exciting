!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: zmatinp
! !INTERFACE:
!
!
Subroutine zmatinp (tapp, n, alpha, x, y, v, a)
! !INPUT/OUTPUT PARAMETERS:
!   tapp  : .true. if the matrix is to be applied to the input vector v,
!           .false. if the full matrix is to be calculated (in,logical)
!   n     : length of vectors (in,integer)
!   alpha : complex constant (in,complex)
!   x     : first input vector (in,complex(n))
!   y     : second input vector (in,complex(n))
!   v     : input vector to which matrix is applied if tapp is .true., otherwise
!           not referenced (in,complex(n))
!   a     : matrix applied to v if tapp is .true., otherwise the full matrix in
!           packed form (inout,complex(n+(n-1)*n/2))
! !DESCRIPTION:
!   Performs the rank-2 operation
!   $$ A_{ij}\rightarrow\alpha{\bf x}_i^*{\bf y}_j+\alpha^*{\bf y}_i^*{\bf x}_j
!    +A_{ij}, $$
!   where $A$ is stored in packed form. This is similar to the {\tt BLAS}
!   routine {\tt zhpr2}, except that here a matrix of inner products is formed
!   instead of an outer product of vectors. If {\tt tapp} is {\tt .true.} then
!   the matrix is applied to an input vector, rather than calculated explicitly.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: tapp
      Integer, Intent (In) :: n
      Complex (8), Intent (In) :: alpha
      Complex (8), Intent (In) :: x (n)
      Complex (8), Intent (In) :: y (n)
      Complex (8), Intent (In) :: v (n)
      Complex (8), Intent (Inout) :: a (*)
! local variables
      Integer :: i, j, k
      Real (8) :: a1, a2
! numbers less than eps are considered to be zero
      Real (8), Parameter :: eps = 1.d-12
      Complex (8) zt1, zt2
      If (tapp) Then
!--------------------------!
!     apply the matrix     !
!--------------------------!
         zt1 = y (1) * v (1)
         zt2 = x (1) * v (1)
         Do j = 2, n
            zt1 = zt1 + y (j) * v (j)
            zt2 = zt2 + x (j) * v (j)
         End Do
         zt1 = conjg (alpha*zt1)
         zt2 = alpha * conjg (zt2)
         If ((Abs(aimag(zt1)) .Gt. eps) .Or. (Abs(aimag(zt2)) .Gt. &
        & eps)) Then
! complex prefactors
            Do i = 1, n
               a (i) = a (i) + conjg (zt1*x(i)+zt2*y(i))
            End Do
         Else
! real prefactors
            a1 = dble (zt1)
            a2 = dble (zt2)
            If ((Abs(a1) .Gt. eps) .Or. (Abs(a2) .Gt. eps)) Then
               Do i = 1, n
                  a (i) = a (i) + conjg (a1*x(i)+a2*y(i))
               End Do
            End If
         End If
      Else
!---------------------------------------!
!     calculate the matrix elements     !
!---------------------------------------!
         k = 0
         Do j = 1, n
            If ((Abs(dble(x(j))) .Gt. eps) .Or. (Abs(aimag(x(j))) .Gt. &
           & eps) .Or. (Abs(dble(y(j))) .Gt. eps) .Or. &
           & (Abs(aimag(y(j))) .Gt. eps)) Then
               zt1 = conjg (alpha*y(j))
               zt2 = alpha * conjg (x(j))
               If ((Abs(aimag(zt1)) .Gt. eps) .Or. (Abs(aimag(zt2)) &
              & .Gt. eps)) Then
! complex prefactors
                  Do i = 1, j - 1
                     k = k + 1
                     a (k) = a (k) + conjg (zt1*x(i)+zt2*y(i))
                  End Do
                  k = k + 1
                  a (k) = dble (a(k)) + 2.d0 * dble (zt1*x(j))
               Else
! real prefactors
                  a1 = dble (zt1)
                  a2 = dble (zt2)
                  Do i = 1, j - 1
                     k = k + 1
                     a (k) = a (k) + conjg (a1*x(i)+a2*y(i))
                  End Do
                  k = k + 1
                  a (k) = dble (a(k)) + 2.d0 * a1 * dble (x(j))
               End If
            Else
               k = k + j
            End If
         End Do
      End If
      Return
End Subroutine
!EOC
