!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_dyson
      Implicit None
Contains
!
!BOP
! !ROUTINE: dyson
! !INTERFACE:
!
!
      Subroutine dyson (n, s0, k, s)
! !USES:
         Use invert
! !INPUT/OUTPUT PARAMETERS:
!   n     : matrix size of local field effects (in,integer)
!   s0    : S0 matrix (in,complex(:,:))
!   k     : kernel matrix (in,complex(:,:))
!   s     : S (solution) matrix (in,complex(:,:))
! !DESCRIPTION:
!   Solve Dyson's equation
!     $$   S = S_0 + S_0 K S  $$
!   for $S$ by inversion;
!     $$ S = \left[ 1 + S_0 K \right]^{-1} S_0. $$
!   The inversion is carried out using the LAPACK routines {\tt zgetrf} and
!   {\tt zgetri}.
!
! !REVISION HISTORY:
!   Created March 2005 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Intent (In) :: n
         Complex (8), Intent (In) :: s0 (:, :), k (:, :)
         Complex (8), Intent (Out) :: s (:, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'dyson'
         Complex (8), Parameter :: zone = (1.d0, 0.d0), zzero = (0.d0, &
        & 0.d0)
         Complex (8), Allocatable :: mt (:, :)
         Integer :: shs0 (2), shk (2), shs (2), nmin, nmax, j
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
         Allocate (mt(n, n))
!
    ! calculate matrix -S0*K
         Call zgemm ('n', 'n', n, n, n,-zone, s0, n, k, n, zzero, mt, &
        & n)
!
    ! calculate matrix T=[1 - S0*K]
         Forall (j=1:n)
           mt (j, j) = mt (j, j) + 1.d0
         end forall
!
    ! invert matrix T
         Call zinvert_lapack (mt, mt)
!
    ! calculate matrix S=T^-1*S0
         Call zgemm ('n', 'n', n, n, n, zone, mt, n, s0, n, zzero, s, &
        & n)
!
         Deallocate (mt)
!
      End Subroutine dyson
!EOC
!
End Module m_dyson
