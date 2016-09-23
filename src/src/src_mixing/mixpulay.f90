
! Copyright (C) 2008 S. Suehara. This file is distributed under the terms of the
! GNU Lesser General Public License. See the file COPYING for license details.

!BOP
! !ROUTINE: mixpulay
! !INTERFACE:
!
!
Subroutine mixpulay (iscl, n, maxsd, nu, mu, f, d)
! !INPUT/OUTPUT PARAMETERS:
!   iscl  : self-consistent loop number (in,integer)
!   n     : vector length (in,integer)
!   maxsd : maximum subspace dimension (in,integer)
!   nu    : current output vector as well as the next input vector of the
!           self-consistent loop (inout,real(n))
!   mu    : used internally (inout,real(n,maxsd))
!   f     : used internally (inout,real(n,maxsd))
!   d     : RMS difference between old and new output vectors (out,real)
! !DESCRIPTION:
!   Pulay's mixing scheme which uses direct inversion in the iterative subspace
!   (DIIS). See {\it Chem. Phys. Lett.} {\bf 73}, 393 (1980).
!
! !REVISION HISTORY:
!   Created October 2008 (S. Suehara, NIMS)
!   Modified, October 2008 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: iscl
      Integer, Intent (In) :: n
      Integer, Intent (In) :: maxsd
      Real (8), Intent (Inout) :: nu (n)
      Real (8), Intent (Inout) :: mu (n, maxsd)
      Real (8), Intent (Inout) :: f (n, maxsd)
      Real (8), Intent (Out) :: d
! local variables
      Integer :: i, j, k, m, jc, jn, info
! initial mixing parameter
      Real (8), Parameter :: beta = 0.1d0
! allocatable arrays
      Integer, Allocatable :: ipiv (:)
      Real (8), Allocatable :: alpha (:), a (:, :), work (:)
! external functions
      Real (8) :: ddot
      External ddot
      If (n .Lt. 1) Then
         Write (*,*)
         Write (*, '("Error(mixpulay): n < 1 : ", I8)') n
         Write (*,*)
         Stop
      End If
      If (maxsd .Lt. 2) Then
         Write (*,*)
         Write (*, '("Error(mixpulay): maxsd < 2 : ", I8)') maxsd
         Write (*,*)
         Stop
      End If
      If (iscl .Le. 1) Then
         mu (:, 1) = nu (:)
         f (:, 1) = 0.d0
         d = 1.d0
         Return
      End If
! current index
      jc = Mod (iscl-1, maxsd) + 1
! next index
      jn = Mod (iscl, maxsd) + 1
      If (iscl .Le. 2) Then
         nu (:) = beta * nu (:) + (1.d0-beta) * mu (:, 1)
         f (:, 2) = nu (:) - mu (:, 1)
         mu (:, 2) = nu (:)
         If (maxsd .Ge. 3) mu (:, 3) = 0.d0
         d = 0.d0
         Do k = 1, n
            d = d + f (k, 2) ** 2
         End Do
         d = Sqrt (d/dble(n))
         Return
      End If
! matrix size
      m = Min (iscl, maxsd) + 1
      Allocate (ipiv(m), alpha(m), a(m, m), work(m))
! compute f and RMS difference
      d = 0.d0
      Do k = 1, n
         f (k, jc) = nu (k) - mu (k, jc)
         d = d + f (k, jc) ** 2
      End Do
      d = Sqrt (d/dble(n))
! solve the linear system
      a (:, :) = 0.d0
      Do i = 1, m - 1
         Do j = i, m - 1
            a (i, j) = a (i, j) + ddot (n, f(:, i), 1, f(:, j), 1)
         End Do
         a (i, m) = 1.d0
      End Do
      alpha (:) = 0.d0
      alpha (m) = 1.d0
      Call dsysv ('U', m, 1, a, m, ipiv, alpha, m, work, m, info)
      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(mixpulay): could not solve linear system")')
         Write (*, '(" DSYSV returned INFO = ", I8)') info
         Write (*,*)
         Stop
      End If
      nu (:) = 0.d0
      Do i = 1, m - 1
         nu (:) = nu (:) + alpha (i) * (mu(:, i)+f(:, i))
      End Do
      mu (:, jn) = nu (:)
      Deallocate (ipiv, alpha, a, work)
      Return
End Subroutine
!EOC
