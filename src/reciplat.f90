!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: reciplat
! !INTERFACE:
!
!
Subroutine reciplat
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Generates the reciprocal lattice vectors from the real-space lattice vectors
!   \begin{align*}
!     {\bf b}_1&=\frac{2\pi}{s}({\bf a}_2\times{\bf a}_3)\\
!     {\bf b}_2&=\frac{2\pi}{s}({\bf a}_3\times{\bf a}_1)\\
!     {\bf b}_3&=\frac{2\pi}{s}({\bf a}_1\times{\bf a}_2)
!   \end{align*}
!   and finds the unit cell volume $\Omega=|s|$, where
!   $s={\bf a}_1\cdot({\bf a}_2\times{\bf a}_3)$.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Real (8) :: t1
      Call r3cross (input%structure%crystal%basevect(:, 2), &
     & input%structure%crystal%basevect(:, 3), bvec(:, 1))
      Call r3cross (input%structure%crystal%basevect(:, 3), &
     & input%structure%crystal%basevect(:, 1), bvec(:, 2))
      Call r3cross (input%structure%crystal%basevect(:, 1), &
     & input%structure%crystal%basevect(:, 2), bvec(:, 3))
      t1 = input%structure%crystal%basevect(1, 1) * bvec (1, 1) + &
     & input%structure%crystal%basevect(2, 1) * bvec (2, 1) + &
     & input%structure%crystal%basevect(3, 1) * bvec (3, 1)
! compute unit cell volume
      omega = Abs (t1)
      If (omega .Lt. 1.d-6) Then
         Write (*,*)
         Write (*, '("Error(reciplat) omega too small : ", G18.10)') &
        & omega
         Write (*, '(" Lattice vectors may be collinear")')
         Write (*,*)
         Stop
      End If
      bvec (:, :) = (twopi/t1) * bvec (:, :)
      Return
End Subroutine
!EOC
