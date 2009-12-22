!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: sumrule
! !INTERFACE:
!
!
Subroutine sumrule (dynq)
! !INPUT/OUTPUT PARAMETERS:
      Use modinput
!   dynq : dynamical matrices on q-point set (in,real(3*natmtot,3*natmtot,nqpt))
! !DESCRIPTION:
!   Applies the same correction to all the dynamical matrices such that the
!   matrix for ${\bf q}=0$ satisfies the acoustic sum rule. In other words, the
!   matrices are updated with
!   $$ D_{ij}^{\bf q}\rightarrow D_{ij}^{\bf q}-\sum_{k=1}^3 \omega_k^0
!    v_{k;i}^0 v_{k;j}^0 $$
!   for all ${\bf q}$, where $\omega_k^0$ is the $k$th eigenvalue of the
!   ${\bf q}=0$ dynamical matrix and $v_{k;i}^0$ the $i$th component of its
!   eigenvector. The eigenvalues are assumed to be arranged in ascending order.
!   This ensures that the ${\bf q}=0$ dynamical matrix has 3 zero eigenvalues,
!   which the uncorrected matrix may not have due to the finite
!   exchange-correlation grid.
!
! !REVISION HISTORY:
!   Created May 2005 (JKD)
!EOP
!BOC
      Use modmain
      Implicit None
! arguments
      Complex (8), Intent (Inout) :: dynq (3*natmtot, 3*natmtot, nqpt)
! local variables
      Integer :: n, iq, iq0, i, j, k
      Integer :: lwork, info
      Real (8) :: t1
! allocatable arrays
      Real (8), Allocatable :: w (:)
      Real (8), Allocatable :: rwork (:)
      Complex (8), Allocatable :: work (:)
      Complex (8), Allocatable :: ev (:, :)
      n = 3 * natmtot
      Allocate (w(n))
      Allocate (rwork(3*n))
      lwork = 2 * n
      Allocate (work(lwork))
      Allocate (ev(n, n))
! find the Gamma-point phonon
      Do iq0 = 1, nqpt
         t1 = Sqrt (vql(1, iq0)**2+vql(2, iq0)**2+vql(3, iq0)**2)
         If (t1 .Lt. input%structure%epslat) Go To 10
      End Do
      Write (*,*)
      Write (*, '("Error(sumrule): no zero length q-vector")')
      Write (*,*)
      Stop
10    Continue
! compute the eigenvalues and vectors of the q=0 dynamical matrix
      Do i = 1, n
         Do j = i, n
            ev (i, j) = 0.5d0 * (dynq(i, j, iq0)+conjg(dynq(j, i, &
           & iq0)))
         End Do
      End Do
      Call zheev ('V', 'U', n, ev, n, w, work, lwork, rwork, info)
! subtract outer products of 3 lowest eigenvectors for q=0 from all the
! dynamical matrices
      Do iq = 1, nqpt
         Do i = 1, n
            Do j = 1, n
               Do k = 1, 3
                  dynq (i, j, iq) = dynq (i, j, iq) - w (k) * ev (i, k) &
                 & * conjg (ev(j, k))
               End Do
            End Do
         End Do
      End Do
      Deallocate (w, rwork, work, ev)
      Return
End Subroutine
!EOC
