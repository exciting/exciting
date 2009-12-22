!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rdirac
! !INTERFACE:
!
!
Subroutine rdirac (n, l, k, np, nr, r, vr, eval, g0, f0)
! !INPUT/OUTPUT PARAMETERS:
!   n    : principal quantum number (in,integer)
!   l    : quantum number l (in,integer)
!   k    : quantum number k (l or l+1) (in,integer)
!   np   : order of predictor-corrector polynomial (in,integer)
!   nr   : number of radial mesh points (in,integer)
!   r    : radial mesh (in,real(nr))
!   vr   : potential on radial mesh (in,real(nr))
!   eval : eigenvalue without rest-mass energy (inout,real)
!   g0   : major component of the radial wavefunction (out,real(nr))
!   f0   : minor component of the radial wavefunction (out,real(nr))
! !DESCRIPTION:
!   Finds the solution to the radial Dirac equation for a given potential $v(r)$
!   and quantum numbers $n$, $k$ and $l$. The method involves integrating the
!   equation using the predictor-corrector method and adjusting $E$ until the
!   number of nodes in the wavefunction equals $n-l-1$. The calling routine must
!   provide an initial estimate for the eigenvalue. Note that the arrays
!   {\tt g0} and {\tt f0} represent the radial functions multiplied by $r$.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Integer, Intent (In) :: l
      Integer, Intent (In) :: k
      Integer, Intent (In) :: np
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: vr (nr)
      Real (8), Intent (Inout) :: eval
      Real (8), Intent (Out) :: g0 (nr)
      Real (8), Intent (Out) :: f0 (nr)
! local variables
      Integer, Parameter :: maxit = 2000
      Integer :: kpa, it, nn, ir, irm, nnd, nndp
! energy convergence tolerance
      Real (8), Parameter :: eps = 1.d-11
      Real (8) :: t1, de
! automatic arrays
      Real (8) :: g1 (nr), f1 (nr), fr (nr), gr (nr), cf (3, nr)
      If (k .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(rdirac): k <= 0 : ", I8)') k
         Write (*,*)
         Stop
      End If
      If (k .Gt. n) Then
         Write (*,*)
         Write (*, '("Error(rdirac): incompatible n and k : ", 2I8)') &
        & n, k
         Write (*,*)
         Stop
      End If
      If ((k .Eq. n) .And. (l .Ne. k-1)) Then
         Write (*,*)
         Write (*, '("Error(rdirac): incompatible n, k and l : ", 3I8)') n, k, l
         Write (*,*)
         Stop
      End If
      If (k .Eq. l) Then
         kpa = k
      Else If (k .Eq. l+1) Then
         kpa = - k
      Else
         Write (*,*)
         Write (*, '("Error(rdirac): incompatible l and k : ", 2I8)') &
        & l, k
         Write (*,*)
         Stop
      End If
      de = 1.d0
      nndp = 0
      Do it = 1, maxit
! integrate the Dirac equation
         Call rdiracdme (0, kpa, eval, np, nr, r, vr, nn, g0, g1, f0, &
        & f1)
! check the number of nodes
         nnd = nn - (n-l-1)
         If (nnd .Gt. 0) Then
            eval = eval - de
         Else
            eval = eval + de
         End If
         If (it .Gt. 1) Then
            If ((nnd .Ne. 0) .Or. (nndp .Ne. 0)) Then
               If (nnd*nndp .Le. 0) Then
                  de = de * 0.5d0
               Else
                  de = de * 1.1d0
               End If
            End If
         End If
         nndp = nnd
         If (de .Lt. eps*(Abs(eval)+1.d0)) Go To 20
      End Do
      Write (*,*)
      Write (*, '("Error(rdirac): maximum iterations exceeded")')
      Write (*,*)
      Stop
20    Continue
! find effective infinity and set wavefunction to zero after that point
! major component
      irm = nr
      Do ir = 2, nr
         If ((g0(ir-1)*g0(ir) .Lt. 0.d0) .Or. (g1(ir-1)*g1(ir) .Lt. &
        & 0.d0)) irm = ir
      End Do
      g0 (irm:nr) = 0.d0
! minor component
      irm = nr
      Do ir = 2, nr
         If ((f0(ir-1)*f0(ir) .Lt. 0.d0) .Or. (f1(ir-1)*f1(ir) .Lt. &
        & 0.d0)) irm = ir
      End Do
      f0 (irm:nr) = 0.d0
! normalise
      Do ir = 1, nr
         fr (ir) = g0 (ir) ** 2 + f0 (ir) ** 2
      End Do
      Call fderiv (-1, nr, r, fr, gr, cf)
      t1 = Sqrt (Abs(gr(nr)))
      If (t1 .Gt. 0.d0) Then
         t1 = 1.d0 / t1
      Else
         Write (*,*)
         Write (*, '("Error(rdirac): zero wavefunction")')
         Write (*,*)
         Stop
      End If
      g0 (:) = t1 * g0 (:)
      f0 (:) = t1 * f0 (:)
      Return
End Subroutine
!EOC
