!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rschroddme
! !INTERFACE:
!
!
Subroutine rschroddme (m, l, k, e, np, nr, r, vr, nn, p0, p1, q0, q1)
! !INPUT/OUTPUT PARAMETERS:
!   m   : order of energy derivative (in,integer)
!   l   : angular momentum quantum number (in,integer)
!   k   : quantum number k, zero if Dirac eqn. is not to be used (in,integer)
!   e   : energy (in,real)
!   np  : order of predictor-corrector polynomial (in,integer)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   nn  : number of nodes (out,integer)
!   p0  : m th energy derivative of P (out,real(nr))
!   p1  : radial derivative of p0 (out,real(nr))
!   q0  : m th energy derivative of Q (out,real(nr))
!   q1  : radial derivative of q0 (out,real(nr))
! !DESCRIPTION:
!   Finds the solution to the $m$th energy derivative of the scalar relativistic
!   radial Schr\"{o}dinger equation using the routine {\tt rschrodint}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
      Integer, Intent (In) :: m
      Integer, Intent (In) :: l
      Integer, Intent (In) :: k
      Real (8), Intent (In) :: e
      Integer, Intent (In) :: np
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: vr (nr)
      Integer, Intent (Out) :: nn
      Real (8), Intent (Out) :: p0 (nr)
      Real (8), Intent (Out) :: p1 (nr)
      Real (8), Intent (Out) :: q0 (nr)
      Real (8), Intent (Out) :: q1 (nr)
! local variables
      Integer :: im, kpa, ir
! fine-structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
      Real (8) :: rm
! allocatable arrays
      Real (8), Allocatable :: p0p (:)
      Real (8), Allocatable :: g0 (:), g1 (:)
      Real (8), Allocatable :: f0 (:), f1 (:)
      Real (8), Allocatable :: cf (:, :)
      If (nr .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(rschroddme): invalid nr : ", I8)') nr
         Write (*,*)
         Stop
      End If
      If ((m .Lt. 0) .Or. (m .Gt. 6)) Then
         Write (*,*)
         Write (*, '("Error(rschroddme): m out of range : ", I8)') m
         Write (*,*)
         Stop
      End If
      If (k .Eq. 0) Then
! use the scalar relativistic Schrodinger equation
         Allocate (p0p(nr))
         If (m .Eq. 0) Then
            Call rschrodint (m, l, e, np, nr, r, vr, nn, p0p, p0, p1, &
           & q0, q1)
         Else
            Do im = 0, m
               Call rschrodint (im, l, e, np, nr, r, vr, nn, p0p, p0, &
              & p1, q0, q1)
               p0p (:) = p0 (:)
            End Do
         End If
         Deallocate (p0p)
      Else
! use the Dirac equation
         Allocate (g0(nr), g1(nr))
         Allocate (f0(nr), f1(nr))
         Allocate (cf(3, nr))
         If (k .Eq. l) Then
            kpa = k
         Else If (k .Eq. l+1) Then
            kpa = - k
         Else
            Write (*,*)
            Write (*, '("Error(rschroddme): incompatible l and k : ", 2&
           &I8)') l, k
            Write (*,*)
            Stop
         End If
         Call rdiracdme (m, kpa, e, np, nr, r, vr, nn, g0, g1, f0, f1)
! determine equivalent scalar relativistic functions
         Do ir = 1, nr
            rm = 1.d0 - 0.5d0 * (alpha**2) * vr (ir)
            p0 (ir) = g0 (ir)
            p1 (ir) = g1 (ir)
            q0 (ir) = (p1(ir)-p0(ir)/r(ir)) / (2.d0*rm)
         End Do
         Call fderiv (1, nr, r, q0, q1, cf)
         Deallocate (g0, g1, f0, f1, cf)
      End If
      Return
End Subroutine
!EOC
