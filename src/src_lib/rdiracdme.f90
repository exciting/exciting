!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rdiracdme
! !INTERFACE:
!
!
Subroutine rdiracdme (m, kpa, e, np, nr, r, vr, nn, g0, g1, f0, f1)
! !INPUT/OUTPUT PARAMETERS:
!   m   : order of energy derivative (in,integer)
!   kpa : quantum number kappa (in,integer)
!   e   : energy (in,real)
!   np  : order of predictor-corrector polynomial (in,integer)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   nn  : number of nodes (out,integer)
!   g0  : m th energy derivative of the major component multiplied by r
!         (out,real(nr))
!   g1  : radial derivative of g0 (out,real(nr))
!   f0  : m th energy derivative of the minor component multiplied by r
!         (out,real(nr))
!   f1  : radial derivative of f0 (out,real(nr))
! !DESCRIPTION:
!   Finds the solution to the $m$th energy derivative of the radial Dirac
!   equation using the routine {\tt rdiracint}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: m
      Integer, Intent (In) :: kpa
      Real (8), Intent (In) :: e
      Integer, Intent (In) :: np
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: vr (nr)
      Integer, Intent (Out) :: nn
      Real (8), Intent (Out) :: g0 (nr)
      Real (8), Intent (Out) :: g1 (nr)
      Real (8), Intent (Out) :: f0 (nr)
      Real (8), Intent (Out) :: f1 (nr)
! local variables
      Integer :: im
! automatic arrays
      Real (8) :: g0p (nr), f0p (nr)
      If (nr .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(rdiracdme): invalid nr : ", I8)') nr
         Write (*,*)
         Stop
      End If
      If ((m .Lt. 0) .Or. (m .Gt. 6)) Then
         Write (*,*)
         Write (*, '("Error(rdiracdme): m out of range : ", I8)') m
         Write (*,*)
         Stop
      End If
      If (m .Eq. 0) Then
         Call rdiracint (m, kpa, e, np, nr, r, vr, nn, g0p, f0p, g0, &
        & g1, f0, f1)
      Else
         Do im = 0, m
            Call rdiracint (im, kpa, e, np, nr, r, vr, nn, g0p, f0p, &
           & g0, g1, f0, f1)
            g0p (:) = g0 (:)
            f0p (:) = f0 (:)
         End Do
      End If
      Return
End Subroutine
!EOC
