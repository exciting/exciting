!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: fderiv
! !INTERFACE:
!
!
Subroutine fderiv (m, n, x, f, g, cf)
! !INPUT/OUTPUT PARAMETERS:
!   m  : order of derivative (in,integer)
!   n  : number of points (in,integer)
!   x  : abscissa array (in,real(n))
!   f  : function array (in,real(n))
!   g  : (anti-)derivative of f (out,real(n))
!   cf : spline coefficients (out,real(3,n))
! !DESCRIPTION:
!   Given function $f$ defined on a set of points $x_i$ then if $m\ge 0$ this
!   routine computes the $m$th derivative of $f$ at each point. If $m<0$ the
!   anti-derivative of $f$ given by
!   $$ g(x_i)=\int_{x_1}^{x_i} f(x)\,dx $$
!   is calculated by fitting the function to a clamped cubic spline. See routine
!   {\tt spline}.
!
! !REVISION HISTORY:
!   Created May 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: m
      Integer, Intent (In) :: n
      Real (8), Intent (In) :: x (n)
      Real (8), Intent (In) :: f (n)
      Real (8), Intent (Out) :: g (n)
      Real (8), Intent (Out) :: cf (3, n)
! local variables
      Integer :: i
      Real (8) :: dx
      If (n .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(fderiv): invalid number of points : ", I8)') n
         Write (*,*)
         Stop
      End If
      If (m .Eq. 0) Then
         g (:) = f (:)
         Return
      End If
      If (m .Ge. 4) Then
         g (:) = 0.d0
         Return
      End If
! high accuracy (anti-)derivatives from a clamped spline fit to the data
      Call spline (n, x, 1, f, cf)
      Select Case (m)
      Case (:-1)
         g (1) = 0.d0
         Do i = 1, n - 1
            dx = x (i+1) - x (i)
            g (i+1) = g (i) + (((0.25d0*cf(3, &
           & i)*dx+0.3333333333333333333d0*cf(2, i))*dx+0.5d0*cf(1, &
           & i))*dx+f(i)) * dx
         End Do
      Case (1)
         g (:) = cf (1, :)
      Case (2)
         g (:) = 2.d0 * cf (2, :)
      Case (3)
         g (:) = 6.d0 * cf (3, :)
      End Select
      Return
End Subroutine
!EOC
