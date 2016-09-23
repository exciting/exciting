!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: sdelta_mp
! !INTERFACE:
Real (8) Function sdelta_mp (n, x)
! !INPUT/OUTPUT PARAMETERS:
!   n : order (in,integer)
!   x : real argument (in,real)
! !DESCRIPTION:
!   Returns the smooth approximation to the Dirac delta function of order $N$
!   given by Methfessel and Paxton, {\it Phys. Rev. B} {\bf 40}, 3616 (1989),
!   $$ \tilde\delta(x)=\sum_{i=0}^N \frac{(-1)^i}{i!4^n\sqrt\pi} H_{2i}(x)
!    e^{-x^2},$$
!   where $H_j$ is the $j$th-order Hermite polynomial. This function has the
!   property
!   $$ \int_{-\infty}^{\infty}\tilde\delta(x)P(x)=P(0), $$
!   where $P(x)$ is any polynomial of degree $2N+1$ or less. The case $N=0$
!   corresponds to Gaussian smearing. This procedure is numerically stable
!   and accurate to near machine precision for $N\le 10$.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Real (8), Intent (In) :: x
! local variables
      Integer :: i
      Real (8), Parameter :: sqpi = 1.7724538509055160273d0
      Real (8) :: sum, t1
! external functions
      Real (8) :: factnm, hermite
      External factnm, hermite
      If (n .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(sdelta_mp): n < 0 : ", I8)') n
         Write (*,*)
         Stop
      End If
      If (n .Gt. 10) Then
         Write (*,*)
         Write (*, '("Error(sdelta_mp): n out of range : ", I8)') n
         Write (*,*)
         Stop
      End If
      If (Abs(x) .Gt. 12.d0) Then
         sdelta_mp = 0.d0
         Return
      End If
      sum = 0.d0
      Do i = 0, n
         t1 = 1.d0 / (factnm(i, 1)*dble(4**i)*sqpi)
         If (Mod(i, 2) .Ne. 0) t1 = - t1
         sum = sum + t1 * hermite (2*i, x) * Exp (-x**2)
      End Do
      sdelta_mp = sum
      Return
End Function
!EOC
