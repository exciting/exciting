!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: sphcover
! !INTERFACE:
!
!
Subroutine sphcover (n, tp)
! !INPUT/OUTPUT PARAMETERS:
!   n  : number of required points (in,integer)
!   tp : (theta, phi) coordinates (out,real(2,n))
! !DESCRIPTION:
!   Produces a set of $N$ points which cover the unit sphere nearly optimally.
!   The points in $(\theta,\phi)$ coordinates are generated using the explicit
!   formula
!   \begin{align*}
!    \theta_k&=\arccos(h_k), \qquad h_k=\frac{2(k-1)}{N-1}-1, \qquad
!    1\le k \le N \\
!    \phi_k&=\left(\phi_{k-1}+C/\sqrt{N(1-h_k^2)}\right)({\rm mod}\;2\pi),
!    \qquad 2\le k\le N-1, \qquad \phi_1=\phi_N=0,
!   \end{align*}
!   where $C=(8\pi/\sqrt{3})^{1/2}$. See E. B. Saff and A. B. J. Kuijlaars,
!   {\it Math. Intell.} {\bf 19}, 5 (1997).
!
! !REVISION HISTORY:
!   Created April 2008 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Real (8), Intent (Out) :: tp (2, n)
! local variables
      Integer :: k
      Real (8), Parameter :: pi = 3.1415926535897932385d0
      Real (8), Parameter :: twopi = 6.2831853071795864769d0
      Real (8) :: c, h, t1, t2
      If (n .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(sphcover): n <= 0 : ", I8)') n
         Write (*,*)
         Stop
      End If
      c = Sqrt (8.d0*pi/Sqrt(3.d0))
      t1 = c / Sqrt (dble(n))
      tp (1, 1) = pi
      tp (2, 1) = 0.d0
      Do k = 2, n - 1
         h = dble (2*(k-1)) / dble (n-1) - 1.d0
         tp (1, k) = Acos (h)
         t2 = tp (2, k-1) + t1 / Sqrt (1.d0-h**2)
         tp (2, k) = Mod (t2, twopi)
      End Do
      tp (1, n) = 0.d0
      tp (2, n) = 0.d0
      Return
End Subroutine
!EOC
