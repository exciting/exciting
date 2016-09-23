!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genwiq2
! !INTERFACE:
!
!
Subroutine genwiq2
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   The Fock matrix elements
!   $$ V^{\rm NL}_{ij{\bf k}}\equiv\sum_{l{\bf k'}}\int
!    \frac{\Psi^{\dag}_{i{\bf k}}({\bf r})\cdot\Psi_{l{\bf k}'}({\bf r})
!    \Psi^{\dag}_{l{\bf k}'}({\bf r}')\cdot\Psi_{j{\bf k}}({\bf r}')}
!    {|{\bf r}-{\bf r'}|}\,d{\bf r}\,d{\bf r'} $$
!   contain a divergent term in the sum over ${\bf k}'$ which behaves as
!   $1/q^2$, where ${\bf q}\equiv{\bf k}-{\bf k}'$ is in the first Brillouin
!   zone. The resulting convergence with respect to the number of discrete
!   $q$-points is very slow. This routine computes the weights
!   \begin{align}\label{genwiq2_1}
!    w_{{\bf q}_i}\equiv\int_{V_i}\frac{1}{q^2}\,d{\bf q}\;,
!   \end{align}
!   where the integral is over the small parallelepiped centered on ${\bf q}_i$,
!   so that integrals over the first Brillouin zone of the form
!   $$ I=\int_{\rm BZ}\frac{f({\bf q})}{q^2}\,d{\bf q}\;, $$
!   can be approximated by the sum
!   $$ I\approx\sum_i w_{{\bf q}_i}f({\bf q}_i) $$
!   which converges rapidly with respect to the number of $q$-points for smooth
!   functions $f$. The integral in (\ref{genwiq2_1}) is determined by evaluating
!   it numerically on increasingly finer grids and extrapolating to the
!   continuum. Agreement with Mathematica to at least 10 significant figures.
!
! !REVISION HISTORY:
!   Created August 2004 (JKD,SS)
!EOP
!BOC
      Implicit None
! local variables
      Integer, Parameter :: np = 5
      Integer, Parameter :: ns0 = 10, nss = 20
      Integer :: ns, iq, i1, i2, i3, i, ip
      Real (8) :: d (3), dv, sum, t1, t2
      Real (8) :: v1 (3), v2 (3), v3 (3)
      Real (8) :: xa (np), ya (np), c (np)
! external functions
      Real (8) :: polynom
      External polynom
! allocate global wiq2 array
      If (allocated(wiq2)) deallocate (wiq2)
      Allocate (wiq2(ngridq(1)*ngridq(2)*ngridq(3)))
! begin loop over q-points, note that the vectors vqc are assumed to be in the
! first Brillouin zone
      Do iq = 1, nqpt
! loop over different subdivisions
         ns = ns0
         Do ip = 1, np
! subdivision vectors in lattice coordinates
            Do i = 1, 3
               d (i) = 1.d0 / (dble(input%groundstate%ngridk(i)*2*ns))
            End Do
! smallest volume element
            dv = ((twopi**3)/omega) * d (1) * d (2) * d (3)
! compute the integral of 1/q^2
            sum = 0.d0
            Do i1 = - ns, ns - 1
               t1 = dble (i1) * d (1)
               v1 (:) = vqc (:, iq) + t1 * bvec (:, 1)
               Do i2 = - ns, ns - 1
                  t1 = dble (i2) * d (2)
                  v2 (:) = v1 (:) + t1 * bvec (:, 2)
                  Do i3 = - ns, ns - 1
                     t1 = dble (i3) * d (3)
                     v3 (:) = v2 (:) + t1 * bvec (:, 3)
                     t2 = v3 (1) ** 2 + v3 (2) ** 2 + v3 (3) ** 2
                     If (t2 .Gt. 1.d-14) Then
                        sum = sum + 1.d0 / t2
                     End If
                  End Do
               End Do
            End Do
            sum = sum * dv
            xa (ip) = dv ** (1.d0/3.d0)
            ya (ip) = sum
! increment number of subdivisions
            ns = ns + nss
         End Do
! extrapolate the volume element to zero with a polynomial
         wiq2 (iq) = polynom (0, np, xa, ya, c, 0.d0)
      End Do
      Return
End Subroutine
!EOC
