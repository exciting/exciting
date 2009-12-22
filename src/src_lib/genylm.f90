!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: genylm
! !INTERFACE:
!
!
Subroutine genylm (lmax, tp, ylm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   tp   : (theta, phi) coordinates (in,real(2))
!   ylm  : array of spherical harmonics (out,complex((lmax+1)**2))
! !DESCRIPTION:
!   Generates a sequence of spherical harmonics, including the Condon-Shortley
!   phase, evaluated at angles $(\theta,\phi)$ for $0<l<l_{\rm max}$. The values
!   are returned in a packed array {\tt ylm} indexed with $j=l(l+1)+m+1$. The
!   algorithm of Masters and Richards-Dinger is used, {\it Geophys. J. Int.}
!   {\bf 135}, 307 (1998). This routine is numerically stable and accurate to
!   near machine precision for $l\le 50$.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!   Improved stability, December 2005 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Real (8), Intent (In) :: tp (2)
      Complex (8), Intent (Out) :: ylm (*)
! local variables
      Integer :: l, m, lm1, lm2
      Real (8), Parameter :: fourpi = 12.566370614359172954d0
      Real (8) :: sn, cs, dx, sum, t1
! automatic arrays
      Real (8) :: x (0:lmax)
      Complex (8) z (lmax)
      If ((lmax .Lt. 0) .Or. (lmax .Gt. 50)) Then
         Write (*,*)
         Write (*, '("Error(genylm): lmax out of range : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      ylm (1) = 0.28209479177387814347d0
      If (lmax .Eq. 0) Return
      sn = Sin (tp(1))
      cs = Cos (tp(1))
! phase factors exp(i*m*phi)
      Do m = 1, lmax
         t1 = dble (m) * tp (2)
         z (m) = cmplx (Cos(t1), Sin(t1), 8)
      End Do
      Do l = 1, lmax
         If (Mod(l, 2) .Eq. 0) Then
            x (l) = 1.d0
         Else
            x (l) = - 1.d0
         End If
! recursion loop
         dx = 0.d0
         Do m = l, 1, - 1
            t1 = Sqrt (dble((l+m)*(l-m+1)))
            x (m-1) = - (sn*dx+dble(2*m)*cs*x(m)) / t1
            dx = sn * x (m) * t1
         End Do
! rescale values and multiply with phase factors
         t1 = sn
         sum = 0.d0
         Do m = 1, l
            x (m) = t1 * x (m)
            sum = sum + x (m) ** 2
            t1 = t1 * sn
         End Do
         sum = 2.d0 * sum + x (0) ** 2
         t1 = Sqrt (dble(2*l+1)/(fourpi*sum))
         lm1 = l * (l+1) + 1
         lm2 = lm1
         ylm (lm1) = t1 * x (0)
         Do m = 1, l
            lm1 = lm1 + 1
            lm2 = lm2 - 1
            ylm (lm1) = t1 * x (m) * z (m)
            ylm (lm2) = conjg (ylm(lm1))
            If (Mod(m, 2) .Ne. 0) ylm (lm2) = - ylm (lm2)
         End Do
      End Do
      Return
End Subroutine
!EOC
