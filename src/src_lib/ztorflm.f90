!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: ztorflm
! !INTERFACE:
!
!
Subroutine ztorflm (lmax, zflm, rflm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   zflm : coefficients of complex spherical harmonic expansion
!          (in,complex((lmax+1)**2)))
!   rflm : coefficients of real spherical harmonic expansion
!          (out,real((lmax+1)**2)))
! !DESCRIPTION:
!   Converts a real function, $z_{lm}$, expanded in terms of complex spherical
!   harmonics into a real spherical harmonic expansion, $r_{lm}$:
!   $$ r_{lm}=\begin{cases}\frac{1}{\sqrt{2}}\Re(z_{lm}+(-1)^m z_{l-m}) & m>0 \\
!    \frac{1}{\sqrt{2}}\Im(-z_{lm}+(-1)^m z_{l-m}) & m<0 \\
!    \Re(z_{lm}) & m=0 \end{cases}\;. $$
!   See routine {\tt genrlm}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Complex (8), Intent (In) :: zflm (*)
      Real (8), Intent (Out) :: rflm (*)
! local variables
      Integer :: l, m, lm1, lm2
! real constant 1/sqrt(2)
      Real (8), Parameter :: c1 = 0.7071067811865475244d0
      If (lmax .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(ztorflm): lmax < 0 : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      lm1 = 0
      Do l = 0, lmax
         lm2 = lm1 + 2 * (l+1)
         Do m = - l, l
            lm1 = lm1 + 1
            lm2 = lm2 - 1
            If (m .Gt. 0) Then
               If (Mod(m, 2) .Ne. 0) Then
                  rflm (lm1) = c1 * (dble(zflm(lm1))-dble(zflm(lm2)))
               Else
                  rflm (lm1) = c1 * (dble(zflm(lm1))+dble(zflm(lm2)))
               End If
            Else If (m .Lt. 0) Then
               If (Mod(m, 2) .Ne. 0) Then
                  rflm (lm1) = - c1 * &
                 & (aimag(zflm(lm1))+aimag(zflm(lm2)))
               Else
                  rflm (lm1) = c1 * (aimag(zflm(lm2))-aimag(zflm(lm1)))
               End If
            Else
               rflm (lm1) = dble (zflm(lm1))
            End If
         End Do
      End Do
      Return
End Subroutine
!EOC
