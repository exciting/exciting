!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rtozflm
! !INTERFACE:
!
!
Subroutine rtozflm (lmax, rflm, zflm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   rflm : coefficients of real spherical harmonic expansion
!          (in,real((lmax+1)**2)))
!   zflm : coefficients of complex spherical harmonic expansion
!          (out,complex((lmax+1)**2)))
! !DESCRIPTION:
!   Converts a real function, $r_{lm}$, expanded in terms of real spherical
!   harmonics into a complex spherical harmonic expansion, $z_{lm}$:
!   $$ z_{lm}=\begin{cases} \frac{1}{\sqrt{2}}(r_{lm}+i(-1)^mr_{l-m}) & m>0 \\
!    \frac{1}{\sqrt{2}}((-1)^mr_{l-m}-ir_{lm}) & m<0 \\
!    r_{lm} & m=0 \end{cases}\;. $$
!   See routine {\tt genrlm}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Real (8), Intent (In) :: rflm (*)
      Complex (8), Intent (Out) :: zflm (*)
! local variables
      Integer :: l, m, lm1, lm2
! real constant 1/sqrt(2)
      Real (8), Parameter :: c1 = 0.7071067811865475244d0
      If (lmax .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(rtozflm): lmax < 0 : ", I8)') lmax
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
                  zflm (lm1) = c1 * cmplx (rflm(lm1),-rflm(lm2), 8)
               Else
                  zflm (lm1) = c1 * cmplx (rflm(lm1), rflm(lm2), 8)
               End If
            Else If (m .Lt. 0) Then
               If (Mod(m, 2) .Ne. 0) Then
                  zflm (lm1) = c1 * cmplx (-rflm(lm2),-rflm(lm1), 8)
               Else
                  zflm (lm1) = c1 * cmplx (rflm(lm2),-rflm(lm1), 8)
               End If
            Else
               zflm (lm1) = rflm (lm1)
            End If
         End Do
      End Do
      Return
End Subroutine
!EOC
