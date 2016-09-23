!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: zflmconj
! !INTERFACE:
!
!
Subroutine zflmconj (lmax, zflm1, zflm2)
! !INPUT/OUTPUT PARAMETERS:
!   lmax  : maximum angular momentum (in,integer)
!   zflm1 : coefficients of input complex spherical harmonic expansion
!           (in,complex((lmax+1)**2)))
!   zflm2 : coefficients of output complex spherical harmonic expansion
!           (out,complex((lmax+1)**2)))
! !DESCRIPTION:
!   Returns the complex conjugate of a function expanded in spherical harmonics.
!   In other words, given the input function coefficients $z_{lm}$, the routine
!   returns  $z'_{lm}=(-1)^m z^*_{l-m}$ so that
!   $$ \sum_{lm}z'_{lm}Y_{lm}(\theta,\phi)=\left(\sum_{lm}z_{lm}Y_{lm}
!    (\theta,\phi)\right)^* $$
!   for all $(\theta,\phi)$. Note that {\tt zflm1} and {\tt zflm2} can refer to
!   the same array.
!
! !REVISION HISTORY:
!   Created April 2004 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Complex (8), Intent (In) :: zflm1 (*)
      Complex (8), Intent (Out) :: zflm2 (*)
! local variables
      Integer :: l, m, lm1, lm2
      Complex (8) zt1
      Do l = 0, lmax
         lm1 = l ** 2
         lm2 = (l+1) ** 2 + 1
         Do m = - l, - 1
            lm1 = lm1 + 1
            lm2 = lm2 - 1
            zt1 = zflm1 (lm1)
            zflm2 (lm1) = conjg (zflm1(lm2))
            zflm2 (lm2) = conjg (zt1)
            If (Mod(m, 2) .Ne. 0) Then
               zflm2 (lm1) = - zflm2 (lm1)
               zflm2 (lm2) = - zflm2 (lm2)
            End If
         End Do
! m=0 case
         lm1 = lm1 + 1
         zflm2 (lm1) = conjg (zflm1(lm1))
      End Do
      Return
End Subroutine
!EOC
