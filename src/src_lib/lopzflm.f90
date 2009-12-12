!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: lopzflm
! !INTERFACE:
!
!
Subroutine lopzflm (lmax, zflm, ld, zlflm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax  : maximum angular momentum (in,integer)
!   zflm  : coefficients of input spherical harmonic expansion
!           (in,complex((lmax+1)**2))
!   ld    : leading dimension (in,integer)
!   zlflm : coefficients of output spherical harmonic expansion
!           (out,complex(ld,3))
! !DESCRIPTION:
!   Applies the angular momentum operator ${\bf L}$ to a function expanded in
!   terms of complex spherical harmonics. This makes use of the identities
!   \begin{align*}
!    (L_x+iL_y)Y_{lm}(\theta,\phi)&=\sqrt{(l-m)(l+m+1)}Y_{lm+1}(\theta,\phi)\\
!    (L_x-iL_y)Y_{lm}(\theta,\phi)&=\sqrt{(l+m)(l-m+1)}Y_{lm-1}(\theta,\phi)\\
!    L_zY_{lm}(\theta,\phi)&=mY_{lm}(\theta,\phi).
!   \end{align*}
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Complex (8), Intent (In) :: zflm (*)
      Integer, Intent (In) :: ld
      Complex (8), Intent (Out) :: zlflm (ld, 3)
! local variables
      Integer :: l, m, lm
      Real (8) :: t1
      Complex (8) zt1
      If (lmax .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(lopzflm): lmax < 0 : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      lm = 0
      Do l = 0, lmax
         Do m = - l, l
            lm = lm + 1
            If (m .Eq.-l) Then
               zlflm (lm, 1) = 0.d0
               zlflm (lm, 2) = 0.d0
            End If
            If (m .Lt. l) Then
               t1 = 0.5d0 * Sqrt (dble((l-m)*(l+m+1)))
               zt1 = t1 * zflm (lm)
               zlflm (lm+1, 1) = zt1
               zlflm (lm+1, 2) = cmplx (aimag(zt1),-dble(zt1), 8)
            End If
            If (m .Gt.-l) Then
               t1 = 0.5d0 * Sqrt (dble((l+m)*(l-m+1)))
               zt1 = t1 * zflm (lm)
               zlflm (lm-1, 1) = zlflm (lm-1, 1) + zt1
               zlflm (lm-1, 2) = zlflm (lm-1, 2) + cmplx (-aimag(zt1), &
              & dble(zt1), 8)
            End If
            If (m .Ne. 0) Then
               zlflm (lm, 3) = dble (m) * zflm (lm)
            Else
               zlflm (lm, 3) = 0.d0
            End If
         End Do
      End Do
      Return
End Subroutine
!EOC
