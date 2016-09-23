!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: genrlm
! !INTERFACE:
!
!
Subroutine genrlm (lmax, tp, rlm)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum angular momentum (in,integer)
!   tp   : (theta, phi) coordinates (in,real(2))
!   rlm  : array of real spherical harmonics (out,real((lmax+1)**2))
! !DESCRIPTION:
!   Generates a sequence of real spherical harmonics evaluated at angles
!   $(\theta,\phi)$ for $0<l<l_{\rm max}$. The values are returned in a packed
!   array {\tt rlm} indexed with $j=l(l+1)+m+1$. Real spherical harmonics are
!   defined by
!   $$ R_{lm}(\theta,\phi)= \begin{cases}
!     \sqrt{2}\,\Re\{Y_{lm}(\theta,\phi)\} & m>0 \\
!     \sqrt{2}\,\Im\{Y_{lm}(\theta,\phi)\} & m<0 \\
!     \Re\{Y_{lm}(\theta,\phi)\} & m=0
!    \end{cases}, $$
!   where $Y_{lm}$ are the complex spherical harmonics. These functions are
!   orthonormal and complete and may be used for expanding real-valued functions
!   on the sphere. This routine is numerically stable and accurate to near
!   machine precision for $l\le 50$. See routine {\tt genylm}.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Real (8), Intent (In) :: tp (2)
      Real (8), Intent (Out) :: rlm (*)
! local variables
      Integer :: lmmax, l, m, lm
      Real (8), Parameter :: sqtwo = 1.4142135623730950488d0
! allocatable arrays
      Complex (8), Allocatable :: ylm (:)
      If ((lmax .Lt. 0) .Or. (lmax .Gt. 50)) Then
         Write (*,*)
         Write (*, '("Error(genrlm): lmax out of range : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      lmmax = (lmax+1) ** 2
      Allocate (ylm(lmmax))
! generate complex spherical harmonics
      Call genylm (lmax, tp, ylm)
! convert to real spherical harmonics
      lm = 0
      Do l = 0, lmax
         Do m = - l, - 1
            lm = lm + 1
            rlm (lm) = sqtwo * aimag (ylm(lm))
         End Do
         lm = lm + 1
         rlm (lm) = dble (ylm(lm))
         Do m = 1, l
            lm = lm + 1
            rlm (lm) = sqtwo * dble (ylm(lm))
         End Do
      End Do
      Deallocate (ylm)
      Return
End Subroutine
!EOC
