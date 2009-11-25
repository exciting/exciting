!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: sbessel
! !INTERFACE:
!
!
Subroutine sbessel (lmax, x, jl)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum order of Bessel function (in,integer)
!   x    : real argument (in,real)
!   jl   : array of returned values (out,real(0:lmax))
! !DESCRIPTION:
!   Computes the spherical Bessel functions of the first kind, $j_l(x)$, for
!   argument $x$ and $l=0,1,\ldots,l_{\rm max}$. The recursion relation
!   $$ j_{l+1}(x)=\frac{2l+1}{x}j_l(x)-j_{l-1}(x) $$
!   is used either downwards for $x<l$ or upwards for $x\ge l$. For $x\ll 1$
!   the asymtotic form is used
!   $$ j_l(x)\approx\frac{x^l}{(2l+1)!!}. $$
!   This procedure is numerically stable and accurate to near machine precision
!   for $l\le 50$.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!   Modified to return an array of values, October 2004 (JKD)
!   Improved stability, August 2006 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Real (8), Intent (In) :: x
      Real (8), Intent (Out) :: jl (0:lmax)
! local variables
! staring value for l above lmax (suitable for lmax < 50)
      Integer, Parameter :: lst = 25
      Integer :: l
! rescale limit
      Real (8), Parameter :: rsc = 1.d150
      Real (8), Parameter :: rsci = 1.d0 / rsc
      Real (8) :: xi, j0, j1, jt, t1, t2
      If ((lmax .Lt. 0) .Or. (lmax .Gt. 50)) Then
         Write (*,*)
         Write (*, '("Error(sbessel): lmax out of range : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      If ((x .Lt. 0.d0) .Or. (x .Gt. 1.d5)) Then
         Write (*,*)
         Write (*, '("Error(sbessel): x out of range : ", G18.10)') x
         Write (*,*)
         Stop
      End If
! treat x << 1
      If (x .Lt. 1.d-8) Then
         jl (0) = 1.d0
         t1 = 1.d0
         t2 = 1.d0
         Do l = 1, lmax
            t1 = t1 / dble (2*l+1)
            t2 = t2 * x
            jl (l) = t2 * t1
         End Do
         Return
      End If
      xi = 1.d0 / x
! for x < lmax recurse down
      If (x .Lt. dble(lmax)) Then
         If (lmax .Eq. 0) Then
            jl (0) = Sin (x) / x
            Return
         End If
! start from truly random numbers
         j0 = 0.6370354636449841609d0 * rsci
         j1 = 0.3532702964695481204d0 * rsci
         Do l = lmax + lst, lmax + 1, - 1
            jt = j0 * dble (2*l+1) * xi - j1
            j1 = j0
            j0 = jt
! check for overflow
            If (Abs(j0) .Gt. rsc) Then
! rescale
               jt = jt * rsci
               j1 = j1 * rsci
               j0 = j0 * rsci
            End If
         End Do
         Do l = lmax, 0, - 1
            jt = j0 * dble (2*l+1) * xi - j1
            j1 = j0
            j0 = jt
! check for overflow
            If (Abs(j0) .Gt. rsc) Then
! rescale
               jt = jt * rsci
               j1 = j1 * rsci
               j0 = j0 * rsci
               jl (l+1:lmax) = jl (l+1:lmax) * rsci
            End If
            jl (l) = j1
         End Do
! rescaling constant
         t1 = 1.d0 / ((jl(0)-x*jl(1))*Cos(x)+x*jl(0)*Sin(x))
         jl (:) = t1 * jl (:)
         Return
      Else
! for large x recurse up
         jl (0) = Sin (x) * xi
         If (lmax .Eq. 0) Return
         jl (1) = (jl(0)-Cos(x)) * xi
         If (lmax .Eq. 1) Return
         j0 = jl (0)
         j1 = jl (1)
         Do l = 2, lmax
            jt = dble (2*l-1) * j1 * xi - j0
            j0 = j1
            j1 = jt
            jl (l) = j1
         End Do
         Return
      End If
End Subroutine
!EOC
