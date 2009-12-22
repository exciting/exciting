!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: sbesseldm
! !INTERFACE:
!
!
Subroutine sbesseldm (m, lmax, x, djl)
! !INPUT/OUTPUT PARAMETERS:
!   m    : order of derivatve (in,integer)
!   lmax : maximum order of Bessel function (in,integer)
!   x    : real argument (in,real)
!   djl  : array of returned values (out,real(0:lmax))
! !DESCRIPTION:
!   Computes the $m$th derivative of the spherical Bessel function of the first
!   kind, $j_l(x)$, for argument $x$ and $l=0,1,\ldots,l_{\rm max}$. For
!   $x\ge 1$ this is done by repeatedly using the relations
!   \begin{align*}
!    \frac{d}{dx}j_l(x)&=\frac{l}{x}j_l(x)-j_{l+1}(x) \\
!    j_{l+1}(x)&=\frac{2l+1}{x}j_l(x)-j_{l-1}(x).
!   \end{align*}
!   While for $x<1$ the series expansion of the Bessel function is used
!   $$ \frac{d^m}{dx^m}j_l(x)=\sum_{i=0}^{\infty}
!    \frac{(2i+l)!}{(-2)^ii!(2i+l-m)!(2i+2l+1)!!}x^{2i+l-m}. $$
!   This procedure is numerically stable and accurate to near machine precision
!   for $l\le 30$ and $m\le 6$.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!   Modified to return an array of values, October 2004 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: m
      Integer, Intent (In) :: lmax
      Real (8), Intent (In) :: x
      Real (8), Intent (Out) :: djl (0:lmax)
! local variables
      Integer, Parameter :: maxm = 6
      Integer, Parameter :: maxns = 20
      Integer :: i, j, l, i0
      Real (8) :: t1, sum, x2
      Integer :: a (0:maxm+1), a1 (0:maxm+1)
      Integer :: b (0:maxm+1), b1 (0:maxm+1)
! automatic arrays
      Real (8) :: jl (0:lmax+1)
! external functions
      Real (8) :: factnm, factr
      External factnm, factr
      If ((m .Lt. 0) .Or. (m .Gt. maxm)) Then
         Write (*,*)
         Write (*, '("Error(sbesseldm): m out of range : ", I8)') m
         Write (*,*)
         Stop
      End If
      If ((lmax .Lt. 0) .Or. (lmax .Gt. 30)) Then
         Write (*,*)
         Write (*, '("Error(sbesseldm): lmax out of range : ", I8)') &
        & lmax
         Write (*,*)
         Stop
      End If
      If ((x .Lt. 0.d0) .Or. (x .Gt. 1.d5)) Then
         Write (*,*)
         Write (*, '("Error(sbesseldm): x out of range : ", G18.10)') x
         Write (*,*)
         Stop
      End If
      If (m .Eq. 0) Then
         Call sbessel (lmax, x, djl)
         Return
      End If
      If (x .Gt. 1.d0) Then
         Call sbessel (lmax+1, x, jl)
         Do l = 0, lmax
            a (1:m+1) = 0
            a (0) = 1
            a1 (0:m+1) = 0
            Do i = 1, m
               b (0) = 0
               b1 (0) = 0
               Do j = 0, i
                  b (j+1) = a (j) * (l-j)
                  b1 (j+1) = - a1 (j) * (j+l+2)
               End Do
               Do j = 0, i
                  b1 (j) = b1 (j) - a (j)
                  b (j) = b (j) + a1 (j)
               End Do
               a (0:i+1) = b (0:i+1)
               a1 (0:i+1) = b1 (0:i+1)
            End Do
            t1 = 1.d0
            sum = dble (a(0)) * jl (l) + dble (a1(0)) * jl (l+1)
            Do i = 1, m + 1
               t1 = t1 * x
               sum = sum + (dble(a(i))*jl(l)+dble(a1(i))*jl(l+1)) / t1
            End Do
            djl (l) = sum
         End Do
      Else
         x2 = x ** 2
         Do l = 0, lmax
            i0 = Max ((m-l+1)/2, 0)
            j = 2 * i0 + l - m
            If (j .Eq. 0) Then
               t1 = 1.d0
            Else
               t1 = x ** j
            End If
            t1 = factr (j+m, j) * t1 / (factnm(i0, 1)*factnm(j+l+m+1, &
           & 2)*dble((-2)**i0))
            sum = t1
            Do i = i0 + 1, maxns
               j = 2 * i + l
               t1 = - t1 * dble ((j-1)*j) * x2 / dble &
              & ((j-l)*(j-m-1)*(j-m)*(j+l+1))
               If (Abs(t1) .Le. 1.d-40) Go To 10
               sum = sum + t1
            End Do
10          Continue
            djl (l) = sum
         End Do
      End If
      Return
End Subroutine
!EOC
