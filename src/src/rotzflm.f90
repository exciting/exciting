!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rotzflm
! !INTERFACE:
!
!
Subroutine rotzflm (rot, lmax, n, ld, zflm1, zflm2)
! !INPUT/OUTPUT PARAMETERS:
!   rot   : rotation matrix (in,real(3,3))
!   lmax  : maximum angular momentum (in,integer)
!   n     : number of functions to rotate (in,integer)
!   ld    : leading dimension (in,integer)
!   zflm1 : coefficients of complex spherical harmonic expansion for each
!           function (in,complex(ld,n))
!   zflm2 : coefficients of rotated functions (out,complex(ld,n))
! !DESCRIPTION:
!   Rotates a set of functions
!   $$ f_i({\bf r})=\sum_{lm}f_{lm}^iY_{lm}(\hat{\bf r}) $$
!   for all $i$, given the coefficients $f_{lm}^i$ and a rotation matrix $R$.
!   This is done by first the computing the Euler angles $(\alpha,\beta,\gamma)$
!   of $R^{-1}$ (see routine {\tt euler}) and then generating the rotation
!   matrix for spherical harmonics, $D^l_{mm'}(\alpha,\beta,\gamma)$, with which
!   $$ Y_{lm}(\theta',\phi')=\sum_{m'}D^l_{mm'}(\alpha,\beta,\gamma)Y_{lm'}
!    (\theta,\phi), $$
!   where $(\theta',\phi')$ are the angles $(\theta,\phi)$ rotated by $R$. The
!   matrix $D$ is given explicitly by
!   \begin{align*}
!    D^l_{mm'}(\alpha,\beta,\gamma)=&\sum_i\frac{(-1)^i\sqrt{(l+m)!(l-m)!(l+m')!
!    (l-m')!}}{(l-m'-i)!(l+m-i)!i!(i+m'-m)!}\\
!    &\times\left(\cos\frac{\beta}{2}\right)^{2l+m-m'-2i}\left(\sin\frac{\beta}
!    {2}\right)^{2i+m'-m}e^{-i(m\alpha+m'\gamma)},
!   \end{align*}
!   where the sum runs over all $i$ which make the factorial arguments
!   non-negative. For improper rotations, i.e. those which are a combination of
!   a rotation and inversion, the rotation is first made proper with
!   $R\rightarrow-R$ and $D$ is modified with
!   $D^l_{mm'}\rightarrow(-1)^l D^l_{mm'}$.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: rot (3, 3)
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: n
      Integer, Intent (In) :: ld
      Complex (8), Intent (In) :: zflm1 (ld, n)
      Complex (8), Intent (Out) :: zflm2 (ld, n)
! local variables
      Integer :: lmmax, l, m1, m2, lm1, lm2
      Integer :: i, j, nm, p
      Real (8) :: det, roti (3, 3), ang (3)
      Real (8) :: cb, sb, sum, t1, t2, t3
      Complex (8), Parameter :: zzero = (0.d0, 0.d0)
      Complex (8), Parameter :: zone = (1.d0, 0.d0)
! allocatable arrays
      Complex (8), Allocatable :: d (:, :)
! external functions
      Real (8) :: factnm
      External factnm
      If (lmax .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(rotzflm): lmax < 0 : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      If (n .Eq. 0) Return
      If (n .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(rotzflm): n < 0 : ", I8)') n
         Write (*,*)
         Stop
      End If
      lmmax = (lmax+1) ** 2
      Allocate (d(lmmax, lmmax))
! find the determinant
      det = rot (1, 2) * rot (2, 3) * rot (3, 1) - rot (1, 3) * rot (2, &
     & 2) * rot (3, 1) + rot (1, 3) * rot (2, 1) * rot (3, 2) - rot (1, &
     & 1) * rot (2, 3) * rot (3, 2) + rot (1, 1) * rot (2, 2) * rot (3, &
     & 3) - rot (1, 2) * rot (2, 1) * rot (3, 3)
! invert rot because the function is to be rotated and not the coordinate system
      Call r3minv (rot, roti)
! make inverse rotation proper
      If (det .Gt. 0.d0) Then
         p = 1
      Else
         p = - 1
         roti (:, :) = - roti (:, :)
      End If
! compute Euler angles of rotation matrix
      Call euler (roti, ang)
      cb = Cos (ang(2)/2.d0)
      sb = Sin (ang(2)/2.d0)
      lm1 = 0
      Do l = 0, lmax
! generate rotation operator for m-components of current l
         Do m1 = - l, l
            lm1 = lm1 + 1
            lm2 = l ** 2
            Do m2 = - l, l
               lm2 = lm2 + 1
               sum = 0.d0
               Do i = 0, Min (l+m1, l-m2)
                  If (((l+m1-i) .Ge. 0) .And. ((l-m2-i) .Ge. 0) .And. &
                 & ((i+m2-m1) .Ge. 0)) Then
                     j = 2 * l + m1 - m2 - 2 * i
                     If (j .Eq. 0) Then
                        t1 = 1.d0
                     Else
                        t1 = cb ** j
                     End If
                     j = 2 * i + m2 - m1
                     If (j .Eq. 0) Then
                        t2 = 1.d0
                     Else
                        t2 = sb ** j
                     End If
                     t3 = t1 * t2 / (factnm(l+m1-i, 1)*factnm(l-m2-i, &
                    & 1)*factnm(i, 1)*factnm(i+m2-m1, 1))
                     If (Mod(i, 2) .Ne. 0) t3 = - t3
                     sum = sum + t3
                  End If
               End Do
               t1 = Sqrt (factnm(l+m1, 1)*factnm(l-m1, 1)*factnm(l+m2, &
              & 1)*factnm(l-m2, 1))
               t2 = - dble (m1) * ang (1) - dble (m2) * ang (3)
               d (lm1, lm2) = sum * t1 * cmplx (Cos(t2), Sin(t2), 8)
               If ((p .Eq.-1) .And. (Mod(l, 2) .Ne. 0)) d (lm1, lm2) = &
              & - d (lm1, lm2)
            End Do
         End Do
! apply rotation operator
         nm = 2 * l + 1
         lm2 = l ** 2 + 1
         Call zgemm ('N', 'N', nm, n, nm, zone, d(lm2, lm2), lmmax, &
        & zflm1(lm2, 1), ld, zzero, zflm2(lm2, 1), ld)
      End Do
      Deallocate (d)
      Return
End Subroutine
!EOC
