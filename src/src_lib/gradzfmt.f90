!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: gradzfmt
! !INTERFACE:
!
!
Subroutine gradzfmt (lmax, nr, r, ld1, ld2, zfmt, gzfmt)
! !INPUT/OUTPUT PARAMETERS:
!   lmax  : maximum angular momentum (in,integer)
!   nr    : number of radial mesh points (in,integer)
!   r     : radial mesh (in,real(nr))
!   ld1   : leading dimension 1 (in,integer)
!   ld2   : leading dimension 2 (in,integer)
!   zfmt  : complex muffin-tin function (in,complex(ld1,nr))
!   gzfmt : gradient of zfmt (out,complex(ld1,ld2,3))
! !DESCRIPTION:
!   Calculates the gradient of a complex muffin-tin function. In other words,
!   given the spherical harmonic expansion coefficients, $f_{lm}(r)$, of a
!   function $f({\bf r})$, the routine returns ${\bf F}_{lm}$ where
!   $$ \sum_{lm}{\bf F}_{lm}(r)Y_{lm}(\hat{\bf r})=\nabla f({\bf r}). $$
!   This is done using the identity
!   \begin{align*}
!    \nabla\left[f_{lm}(r)Y_{lm}(\hat{\bf r})\right]=&-\left[\frac{l+1}{2l+1}
!    \right]^{1/2}\left[\frac{d}{dr}-\frac{l}{r}\right]f_{lm}(r)
!    {\bf Y}_{ll+1m}(\hat{\bf r})\\
!    &+\left[\frac{l}{2l+1}\right]^{1/2}\left[\frac{d}{dr}+\frac{l+1}{r}\right]
!    f_{lm}(r){\bf Y}_{ll-1m}(\hat{\bf r}),
!   \end{align*}
!   where the vector spherical harmonics are given by
!   $$ {\bf Y}_{ll'm}(\hat{\bf r})=\sum_{m'=-l'}^{l'}\sum_{m''=-1}^1
!    C(l'1l|m'm''m)\,Y_{lm}(\hat{\bf r})\,\hat{\bf e}_{m''}, $$
!   $C$ is a Clebsch-Gordan coefficient and
!   $$ \hat{\bf e}_{+1}=-\frac{\hat{\bf x}+i\hat{\bf y}}{\sqrt{2}},\quad\quad
!    \hat{\bf e}_0=\hat{\bf z},\quad\quad \hat{\bf e}_{-1}=\frac{\hat{\bf x}-
!    i\hat{\bf y}}{\sqrt{2}} $$
!   are unit vectors. Note that the gradient returned is in terms of the global
!   $(\hat{\bf x},\hat{\bf y},\hat{\bf z})$ coordinate system.
!
! !REVISION HISTORY:
!   Created August 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Integer, Intent (In) :: ld1
      Integer, Intent (In) :: ld2
      Complex (8), Intent (In) :: zfmt (ld1, nr)
      Complex (8), Intent (Out) :: gzfmt (ld1, ld2, 3)
! local variables
      Integer :: ir, l1, l2, m1, m2, lm1, lm2
! square root of two
      Real (8), Parameter :: sqtwo = 1.4142135623730950488d0
      Real (8) :: t1, t2, t3, t4, t5
      Complex (8), Parameter :: zi = (0.d0, 1.d0)
      Complex (8) zt1, zt2
! automatic arrays
      Real (8) :: f (nr), g1 (nr), g2 (nr), cf (3, nr)
! external functions
      Real (8) :: clebgor
      External clebgor
      If (lmax .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(gradzfmt): lmax < 0 : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      Do ir = 1, nr
         gzfmt (:, ir, :) = 0.d0
      End Do
      lm1 = 0
      Do l1 = 0, lmax
         Do m1 = - l1, l1
            lm1 = lm1 + 1
! compute the radial derivatives
            f (1:nr) = dble (zfmt(lm1, 1:nr))
            Call fderiv (1, nr, r, f, g1, cf)
            f (1:nr) = aimag (zfmt(lm1, 1:nr))
            Call fderiv (1, nr, r, f, g2, cf)
            t1 = Sqrt (dble(l1+1)/dble(2*l1+1))
            t2 = Sqrt (dble(l1)/dble(2*l1+1))
            l2 = l1 + 1
            If (l2 .Le. lmax) Then
               lm2 = l2 ** 2
               Do m2 = - l2, l2
                  lm2 = lm2 + 1
                  t3 = clebgor (l2, 1, l1, m2,-1, m1)
                  t4 = clebgor (l2, 1, l1, m2, 0, m1)
                  t5 = clebgor (l2, 1, l1, m2, 1, m1)
                  Do ir = 1, nr
                     zt1 = cmplx (g1(ir), g2(ir), 8)
                     zt2 = t1 * (zfmt(lm1, ir)*dble(l1)/r(ir)-zt1)
                     gzfmt (lm2, ir, 1) = gzfmt (lm2, ir, 1) + &
                    & ((t3-t5)/sqtwo) * zt2
                     gzfmt (lm2, ir, 2) = gzfmt (lm2, ir, 2) + &
                    & ((-t3-t5)/sqtwo) * zi * zt2
                     gzfmt (lm2, ir, 3) = gzfmt (lm2, ir, 3) + t4 * zt2
                  End Do
               End Do
            End If
            l2 = l1 - 1
            If (l2 .Ge. 0) Then
               lm2 = l2 ** 2
               Do m2 = - l2, l2
                  lm2 = lm2 + 1
                  t3 = clebgor (l2, 1, l1, m2,-1, m1)
                  t4 = clebgor (l2, 1, l1, m2, 0, m1)
                  t5 = clebgor (l2, 1, l1, m2, 1, m1)
                  Do ir = 1, nr
                     zt1 = cmplx (g1(ir), g2(ir), 8)
                     zt2 = t2 * (zfmt(lm1, ir)*dble(l1+1)/r(ir)+zt1)
                     gzfmt (lm2, ir, 1) = gzfmt (lm2, ir, 1) + &
                    & ((t3-t5)/sqtwo) * zt2
                     gzfmt (lm2, ir, 2) = gzfmt (lm2, ir, 2) + &
                    & ((-t3-t5)/sqtwo) * zi * zt2
                     gzfmt (lm2, ir, 3) = gzfmt (lm2, ir, 3) + t4 * zt2
                  End Do
               End Do
            End If
         End Do
      End Do
      Return
End Subroutine
!EOC
