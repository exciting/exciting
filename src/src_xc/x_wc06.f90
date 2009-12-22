!
!
!
! Copyright (C) 2006 Zhigang Wu and R. E. Cohen.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine x_wc06 (rho, s, u, v, ex, vx)
      Implicit None
! arguments
      Real (8), Intent (In) :: rho
      Real (8), Intent (In) :: s
      Real (8), Intent (In) :: u
      Real (8), Intent (In) :: v
      Real (8), Intent (Out) :: ex
      Real (8), Intent (Out) :: vx
! local variables
      Real (8), Parameter :: pi = 3.1415926535897932385d0
      Real (8), Parameter :: ax = - 0.7385587663820224059d0
      Real (8), Parameter :: mu = 0.2195149727645171d0
      Real (8), Parameter :: kappa = 0.804d0
      Real (8), Parameter :: b = 10.d0 / 81.d0
      Real (8), Parameter :: c = 0.00793746933516d0
      Real (8), Parameter :: thrd = 1.d0 / 3.d0
      Real (8), Parameter :: thrd4 = 4.d0 / 3.d0
      Real (8) :: dmu, exu
      Real (8) :: s2, s4, es2, x, p0, fxwc
      Real (8) :: fs, fss, t0, t1, t2, t3
! lda exchange energy density
      exu = ax * rho ** thrd
      s2 = s ** 2
      s4 = s2 ** 2
      es2 = Exp (-s2)
      t0 = 1.d0 + c * s4
      dmu = mu - b
      x = b * s2 + dmu * s2 * es2 + Log (t0)
      p0 = 1.d0 + x / kappa
! WC enhancement factor
      fxwc = 1.d0 + kappa - kappa / p0
! exchange energy density
      ex = exu * fxwc
      t1 = b + dmu * (1.d0-s2) * es2 + 2.d0 * c * s2 / t0
      t2 = dmu * s * (s2-2.d0) * es2 + 2.d0 * c / t0 - 4.d0 * (c**2) * &
     & s4 / (t0**2)
      t3 = 1.d0 / (p0**2)
      fs = 2.d0 * t1 * t3
      fss = t3 * (4.d0*t2-8.d0*s*(t1**2)/(kappa*p0))
! exchange potential
      vx = exu * (thrd4*fxwc-(u-thrd4*s2*s)*fss-v*fs)
      Return
End
