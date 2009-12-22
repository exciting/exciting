!
!
!
! This routine is based on code written by K. Burke.
!
!
Subroutine x_pbe (kappa, mu, rho, s, u, v, ex, vx)
      Implicit None
! arguments
      Real (8), Intent (In) :: kappa
      Real (8), Intent (In) :: mu
      Real (8), Intent (In) :: rho
      Real (8), Intent (In) :: s
      Real (8), Intent (In) :: u
      Real (8), Intent (In) :: v
      Real (8), Intent (Out) :: ex
      Real (8), Intent (Out) :: vx
! local variables
      Real (8), Parameter :: pi = 3.1415926535897932385d0
      Real (8), Parameter :: ax = - 0.7385587663820224058d0
      Real (8), Parameter :: thrd = 1.d0 / 3.d0
      Real (8), Parameter :: thrd4 = 4.d0 / 3.d0
      Real (8) :: ul, exu, s2, p0
      Real (8) :: fxpbe, fs, fss
      ul = mu / kappa
! LDA exchange energy density
      exu = ax * rho ** thrd
! PBE enhancement factor
      s2 = s ** 2
      p0 = 1.d0 + ul * s2
      fxpbe = 1.d0 + kappa - kappa / p0
      ex = exu * fxpbe
      fs = 2.d0 * kappa * ul / (p0*p0)
      fss = - 4.d0 * ul * s * fs / p0
! exchange potential
      vx = exu * (thrd4*fxpbe-(u-thrd4*s2*s)*fss-v*fs)
      Return
End Subroutine
