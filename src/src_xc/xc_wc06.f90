!
!
!
! Copyright (C) 2006 Zhigang Wu and R. E. Cohen.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine xc_wc06 (n, rho, grho, g2rho, g3rho, ex, ec, vx, vc)
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Real (8), Intent (In) :: rho (n)
      Real (8), Intent (In) :: grho (n)
      Real (8), Intent (In) :: g2rho (n)
      Real (8), Intent (In) :: g3rho (n)
      Real (8), Intent (Out) :: ex (n)
      Real (8), Intent (Out) :: ec (n)
      Real (8), Intent (Out) :: vx (n)
      Real (8), Intent (Out) :: vc (n)
! local variables
      Integer :: i
      Real (8), Parameter :: pi = 3.1415926535897932385d0
      Real (8), Parameter :: thrd = 1.d0 / 3.d0
      Real (8), Parameter :: thrd2 = 2.d0 / 3.d0
! default PBE beta
      Real (8), Parameter :: beta = 0.06672455060314922d0
! maximum allowed |grad rho|
      Real (8), Parameter :: gmax = 1.d6
! maximum allowed grad^2 rho
      Real (8), Parameter :: g2max = 1.d12
! maximum allowed (grad rho).(grad |grad rho|)
      Real (8), Parameter :: g3max = 1.d14
      Real (8) :: r, grho_, g2rho_, g3rho_
      Real (8) :: kf, s, u, v, rs, z, g
      Real (8) :: ks, ksg, t, uu, vv, ww
      Do i = 1, n
         r = rho (i)
         If (r .Gt. 1.d-12) Then
            grho_ = grho (i)
            g2rho_ = g2rho (i)
            g3rho_ = g3rho (i)
! check gradients are within range
            If (grho_ .Gt. gmax) grho_ = gmax
            If (Abs(g2rho_) .Gt. g2max) g2rho_ = sign (g2max, g2rho_)
            If (Abs(g3rho_) .Gt. g3max) g3rho_ = sign (g3max, g3rho_)
            kf = (r*3.d0*pi**2) ** thrd
            s = grho_ / (2.d0*kf*r)
            u = g3rho_ / ((r**2)*(2.d0*kf)**3)
            v = g2rho_ / (r*(2.d0*kf)**2)
! Wu-Cohen exchange
            Call x_wc06 (r, s, u, v, ex(i), vx(i))
! Perdew-Burke-Ernzerhof correlation
            rs = (3.d0/(4.d0*pi*r)) ** thrd
            z = 0.d0
            g = 1.d0
            ks = Sqrt (4.d0*kf/pi)
            ksg = 2.d0 * ks * g
            t = grho_ / (ksg*r)
            uu = g3rho_ / ((r**2)*ksg**3)
            vv = g2rho_ / (r*ksg**2)
            ww = 0.d0
            Call c_pbe (beta, rs, z, t, uu, vv, ww, ec(i), vc(i), &
           & vc(i))
         Else
            ex (i) = 0.d0
            ec (i) = 0.d0
            vx (i) = 0.d0
            vc (i) = 0.d0
         End If
      End Do
      Return
End Subroutine
