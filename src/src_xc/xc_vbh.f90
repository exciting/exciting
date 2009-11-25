!
!
!
! Copyright (C) 2002-2005 F. Cricchio, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: xc_vbh
! !INTERFACE:
!
!
Subroutine xc_vbh (n, rhoup, rhodn, ex, ec, vxup, vxdn, vcup, vcdn)
! !INPUT/OUTPUT PARAMETERS:
!   n     : number of density points (in,integer)
!   rhoup : spin-up charge density (in,real(n))
!   rhodn : spin-down charge density (in,real(n))
!   ex    : exchange energy density (out,real(n))
!   ec    : correlation energy density (out,real(n))
!   vxup  : spin-up exchange potential (out,real(n))
!   vxdn  : spin-down exchange potential (out,real(n))
!   vcup  : spin-up correlation potential (out,real(n))
!   vcdn  : spin-down correlation potential (out,real(n))
! !DESCRIPTION:
!   Spin-polarised exchange-correlation potential and energy functional of
!   von Barth and Hedin: {\it J. Phys. C} {\bf 5}, 1629 (1972). Note that the
!   implementation is in Rydbergs in order to follow the paper step by step, at
!   the end the potential and energy are converted to Hartree.
!
! !REVISION HISTORY:
!   Created September 2007 (F. Cricchio)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Real (8), Intent (In) :: rhoup (n)
      Real (8), Intent (In) :: rhodn (n)
      Real (8), Intent (Out) :: ex (n)
      Real (8), Intent (Out) :: ec (n)
      Real (8), Intent (Out) :: vxup (n)
      Real (8), Intent (Out) :: vxdn (n)
      Real (8), Intent (Out) :: vcup (n)
      Real (8), Intent (Out) :: vcdn (n)
! local variables
      Integer :: i
      Real (8), Parameter :: pi = 3.1415926535897932385d0
      Real (8), Parameter :: cp = 0.0504d0
      Real (8), Parameter :: cf = 0.0254d0
      Real (8), Parameter :: rp = 30.d0
      Real (8), Parameter :: rf = 75.d0
      Real (8) :: alpha0, a, gamma, r, rs, x, zf, zp
      Real (8) :: fx, fp, ff, eps0_x, epsp_x, mup_x
      Real (8) :: epsp_c, epsf_c, mup_c, muf_c, vc, tau_c
      alpha0 = (4.d0/(9.d0*pi)) ** (1.d0/3.d0)
      eps0_x = (3.d0/2.d0) / (pi*alpha0)
      a = 2.d0 ** (-1.d0/3.d0)
      gamma = (4.d0/3.d0) * a / (1.d0-a)
      Do i = 1, n
         If ((rhoup(i) .Gt. 1.d-12) .And. (rhodn(i) .Gt. 1.d-12)) Then
            r = rhoup (i) + rhodn (i)
! Wigner-Seitz radius in atomic units (a0=1)
            rs = (3.d0/(4.d0*pi*r)) ** (1.d0/3.d0)
            x = rhoup (i) / r
            fx = (1.d0/(1.d0-a)) * &
           & (x**(4.d0/3.d0)+(1.d0-x)**(4.d0/3.d0)-a)
            epsp_x = - eps0_x / rs
            mup_x = (4.d0/3.d0) * epsp_x
! exchange energy
            ex (i) = epsp_x + (1.d0/gamma) * mup_x * fx
            zp = rs / rp
            fp = (1.d0+zp**3) * Log (1.d0+1.d0/zp) + 0.5d0 * zp - zp ** &
           & 2 - 1.d0 / 3.d0
            zf = rs / rf
            ff = (1.d0+zf**3) * Log (1.d0+1.d0/zf) + 0.5d0 * zf - zf ** &
           & 2 - 1.d0 / 3.d0
            epsp_c = - cp * fp
            epsf_c = - cf * ff
            vc = gamma * (epsf_c-epsp_c)
! correlation energy
            ec (i) = epsp_c + (1.d0/gamma) * vc * fx
            mup_c = - cp * Log (1.d0+rp/rs)
            muf_c = - cf * Log (1.d0+rf/rs)
            tau_c = muf_c - mup_c - (4.d0/3.d0) * (epsf_c-epsp_c)
! exchange potential
            vxup (i) = mup_x * (2.d0*x) ** (1.d0/3.d0)
            vxdn (i) = mup_x * (2.d0*(1.d0-x)) ** (1.d0/3.d0)
! correlation potential
            vcup (i) = vc * (2.d0*x) ** (1.d0/3.d0) + mup_c - vc + &
           & tau_c * fx
            vcdn (i) = vc * (2.d0*(1.d0-x)) ** (1.d0/3.d0) + mup_c - vc &
           & + tau_c * fx
! convert from Rybergs to Hartree
            ex (i) = 0.5d0 * ex (i)
            ec (i) = 0.5d0 * ec (i)
            vxup (i) = 0.5d0 * vxup (i)
            vxdn (i) = 0.5d0 * vxdn (i)
            vcup (i) = 0.5d0 * vcup (i)
            vcdn (i) = 0.5d0 * vcdn (i)
         Else
            ex (i) = 0.d0
            ec (i) = 0.d0
            vxup (i) = 0.d0
            vxdn (i) = 0.d0
            vcup (i) = 0.d0
            vcdn (i) = 0.d0
         End If
      End Do
      Return
End Subroutine
!EOC
