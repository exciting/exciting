!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: zpotclmt
! !INTERFACE:
!
!
Subroutine zpotclmt (ptnucl, lmax, nr, r, zn, ld, zrhomt, zvclmt)
! !INPUT/OUTPUT PARAMETERS:
!   ptnucl : .true. if the nucleus is a point particle (in,logical)
!   lmax   : maximum angular momentum (in,integer)
!   nr     : number of radial mesh points (in,integer)
!   r      : radial mesh (in,real(nr))
!   zn     : nuclear charge at the atomic center (in,real)
!   ld     : leading dimension (in,integer)
!   zrhomt : muffin-tin charge density (in,complex(ld,nr))
!   zvclmt : muffin-tin Coulomb potential (out,complex(ld,nr))
! !DESCRIPTION:
!   Solves the Poisson equation for the charge density contained in an isolated
!   muffin-tin using the Green's function approach. In other words, the
!   spherical harmonic expansion of the Coulomb potential, $V_{lm}$, is obtained
!   from the density expansion, $\rho_{lm}$, by
!   $$ V_{lm}(r)=\frac{4\pi}{2l+1}\left(\frac{1}{r^{l+1}}\int_0^r
!    \rho_{lm}(r'){r'}^{l+2}dr'+r^l\int_r^R\frac{\rho_{lm}(r')}{{r'}^{l-1}}dr'
!    \right)+\frac{1}{Y_{00}}\frac{z}{r}\delta_{l,0} $$
!   where the last term is the monopole arising from the point charge $z$, and
!   $R$ is the muffin-tin radius.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: ptnucl
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: zn
      Integer, Intent (In) :: ld
      Complex (8), Intent (In) :: zrhomt (ld, nr)
      Complex (8), Intent (Out) :: zvclmt (ld, nr)
! local variables
      Integer :: l, m, lm, ir
      Real (8), Parameter :: fourpi = 12.566370614359172954d0
! spherical harmonic for l=m=0
      Real (8), Parameter :: y00 = 0.28209479177387814347d0
      Real (8) :: t1, t2, t3, t4, t5, t6, t7
! automatic arrays
      Real (8) :: ri (nr), rl (nr), ril1 (nr), cf (3, nr), vn (nr)
      Real (8) :: fr1 (nr), fr2 (nr), fr3 (nr), fr4 (nr)
      Real (8) :: gr1 (nr), gr2 (nr), gr3 (nr), gr4 (nr)
! initialise r^l and r^(-l-1)
      Do ir = 1, nr
         rl (ir) = 1.d0
         ri (ir) = 1.d0 / r (ir)
         ril1 (ir) = ri (ir)
      End Do
      lm = 0
      Do l = 0, lmax
         t1 = fourpi / dble (2*l+1)
         Do m = - l, l
            lm = lm + 1
            Do ir = 1, nr
               t2 = rl (ir) * r (ir) ** 2
               t3 = ril1 (ir) * r (ir) ** 2
               t4 = dble (zrhomt(lm, ir))
               t5 = aimag (zrhomt(lm, ir))
               fr1 (ir) = t2 * t4
               fr2 (ir) = t2 * t5
               fr3 (ir) = t3 * t4
               fr4 (ir) = t3 * t5
            End Do
            Call fderiv (-1, nr, r, fr1, gr1, cf)
            Call fderiv (-1, nr, r, fr2, gr2, cf)
            Call fderiv (-1, nr, r, fr3, gr3, cf)
            Call fderiv (-1, nr, r, fr4, gr4, cf)
            t2 = gr3 (nr)
            t3 = gr4 (nr)
            Do ir = 1, nr
               t4 = ril1 (ir)
               t5 = rl (ir)
               t6 = t4 * gr1 (ir) + t5 * (t2-gr3(ir))
               t7 = t4 * gr2 (ir) + t5 * (t3-gr4(ir))
               zvclmt (lm, ir) = t1 * cmplx (t6, t7, 8)
            End Do
         End Do
! update r^l and r^(-l-1)
         If (l .Lt. lmax) Then
            Do ir = 1, nr
               rl (ir) = rl (ir) * r (ir)
               ril1 (ir) = ril1 (ir) * ri (ir)
            End Do
         End If
      End Do

! add the nuclear potential
if (.false.) then
      If (zn .Ne. 0.d0) Then
         Call potnucl (ptnucl, nr, r, zn, vn)
         t1 = 1.d0 / y00
         Do ir = 1, nr
            zvclmt (1, ir) = zvclmt (1, ir) + t1 * vn (ir)
         End Do
      End If
endif
      Return
End Subroutine
!EOC
