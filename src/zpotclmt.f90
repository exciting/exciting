
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zpotclmt
! !INTERFACE:
subroutine zpotclmt(lmax,nr,r,zpchg,ld,zrhomt,zvclmt)
! !INPUT/OUTPUT PARAMETERS:
!   lmax   : maximum angular momentum (in,integer)
!   nr     : number of radial mesh points (in,integer)
!   r      : radial mesh (in,real(nr))
!   zpchg  : point charge  at the atomic center (in,complex(natmtot))
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
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
complex(8), intent(in) :: zpchg
integer, intent(in) :: ld
complex(8), intent(in) :: zrhomt(ld,nr)
complex(8), intent(out) :: zvclmt(ld,nr)
! local variables
integer l,m,lm,ir
real(8), parameter :: fourpi=12.566370614359172954d0
! spherical harmonic for l=m=0
real(8), parameter :: y00=0.28209479177387814347d0
real(8) t1,t2,t3,t4,t5,t6,t7
complex(8) zt1
! automatic arrays
real(8) ri(nr),rl(nr),ril1(nr),cf(3,nr)
real(8) fr1(nr),fr2(nr),fr3(nr),fr4(nr)
real(8) gr1(nr),gr2(nr),gr3(nr),gr4(nr)
! initialise r^l and r^(-l-1)
do ir=1,nr
  rl(ir)=1.d0
  ri(ir)=1.d0/r(ir)
  ril1(ir)=ri(ir)
end do
lm=0
do l=0,lmax
  t1=fourpi/dble(2*l+1)
  do m=-l,l
    lm=lm+1
    do ir=1,nr
      t2=rl(ir)*r(ir)**2
      t3=ril1(ir)*r(ir)**2
      t4=dble(zrhomt(lm,ir))
      t5=aimag(zrhomt(lm,ir))
      fr1(ir)=t2*t4
      fr2(ir)=t2*t5
      fr3(ir)=t3*t4
      fr4(ir)=t3*t5
    end do
    call fderiv(-1,nr,r,fr1,gr1,cf)
    call fderiv(-1,nr,r,fr2,gr2,cf)
    call fderiv(-1,nr,r,fr3,gr3,cf)
    call fderiv(-1,nr,r,fr4,gr4,cf)
    t2=gr3(nr)
    t3=gr4(nr)
    do ir=1,nr
      t4=ril1(ir)
      t5=rl(ir)
      t6=t4*gr1(ir)+t5*(t2-gr3(ir))
      t7=t4*gr2(ir)+t5*(t3-gr4(ir))
      zvclmt(lm,ir)=t1*cmplx(t6,t7,8)
    end do
  end do
! update r^l and r^(-l-1)
  if (l.lt.lmax) then
    do ir=1,nr
      rl(ir)=rl(ir)*r(ir)
      ril1(ir)=ril1(ir)*ri(ir)
    end do
  end if
end do
! add the point charge potential
if (zpchg.ne.(0.d0,0.d0)) then
  zt1=zpchg/y00
  do ir=1,nr
    zvclmt(1,ir)=zvclmt(1,ir)+zt1*ri(ir)
  end do
end if
return
end subroutine
!EOC

