
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: energy
! !INTERFACE:
subroutine energy
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the total energy and its individual contributions. The kinetic
!   energy is given by
!   $$ T_{\rm s}=\sum_i n_i\epsilon_i-\int\rho({\bf r})[v_{\rm C}({\bf r})
!    +v_{\rm xc}({\bf r})]d{\bf r}-\int {\bf m}({\bf r})\cdot
!    ({\bf B}_{\rm xc}({\bf r})+{\bf B}_{\rm ext}({\bf r}))d{\bf r}, $$
!   where $n_i$ are the occupancies and $\epsilon_i$ are the eigenvalues of both
!   the core and valence states; $\rho$ is the density; ${\bf m}$ is the
!   magnetisation density; $v_{\rm C}$ is the Coulomb potential; $v_{\rm xc}$
!   and ${\bf B}_{\rm xc}$ are the exchange-correlation potential and effective
!   magnetic field, respectively; and ${\bf B}_{\rm ext}$ is the external
!   magnetic field. The Hartree, electron-nuclear and nuclear-nuclear
!   electrostatic energies are combined into the Coulomb energy:
!   \begin{align*}
!    E_{\rm C}&=E_{\rm H}+E_{\rm en}+E_{\rm nn} \\
!             &=\frac{1}{2}V_{\rm C}+E_{\rm Mad},
!   \end{align*}
!   where
!   $$ V_{\rm C}=\int\rho({\bf r})v_{\rm C}({\bf r})d{\bf r} $$
!   is the Coulomb potential energy. The Madelung energy is given by
!   $$ E_{\rm Mad}=\frac{1}{2}\sum_{\alpha}z_{\alpha}R_{\alpha}, $$
!   where
!   $$ R_{\alpha}=\lim_{r\rightarrow 0}\left(v^{\rm C}_{\alpha;00}(r)Y_{00}
!    +\frac{z_{\alpha}}{r}\right) $$
!   for atom $\alpha$, with $v^{\rm C}_{\alpha;00}$ being the $l=0$ component of
!   the spherical harmonic expansion of $v_{\rm C}$ in the muffin-tin, and
!   $z_{\alpha}$ is the nuclear charge. Using the nuclear-nuclear energy
!   determined at the start of the calculation, the electron-nuclear and Hartree
!   energies can be isolated with
!   $$ E_{\rm en}=2\left(E_{\rm Mad}-E_{\rm nn}\right) $$
!   and
!   $$ E_{\rm H}=\frac{1}{2}(E_{\rm C}-E_{\rm en}). $$
!   Finally, the total energy is
!   $$ E=T_{\rm s}+E_{\rm C}+E_{\rm xc}, $$
!   where $E_{\rm xc}$ is obtained either by integrating the
!   exchange-correlation energy density, or in the case of exact exchange, the
!   explicit calculation of the Fock exchange integral. If the global variable
!   {\tt bfinite} is {\tt .false.} then the external magnetic fields are taken
!   to be infinitesimal (i.e. their only purpose is to break the spin symmetry).
!   In this case the external field contribution to the total energy is set to
!   zero. When {\tt bfinite} is {\tt .true.} the integral over the unit cell
!   $$ E_{{\bf B}\rm ext}=\int{\bf m}({\bf r})\cdot{\bf B}_{\rm ext}({\bf r})
!    d{\bf r} $$
!   is added to the total energy. See {\tt potxc}, {\tt exxengy} and related
!   subroutines.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ik,ist,idm,jdm
! fine structure constant
real(8), parameter :: alpha=1.d0/137.03599911d0
! electron g factor
real(8), parameter :: ge=2.0023193043718d0
real(8), parameter :: ga4=ge*alpha/4.d0
! external functions
real(8) rfinp
external rfinp
!-----------------------------------------------!
!     exchange-correlation potential energy     !
!-----------------------------------------------!
engyvxc=rfinp(1,rhomt,vxcmt,rhoir,vxcir)
!-----------------------------------------------------!
!     exchange-correlation effective field energy     !
!-----------------------------------------------------!
engybxc=0.d0
do idm=1,ndmag
  engybxc=engybxc+rfinp(1,magmt(1,1,1,idm),bxcmt(1,1,1,idm),magir(1,idm), &
   bxcir(1,idm))
end do
!----------------------------------------!
!     external magnetic field energy     !
!----------------------------------------!
engybext=0.d0
do idm=1,ndmag
  if (ndmag.eq.3) then
    jdm=idm
  else
    jdm=3
  end if
  engybext=engybext+ga4*momir(idm)*bfieldc(jdm)
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      engybext=engybext+ga4*mommt(idm,ias)*bfcmt(jdm,ia,is)
    end do
  end do
end do
!----------------------------------!
!     Coulomb potential energy     !
!----------------------------------!
engyvcl=rfinp(1,rhomt,vclmt,rhoir,vclir)
!-----------------------!
!     Madelung term     !
!-----------------------!
engymad=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    engymad=engymad+0.5d0*spzn(is)*(vclmt(1,1,ias)*y00-spzn(is)/spr(1,is))
  end do
end do
!---------------------------------------------!
!     electron-nuclear interaction energy     !
!---------------------------------------------!
engyen=2.d0*(engymad-engynn)
!------------------------!
!     Hartree energy     !
!------------------------!
engyhar=0.5d0*(engyvcl-engyen)
!------------------------!
!     Coulomb energy     !
!------------------------!
engycl=engynn+engyen+engyhar
!-------------------------!
!     exchange energy     !
!-------------------------!
if ((xctype.lt.0).or.(hartfock)) then
! exact exchange calculation using Kohn-Sham orbitals
  if (tlast) call exxengy
else
! exchange energy from the density
  engyx=rfinp(1,rhomt,exmt,rhoir,exir)
end if
!----------------------------!
!     correlation energy     !
!----------------------------!
engyc=rfinp(1,rhomt,ecmt,rhoir,ecir)
!----------------------------!
!     sum of eigenvalues     !
!----------------------------!
! core eigenvalues
evalsum=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist=1,spnst(is)
      if (spcore(ist,is)) evalsum=evalsum+spocc(ist,is)*evalcr(ist,ias)
    end do
  end do
end do
! valence eigenvalues
do ik=1,nkpt
  do ist=1,nstsv
    evalsum=evalsum+wkpt(ik)*occsv(ist,ik)*evalsv(ist,ik)
  end do
end do
!------------------------!
!     kinetic energy     !
!------------------------!
engyekn=evalsum-engyvcl-engyvxc-engybxc-engybext
!----------------------!
!     total energy     !
!----------------------!
engytot=engyekn+0.5d0*engyvcl+engymad+engyx+engyc
! add external field energy if required
if (bfinite) then
  engytot=engytot+engybext
end if
return
end subroutine
!EOC

