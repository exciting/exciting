
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potcoul
! !INTERFACE:
subroutine potcoul
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the Coulomb potential of the real charge density stored in the
!   global variables {\tt rhomt} and {\tt rhoir} by solving Poisson's equation.
!   These variables are coverted to complex representations and passed to the
!   routine {\tt zpotcoul}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ir,lmax
complex(8) zrho0
! allocatable arrays
real(8), allocatable :: jlgr(:,:,:)
complex(8), allocatable :: zpchg(:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
allocate(jlgr(0:lmaxvr+npsden+1,ngvec,nspecies))
allocate(zpchg(natmtot))
allocate(zrhomt(lmmaxvr,nrmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrmtmax,natmtot))
allocate(zvclir(ngrtot))
! convert real muffin-tin charge density to complex spherical harmonic expansion
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      call rtozflm(lmaxvr,rhomt(1,ir,ias),zrhomt(1,ir,ias))
    end do
  end do
end do
! store real interstitial charge density in complex array
zrhoir(:)=rhoir(:)
! set up the point charge array
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    zpchg(ias)=spzn(is)
  end do
end do
! compute the required spherical Bessel functions
lmax=lmaxvr+npsden+1
call genjlgpr(lmax,gc,jlgr)
! solve the complex Poisson's equation
call zpotcoul(nrmt,nrmtmax,spnrmax,spr,1,gc,jlgr,ylmg,sfacg,zpchg,zrhomt, &
 zrhoir,zvclmt,zvclir,zrho0)
! convert complex muffin-tin potential to real spherical harmonic expansion
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      call ztorflm(lmaxvr,zvclmt(1,ir,ias),vclmt(1,ir,ias))
    end do
  end do
end do
! store complex interstitial potential in real array
vclir(:)=dble(zvclir(:))
! symmetrise the potential
call symrf(1,vclmt,vclir)
deallocate(jlgr,zpchg,zrhomt,zrhoir,zvclmt,zvclir)
return
end subroutine
!EOC

