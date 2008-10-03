
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: updatpos
! !INTERFACE:
subroutine updatpos
! !USES:
use modmain
! !DESCRIPTION:
!   Updates the current atomic positions according to the force on each atom. If
!   ${\bf r}_{ij}^m$ is the position and ${\bf F}_{ij}^m$ is the force acting on
!   it for atom $j$ of species $i$ and after time step $m$, then the new
!   position is calculated by
!   $$ {\bf r}_{ij}^{m+1}={\bf r}_{ij}^m+\tau_{ij}^m\left({\bf F}_{ij}^m
!    +{\bf F}_{ij}^{m-1}\right), $$
!   where $\tau_{ij}^m$ is a parameter governing the size of the displacement.
!   If ${\bf F}_{ij}^m\cdot{\bf F}_{ij}^{m-1}>0$ then $\tau_{ij}^m$ is
!   increased, otherwise it is decreased.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,ispn,is,ia,ias
real(8) t1
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the dot-product between the current and previous total force
    t1=dot_product(forcetot(:,ias),forcetp(:,ias))
! if the force is in the same direction then increase step size parameter
    if (t1.gt.0.d0) then
      tauatm(ias)=tauatm(ias)+tau0atm
    else
      tauatm(ias)=tau0atm
    end if
! check for negative mass
    if (spmass(is).gt.0.d0) then
      atposc(:,ia,is)=atposc(:,ia,is)+tauatm(ias)*(forcetot(:,ias) &
       +forcetp(:,ias))
    end if
  end do
end do
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the lattice coordinates of the atomic positions
    call r3mv(ainv,atposc(:,ia,is),atposl(:,ia,is))
! set the previous to the current total force
    forcetp(:,ias)=forcetot(:,ias)
  end do
end do
! write lattice vectors and optimised atomic positions to file
call writegeom(.true.)
! write the optimised interatomic distances to file
call writeiad(.true.)
! check for overlapping muffin-tins
call checkmt
! generate structure factors for G-vectors
call gensfacgp(ngvec,vgc,ngvec,sfacg)
! generate the characteristic function
call gencfun
! generate structure factors for G+k-vectors
do ik=1,nkpt
  do ispn=1,nspnfv
    call gensfacgp(ngk(ispn,ik),vgkc(:,:,ispn,ik),ngkmax,sfacgk(:,:,ispn,ik))
  end do
end do
! determine the new nuclear-nuclear energy
call energynn
return
end subroutine
!EOC
