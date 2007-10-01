
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: olprad
! !INTERFACE:
subroutine olprad
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the radial overlap integrals of the APW and local-orbital basis
!   functions. In other words, for spin $\sigma$ and atom $j$ of species $i$, it
!   computes integrals of the form
!   $$ o^{\sigma;ij}_{qp}=\int_0^{R_i}u^{\sigma;ij}_{q;l_p}(r)v^{\sigma;ij}_p(r)
!    r^2dr $$
!   and
!   $$ o^{\sigma;ij}_{pp'}=\int_0^{R_i}v^{\sigma;ij}_p(r)v^{\sigma;ij}_{p'}(r)
!    r^2dr,\quad l_p=l_{p'} $$
!   where $u^{\sigma;ij}_{q;l}$ is the $q$th APW radial function for angular
!   momentum $l$; and $v^{\sigma;ij}_p$ is the $p$th local-orbital radial
!   function and has angular momentum $l_p$.
!
! !REVISION HISTORY:
!   Created November 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ir,nr
integer l,ilo,ilo1,ilo2,io
! automatic arrays
real(8) r2(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
do is=1,nspecies
  nr=nrmt(is)
  do ir=1,nr
    r2(ir)=spr(ir,is)**2
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!--------------------------------------!
!     APW-local-orbital integtrals     !
!--------------------------------------!
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      do io=1,apword(l,is)
        do ir=1,nr
          fr(ir)=apwfr(ir,1,io,l,ias)*lofr(ir,1,ilo,ias)*r2(ir)
        end do
        call fderiv(-1,nr,spr(1,is),fr,gr,cf)
        oalo(io,ilo,ias)=gr(nr)
      end do
    end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
    do ilo1=1,nlorb(is)
      l=lorbl(ilo1,is)
      do ilo2=1,nlorb(is)
        if (lorbl(ilo2,is).eq.l) then
          do ir=1,nr
            fr(ir)=lofr(ir,1,ilo1,ias)*lofr(ir,1,ilo2,ias)*r2(ir)
          end do
          call fderiv(-1,nr,spr(1,is),fr,gr,cf)
          ololo(ilo1,ilo2,ias)=gr(nr)
        end if
      end do
    end do
  end do
end do
return
end subroutine
!EOC
