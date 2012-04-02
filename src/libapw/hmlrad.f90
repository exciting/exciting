
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlrad
! !INTERFACE:
subroutine hmlrad
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the radial Hamiltonian integrals of the APW and local-orbital
!   basis functions. In other words, for spin $\sigma$ and atom $j$ of species
!   $i$, it computes integrals of the form
!   $$ h^{\sigma;ij}_{qq';ll'l''m''}=\begin{cases}
!    \int_0^{R_i}u^{\sigma;ij}_{q;l}(r)Hu^{\sigma;ij}_{q';l'}(r)r^2dr & l''=0 \\
!    \int_0^{R_i}u^{\sigma;ij}_{q;l}(r)V^{\sigma;ij}_{l''m''}(r)
!    u^{\sigma;ij}_{q';l'}(r)r^2dr & l''>0 \end{cases}, $$
!   where $u^{\sigma;ij}_{q;l}$ is the $q$th APW radial function for angular
!   momentum $l$; $H$ is the Hamiltonian of the radial Schr\"{o}dinger equation;
!   and $V^{\sigma;ij}_{l''m''}$ is the effective muffin-tin potential. Similar
!   integrals are calculated for APW-local-orbital and
!   local-orbital-local-orbital contributions.
!
! !REVISION HISTORY:
!   Created December 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,nr,ir
integer l1,l2,l3,m2,m3,lm3
integer ilo,ilo1,ilo2,io,io1,io2
real(8) t1
! automatic arrays
real(8) r2(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
! begin loops over atoms and species
do is=1,nspecies
  nr=nrmt(is)
  do ir=1,nr
    r2(ir)=spr(ir,is)**2
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
    do l1=0,lmaxapw
      do io1=1,apword(l1,is)
        do l2=0,lmaxapw
          do io2=1,apword(l2,is)
            if (l1.eq.l2) then
              do ir=1,nr
                fr(ir)=apwfr(ir,1,io1,l1,ias)*apwfr(ir,2,io2,l2,ias)*r2(ir)
              end do
              call fderiv(-1,nr,spr(:,is),fr,gr,cf)
              haa(1,io1,l1,io2,l2,ias)=gr(nr)/y00
            else
              haa(1,io1,l1,io2,l2,ias)=0.d0
            end if
            do l3=1,lmaxvr
              do m3=-l3,l3
                lm3=idxlm(l3,m3)
                do ir=1,nr
                  t1=apwfr(ir,1,io1,l1,ias)*apwfr(ir,1,io2,l2,ias)*r2(ir)
                  fr(ir)=t1*veffmt(lm3,ir,ias)
                end do
                call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                haa(lm3,io1,l1,io2,l2,ias)=gr(nr)
              end do
            end do
          end do
        end do
      end do
    end do
!-------------------------------------!
!     local-orbital-APW integrals     !
!-------------------------------------!
    do ilo=1,nlorb(is)
      l1=lorbl(ilo,is)
      do l2=0,lmaxapw
        do io2=1,apword(l2,is)
          if (l1.eq.l2) then
            do ir=1,nr
              fr(ir)=lofr(ir,1,ilo,ias)*apwfr(ir,2,io2,l2,ias)*r2(ir)
            end do
            call fderiv(-1,nr,spr(:,is),fr,gr,cf)
            hloa(1,ilo,io2,l2,ias)=gr(nr)/y00
          else
            hloa(1,ilo,io2,l2,ias)=0.d0
          end if
          do l3=1,lmaxvr
            do m3=-l3,l3
              lm3=idxlm(l3,m3)
              do ir=1,nr
                t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io2,l2,ias)*r2(ir)
                fr(ir)=t1*veffmt(lm3,ir,ias)
              end do
              call fderiv(-1,nr,spr(:,is),fr,gr,cf)
              hloa(lm3,ilo,io2,l2,ias)=gr(nr)
            end do
          end do
        end do
      end do
    end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
    do ilo1=1,nlorb(is)
      l1=lorbl(ilo1,is)
      do ilo2=1,nlorb(is)
        l2=lorbl(ilo2,is)
        if (l1.eq.l2) then
          do ir=1,nr
            fr(ir)=lofr(ir,1,ilo1,ias)*lofr(ir,2,ilo2,ias)*r2(ir)
          end do
          call fderiv(-1,nr,spr(:,is),fr,gr,cf)
          hlolo(1,ilo1,ilo2,ias)=gr(nr)/y00
        else
          hlolo(1,ilo1,ilo2,ias)=0.d0
        end if
        do l3=1,lmaxvr
          do m3=-l3,l3
            lm3=idxlm(l3,m3)
            do ir=1,nr
              t1=lofr(ir,1,ilo1,ias)*lofr(ir,1,ilo2,ias)*r2(ir)
              fr(ir)=t1*veffmt(lm3,ir,ias)
            end do
            call fderiv(-1,nr,spr(:,is),fr,gr,cf)
            hlolo(lm3,ilo1,ilo2,ias)=gr(nr)
          end do
        end do
      end do
    end do
! end loops over atoms and species
  end do
end do
return
end subroutine
!EOC

