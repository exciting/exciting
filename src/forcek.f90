
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: forcek
! !INTERFACE:
subroutine forcek(ik,ff)
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the {\bf k}-dependent contribution to the incomplete basis set
!   (IBS) force. See the calling routine {\tt force} for a full description.
!
! !REVISION HISTORY:
!   Created June 2006 (JKD)
!   Updated for spin-spiral case, May 2007 (Francesco Cricchio and JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: ff(ngvec,nspecies)
! local variables
integer ist,i,j,k,iv(3),ig
integer ispn,jspn,is,ia,ias
real(8) eval,sum,v2(3),t1,t2,t3
complex(8) zt1
! automatic arrays
real(8) v1(3,natmtot)
! allocatable arrays
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: h(:)
complex(8), allocatable :: o(:)
allocate(evalfv(nstfv,nspnfv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(h(nmatmax))
allocate(o(nmatmax))
! get the second-variational eigenvalues/vectors and occupancies from file
call getevalfv(vkl(1,ik),evalfv)
call getoccsv(vkl(1,ik),occsv(1,ik))
call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
call getevecsv(vkl(1,ik),evecsv)
! begin loop over first-variational spin components
do ispn=1,nspnfv
! find the matching coefficients
  call match(ngk(ik,ispn),gkc(1,ik,ispn),tpgkc(1,1,ik,ispn), &
   sfacgk(1,1,ik,ispn),apwalm)
! begin loop over first-variational states
  do ist=1,nstfv
    eval=evalfv(ist,ispn)
    v1(:,:)=0.d0
! APW-APW interstitial contribution
    do i=1,ngk(ik,ispn)
      do j=i+1,ngk(ik,ispn)
        iv(:)=ivg(:,igkig(i,ik,ispn))-ivg(:,igkig(j,ik,ispn))
        iv(:)=modulo(iv(:)-intgv(:,1),ngrid(:))+intgv(:,1)
        ig=ivgig(iv(1),iv(2),iv(3))
        if ((ig.gt.0).and.(ig.le.ngvec)) then
! factor of -2 from i times sum over G and G'
          t1=dot_product(vgkc(:,i,ik,ispn),vgkc(:,j,ik,ispn))-2.d0*eval
          zt1=conjg(evecfv(i,ist,ispn))*evecfv(j,ist,ispn)
          do is=1,nspecies
            do ia=1,natoms(is)
              ias=idxas(ia,is)
              t2=ff(ig,is)*aimag(zt1*conjg(sfacg(ig,ias)))
              t3=t1*t2
              v1(:,ias)=v1(:,ias)+t3*vgc(:,ig)
            end do
          end do
        end if
      end do
    end do
! APW-APW and APW-local-orbital muffin-tin contributions
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        h(1:nmat(ik,ispn))=0.d0
        o(1:nmat(ik,ispn))=0.d0
        call hmlaa(.true.,is,ia,ngk(ik,ispn),apwalm,evecfv(1,ist,ispn),h)
        call olpaa(.true.,is,ia,ngk(ik,ispn),apwalm,evecfv(1,ist,ispn),o)
        call hmlalo(.true.,is,ia,ngk(ik,ispn),apwalm,evecfv(1,ist,ispn),h)
        call olpalo(.true.,is,ia,ngk(ik,ispn),apwalm,evecfv(1,ist,ispn),o)
        do k=1,3
          sum=0.d0
          do i=1,ngk(ik,ispn)
! factor of -2 from i times sum over G and G'
            zt1=cmplx(aimag(evecfv(i,ist,ispn)),dble(evecfv(i,ist,ispn)),8)
            sum=sum+vgkc(k,i,ik,ispn)*dble(zt1*(h(i)-eval*o(i)))
          end do
          v1(k,ias)=v1(k,ias)+2.d0*sum
        end do
! end loops over atoms and species
      end do
    end do
! correction to second-variational eigenvalues
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        v2(:)=0.d0
        do jspn=1,nspinor
! check for spin-spiral case
          if ((spinsprl).and.(ispn.ne.jspn)) goto 10
          i=(jspn-1)*nstfv+ist
          if (tevecsv) then
            do j=1,nstsv
              t1=dble(evecsv(i,j))**2+aimag(evecsv(i,j))**2
              t2=t1*wkpt(ik)*occsv(j,ik)
              v2(:)=v2(:)+t2*v1(:,ias)
            end do
          else
            t2=wkpt(ik)*occsv(i,ik)
            v2(:)=v2(:)+t2*v1(:,ias)
          end if
10 continue
        end do
!$OMP CRITICAL
        forceibs(:,ias)=forceibs(:,ias)+v2(:)
!$OMP END CRITICAL
      end do
    end do
! end loop over first-variational states
  end do
! end loop over first-variational spin components
end do
deallocate(evalfv,apwalm,evecfv,evecsv,h,o)
return
end subroutine
!EOC

