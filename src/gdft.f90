
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gdft(ik,delta)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(out) :: delta(nstsv,nstsv)
! local variables
integer ist,jst,is,ia,ias
integer ir,irc,itp
real(8) sum,t1,t2
! automatic arrays
real(8) fr(nrcmtmax),gr(nrcmtmax),cf(3,nrcmtmax)
real(8) rflm(lmmaxvr),rftp(lmmaxvr)
! allocatable arrays
real(8), allocatable :: rfmt(:,:,:)
real(8), allocatable :: rfir(:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
allocate(rfmt(lmmaxvr,nrcmtmax,natmtot))
allocate(rfir(ngrtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor,nstsv))
! get the eigenvectors from file for input k-point
call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(vkl(:,ik),evecsv)
! find the matching coefficients
call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! calculate the wavefunctions for all second-variational states
call genwfsv(.false.,ngk(1,ik),igkig(:,1,ik),evalsv(:,ik),apwalm,evecfv, &
 evecsv,wfmt,wfir)
do ist=1,nstsv
  delta(ist,ist)=0.d0
  do jst=ist+1,nstsv
! muffin-tin part
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        irc=0
        do ir=1,nrmt(is),lradstp
          irc=irc+1
          rflm(:)=2.d0*(exmt(:,ir,ias)+ecmt(:,ir,ias))-vxcmt(:,ir,ias)
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,rflm,1,0.d0, &
           rftp,1)
          do itp=1,lmmaxvr
            t1=dble(wfmt(itp,irc,ias,1,ist))**2 &
             +aimag(wfmt(itp,irc,ias,1,ist))**2
            t2=dble(wfmt(itp,irc,ias,1,jst))**2 &
             +aimag(wfmt(itp,irc,ias,1,jst))**2
            if (spinpol) then
              t1=t1+dble(wfmt(itp,irc,ias,2,ist))**2 &
               +aimag(wfmt(itp,irc,ias,2,ist))**2
              t2=t2+dble(wfmt(itp,irc,ias,2,jst))**2 &
               +aimag(wfmt(itp,irc,ias,2,jst))**2
            end if
            rftp(itp)=rftp(itp)*(t1-t2)
          end do
          call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp,1,0.d0, &
           rfmt(1,irc,ias),1)
        end do
      end do
    end do
! interstitial part
    do ir=1,ngrtot
      t1=dble(wfir(ir,1,ist))**2+aimag(wfir(ir,1,ist))**2
      t2=dble(wfir(ir,1,jst))**2+aimag(wfir(ir,1,jst))**2
      if (spinpol) then
        t1=t1+dble(wfir(ir,2,ist))**2+aimag(wfir(ir,2,ist))**2
        t2=t2+dble(wfir(ir,2,jst))**2+aimag(wfir(ir,2,jst))**2
      end if
      rfir(ir)=(2.d0*(exir(ir)+ecir(ir))-vxcir(ir))*(t1-t2)
    end do
! integrate function
    sum=0.d0
    do ir=1,ngrtot
      sum=sum+rfir(ir)*cfunir(ir)
    end do
    sum=sum*omega/dble(ngrtot)
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do irc=1,nrcmt(is)
          fr(irc)=rfmt(1,irc,ias)*rcmt(irc,is)**2
        end do
        call fderiv(-1,nrcmt(is),rcmt(1,is),fr,gr,cf)
        sum=sum+fourpi*y00*gr(nrcmt(is))
      end do
    end do
    delta(ist,jst)=sum
    delta(jst,ist)=-sum
  end do
end do
deallocate(rfmt,rfir,wfmt,wfir)
deallocate(apwalm,evecfv,evecsv)
return
end subroutine

