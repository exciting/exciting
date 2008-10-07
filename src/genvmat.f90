
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmat(vmt,vir,vmat)
! generates potential matrix elements for all states and k-points
use modmain
implicit none
! arguments
real(8), intent(in) :: vmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: vir(ngrtot)
complex(8), intent(out) :: vmat(nstsv,nstsv,nkpt)
! local variables
integer is,ia,ias,irc,ir
integer ik
! local arrays
real(8), allocatable :: rfmt(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
! allocate local arrays
allocate(rfmt(lmmaxvr,nrcmtmax,natmtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor,nstsv))
! convert muffin-tin potential to spherical coordinates
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,vmt(:,ir,ias),1, &
       0.d0,rfmt(:,irc,ias),1)
    end do
  end do
end do
! loop over k-points
do ik=1,nkpt
! get the eigenvectors and values from file
  call getevalsv(vkl(:,ik),evalsv)
  call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! calculate the wavefunctions for all states
  call genwfsv(.false.,ngk(1,ik),igkig(:,1,ik),evalsv,apwalm,evecfv,evecsv, &
   wfmt,wfir)
  call genvmatk(rfmt,vir,wfmt,wfir,vmat(:,:,ik))
end do
deallocate(apwalm,evecfv,evecsv,wfmt,wfir)
return
end subroutine

