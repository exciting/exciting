
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepmain
use modmain
implicit none
! local variables
integer is,ia,ias,ik
integer ir,irc,it,idm
real(8) tau,resp
! allocatable arrays
real(8), allocatable :: rfmt(:,:,:)
real(8), allocatable :: rfir(:)
real(8), allocatable :: rvfmt(:,:,:,:)
real(8), allocatable :: rvfir(:,:)
complex(8), allocatable :: vnlcv(:,:,:,:)
complex(8), allocatable :: vnlvv(:,:,:)
complex(8), allocatable :: zvxmt(:,:,:)
complex(8), allocatable :: zvxir(:)
complex(8), allocatable :: dvxmt(:,:,:)
complex(8), allocatable :: dvxir(:)
complex(8), allocatable :: zbxmt(:,:,:,:)
complex(8), allocatable :: zbxir(:,:)
complex(8), allocatable :: dbxmt(:,:,:,:)
complex(8), allocatable :: dbxir(:,:)
complex(8), allocatable :: zflm(:)
! external functions
real(8) rfinp
external rfinp
if (iscl.lt.1) return
! calculate nonlocal matrix elements
allocate(vnlcv(ncrmax,natmtot,nstsv,nkpt))
allocate(vnlvv(nstsv,nstsv,nkpt))
call oepvnl(vnlcv,vnlvv)
! allocate local arrays
allocate(rfmt(lmmaxvr,nrmtmax,natmtot))
allocate(rfir(ngrtot))
allocate(zvxmt(lmmaxvr,nrcmtmax,natmtot))
allocate(zvxir(ngrtot))
allocate(dvxmt(lmmaxvr,nrcmtmax,natmtot))
allocate(dvxir(ngrtot))
allocate(zflm(lmmaxvr))
if (spinpol) then
  allocate(rvfmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(rvfir(ngrtot,ndmag))
  allocate(zbxmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
  allocate(zbxir(ngrtot,ndmag))
  allocate(dbxmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
  allocate(dbxir(ngrtot,ndmag))
end if
! zero the potential
zvxmt(:,:,:)=0.d0
zvxir(:)=0.d0
if (spinpol) then
  zbxmt(:,:,:,:)=0.d0
  zbxir(:,:)=0.d0
end if
resp=0.d0
! initial step size
tau=tau0oep
! start iteration loop
do it=1,maxitoep
  if (mod(it,10).eq.0) then
    write(*,'("Info(oepmain): done ",I6," iterations of ",I6)') it,maxitoep
  end if
! zero the residues
  dvxmt(:,:,:)=0.d0
  dvxir(:)=0.d0
  if (spinpol) then
    dbxmt(:,:,:,:)=0.d0
    dbxir(:,:)=0.d0
  end if
! calculate the k-dependent residues
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ik=1,nkpt
    call oepresk(ik,vnlcv,vnlvv,zvxmt,zvxir,zbxmt,zbxir,dvxmt,dvxir,dbxmt,dbxir)
  end do
!$OMP END DO
!$OMP END PARALLEL
! compute the real residues
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
        call zflmconj(lmaxvr,dvxmt(1,irc,ias),zflm)
        zflm(:)=zflm(:)+dvxmt(:,irc,ias)
        call ztorflm(lmaxvr,zflm,rfmt(1,ir,ias))
        do idm=1,ndmag
          call zflmconj(lmaxvr,dbxmt(1,irc,ias,idm),zflm)
          zflm(:)=zflm(:)+dbxmt(:,irc,ias,idm)
          call ztorflm(lmaxvr,zflm,rvfmt(1,ir,ias,idm))
        end do
      end do
    end do
  end do
  rfir(:)=2.d0*dble(dvxir(:))
  do idm=1,ndmag
    rvfir(:,idm)=2.d0*dble(dbxir(:,idm))
  end do
! symmetrise the residues
  call symrf(lradstp,rfmt,rfir)
  if (spinpol) call symrvf(lradstp,rvfmt,rvfir)
! magnitude of residues
  resoep=sqrt(abs(rfinp(lradstp,rfmt,rfmt,rfir,rfir)))
  do idm=1,ndmag
    resoep=resoep+sqrt(abs(rfinp(lradstp,rvfmt(1,1,1,idm),rvfmt(1,1,1,idm), &
     rvfir(1,idm),rvfir(1,idm))))
  end do
  resoep=resoep/omega
  if (it.gt.1) then
    if (resoep.lt.resp) then
      tau=tau+dtauoep
    else
      tau=tau0oep
    end if
    tau=min(tau,100.d0)
  end if
  resp=resoep
!--------------------------------------------!
!     update complex potential and field     !
!--------------------------------------------!
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
        call rtozflm(lmaxvr,rfmt(1,ir,ias),zflm)
        zvxmt(:,irc,ias)=zvxmt(:,irc,ias)-tau*zflm(:)
        do idm=1,ndmag
          call rtozflm(lmaxvr,rvfmt(1,ir,ias,idm),zflm)
          zbxmt(:,irc,ias,idm)=zbxmt(:,irc,ias,idm)-tau*zflm(:)
        end do
      end do
    end do
  end do
  zvxir(:)=zvxir(:)-tau*rfir(:)
  do idm=1,ndmag
    zbxir(:,idm)=zbxir(:,idm)-tau*rvfir(:,idm)
  end do
! end iteration loop
end do
! generate the real potential and field
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      call ztorflm(lmaxvr,zvxmt(1,irc,ias),rfmt(1,ir,ias))
      do idm=1,ndmag
        call ztorflm(lmaxvr,zbxmt(1,irc,ias,idm),rvfmt(1,ir,ias,idm))
      end do
    end do
  end do
end do
! convert potential and field from a coarse to a fine radial mesh
call rfmtctof(rfmt)
do idm=1,ndmag
  call rfmtctof(rvfmt(1,1,1,idm))
end do
! add to existing correlation potential and field
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      vxcmt(:,ir,ias)=vxcmt(:,ir,ias)+rfmt(:,ir,ias)
      do idm=1,ndmag
        bxcmt(:,ir,ias,idm)=bxcmt(:,ir,ias,idm)+rvfmt(:,ir,ias,idm)
      end do
    end do
  end do
end do
vxcir(:)=vxcir(:)+dble(zvxir(:))
do idm=1,ndmag
  bxcir(:,idm)=bxcir(:,idm)+dble(zbxir(:,idm))
end do
deallocate(rfmt,rfir,vnlcv,vnlvv)
deallocate(zvxmt,zvxir,dvxmt,dvxir,zflm)
if (spinpol) then
  deallocate(rvfmt,rvfir)
  deallocate(zbxmt,zbxir,dbxmt,dbxir)
end if
return
end subroutine
