
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readstate
! !INTERFACE:
subroutine readstate
! !USES:
use modmain
! !DESCRIPTION:
!   Reads in the charge density and other relevant variables from the file
!   {\tt STATE.OUT}. Checks for version and parameter compatibility.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
logical spinpol_
integer is,ia,ias,lmmax,lm,ir,jr
integer idm,ngm,i1,i2,i3,j1,j2,j3
integer version_(3),nspecies_,lmmaxvr_,nrmtmax_
integer natoms_,nrmt_(maxspecies),ngrid_(3)
integer ngrtot_,ngvec_,ndmag_
real(8) t1
! allocatable arrays
integer, allocatable :: mapir(:)
real(8), allocatable :: spr_(:,:)
real(8), allocatable :: rhomt_(:,:,:)
real(8), allocatable :: rhoir_(:)
real(8), allocatable :: vclmt_(:,:,:)
real(8), allocatable :: vclir_(:)
real(8), allocatable :: vxcmt_(:,:,:)
real(8), allocatable :: vxcir_(:)
real(8), allocatable :: veffmt_(:,:,:)
real(8), allocatable :: veffir_(:)
real(8), allocatable :: magmt_(:,:,:,:)
real(8), allocatable :: magir_(:,:)
real(8), allocatable :: bxcmt_(:,:,:,:)
real(8), allocatable :: bxcir_(:,:)
complex(8), allocatable :: veffig_(:)
open(50,file='STATE'//trim(filext),action='READ',form='UNFORMATTED', &
 status='OLD')
read(50) version_
if ((version(1).ne.version_(1)).or.(version(2).ne.version_(2)) &
 .or.(version(3).ne.version_(3))) then
  write(*,*)
  write(*,'("Warning(readstate): different versions")')
  write(*,'(" current   : ",I3.3,".",I3.3,".",I3.3)') version
  write(*,'(" STATE.OUT : ",I3.3,".",I3.3,".",I3.3)') version_
end if
read(50) spinpol_
read(50) nspecies_
if (nspecies.ne.nspecies_) then
  write(*,*)
  write(*,'("Error(readstate): differing nspecies")')
  write(*,'(" current   : ",I4)') nspecies
  write(*,'(" STATE.OUT : ",I4)') nspecies_
  write(*,*)
  stop
end if
read(50) lmmaxvr_
read(50) nrmtmax_
allocate(spr_(nrmtmax_,nspecies))
do is=1,nspecies
  read(50) natoms_
  if (natoms(is).ne.natoms_) then
    write(*,*)
    write(*,'("Error(readstate): differing natoms for species ",I4)') is
    write(*,'(" current   : ",I4)') natoms(is)
    write(*,'(" STATE.OUT : ",I4)') natoms_
    write(*,*)
    stop
  end if
  read(50) nrmt_(is)
  read(50) spr_(1:nrmt_(is),is)
end do
read(50) ngrid_
read(50) ngvec_
read(50) ndmag_
if ((spinpol_).and.(ndmag_.ne.1).and.(ndmag_.ne.3)) then
  write(*,*)
  write(*,'("Error(readstate): invalid ndmag in STATE.OUT : ",I8)') ndmag_
  write(*,*)
  stop
end if
ngrtot_=ngrid_(1)*ngrid_(2)*ngrid_(3)
allocate(mapir(ngrtot))
allocate(rhomt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(rhoir_(ngrtot_))
allocate(vclmt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(vclir_(ngrtot_))
allocate(vxcmt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(vxcir_(ngrtot_))
allocate(veffmt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(veffir_(ngrtot_))
allocate(veffig_(ngvec_))
! read muffin-tin density
read(50) rhomt_,rhoir_
! read Coulomb potential (spin independent)
read(50) vclmt_,vclir_
! read exchange-correlation potential
read(50) vxcmt_,vxcir_
! read effective potential
read(50) veffmt_,veffir_,veffig_
! read magnetisation and effective field
if (spinpol_) then
  allocate(magmt_(lmmaxvr_,nrmtmax_,natmtot,ndmag_))
  allocate(magir_(ngrtot_,ndmag_))
  allocate(bxcmt_(lmmaxvr_,nrmtmax_,natmtot,ndmag_))
  allocate(bxcir_(ngrtot_,ndmag_))
  read(50) magmt_,magir_
  read(50) bxcmt_,bxcir_
end if
close(50)
!---------------------------!
!     muffin-tin arrays     !
!---------------------------!
rhomt(:,:,:)=0.d0
vclmt(:,:,:)=0.d0
vxcmt(:,:,:)=0.d0
veffmt(:,:,:)=0.d0
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  bxcmt(:,:,:,:)=0.d0
end if
lmmax=min(lmmaxvr,lmmaxvr_)
! interpolate the old arrays on the new radial mesh
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do lm=1,lmmax
      call rfinterp(nrmt_(is),spr_(1,is),lmmaxvr_,rhomt_(lm,1,ias),nrmt(is), &
       spr(1,is),lmmaxvr,rhomt(lm,1,ias))
      call rfinterp(nrmt_(is),spr_(1,is),lmmaxvr_,vclmt_(lm,1,ias),nrmt(is), &
       spr(1,is),lmmaxvr,vclmt(lm,1,ias))
      call rfinterp(nrmt_(is),spr_(1,is),lmmaxvr_,vxcmt_(lm,1,ias),nrmt(is), &
       spr(1,is),lmmaxvr,vxcmt(lm,1,ias))
      call rfinterp(nrmt_(is),spr_(1,is),lmmaxvr_,veffmt_(lm,1,ias),nrmt(is), &
       spr(1,is),lmmaxvr,veffmt(lm,1,ias))
    end do
    if ((spinpol).and.(spinpol_)) then
      if (ndmag.eq.ndmag_) then
        do idm=1,ndmag
          do lm=1,lmmax
            call rfinterp(nrmt_(is),spr_(1,is),lmmaxvr_,magmt_(lm,1,ias,idm), &
             nrmt(is),spr(1,is),lmmaxvr,magmt(lm,1,ias,idm))
            call rfinterp(nrmt_(is),spr_(1,is),lmmaxvr_,bxcmt_(lm,1,ias,idm), &
             nrmt(is),spr(1,is),lmmaxvr,bxcmt(lm,1,ias,idm))
          end do
        end do
      else
        do lm=1,lmmax
          call rfinterp(nrmt_(is),spr_(1,is),lmmaxvr_,magmt_(lm,1,ias,ndmag_), &
           nrmt(is),spr(1,is),lmmaxvr,magmt(lm,1,ias,ndmag))
          call rfinterp(nrmt_(is),spr_(1,is),lmmaxvr_,bxcmt_(lm,1,ias,ndmag_), &
           nrmt(is),spr(1,is),lmmaxvr,bxcmt(lm,1,ias,ndmag))
        end do
      end if
    end if
  end do
end do
!-----------------------------!
!     interstitial arrays     !
!-----------------------------!
rhoir(:)=0.d0
vclir(:)=0.d0
vxcir(:)=0.d0
veffir(:)=0.d0
veffig(:)=0.d0
if (spinpol) then
  magir(:,:)=0.d0
  bxcir(:,:)=0.d0
end if
! map from new grid to old
do i3=0,ngrid(3)-1
  t1=dble(ngrid_(3))/dble(ngrid(3))
  j3=nint(t1*dble(i3))
  j3=min(max(j3,0),ngrid_(3)-1)
  do i2=0,ngrid(2)-1
    t1=dble(ngrid_(2))/dble(ngrid(2))
    j2=nint(t1*dble(i2))
    j2=min(max(j2,0),ngrid_(2)-1)
    do i1=0,ngrid(1)-1
      t1=dble(ngrid_(1))/dble(ngrid(1))
      j1=nint(t1*dble(i1))
      j1=min(max(j1,0),ngrid_(1)-1)
      ir=i3*ngrid(2)*ngrid(1)+i2*ngrid(1)+i1+1
      jr=j3*ngrid_(2)*ngrid_(1)+j2*ngrid_(1)+j1+1
      mapir(ir)=jr
    end do
  end do
end do
do ir=1,ngrtot
  jr=mapir(ir)
  rhoir(ir)=rhoir_(jr)
  vclir(ir)=vclir_(jr)
  vxcir(ir)=vxcir_(jr)
  veffir(ir)=veffir_(jr)
end do
ngm=min(ngvec,ngvec_)
veffig(1:ngm)=veffig_(1:ngm)
if ((spinpol).and.(spinpol_)) then
  do ir=1,ngrtot
    jr=mapir(ir)
    if (ndmag.eq.ndmag_) then
      magir(ir,:)=magir_(jr,:)
      bxcir(ir,:)=bxcir_(jr,:)
    else
      magir(ir,ndmag)=magir_(jr,ndmag_)
      bxcir(ir,ndmag)=bxcir_(jr,ndmag_)
    end if
  end do
end if
deallocate(mapir,spr_,rhomt_,rhoir_,vclmt_,vclir_)
deallocate(vxcmt_,vxcir_,veffmt_,veffir_,veffig_)
if (spinpol_) deallocate(magmt_,magir_,bxcmt_,bxcir_)
return
end subroutine
!EOC
