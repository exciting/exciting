
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genshtmat
! !INTERFACE:
subroutine genshtmat
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the forward and backward spherical harmonic transformation (SHT)
!   matrices using the spherical covering set produced by the routine
!   {\tt sphcover}. These matrices are used to transform a function between its
!   $(l,m)$-expansion coefficients and its values at the $(\theta,\phi)$ points
!   on the sphere.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer itp,lwork,info
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: tp(:,:)
real(8), allocatable :: rlm(:)
real(8), allocatable :: work(:)
complex(8), allocatable :: ylm(:)
complex(8), allocatable :: zwork(:)
allocate(tp(2,lmmaxapw))
allocate(rlm(lmmaxapw))
allocate(ylm(lmmaxapw))
allocate(ipiv(lmmaxapw))
lwork=2*lmmaxapw
allocate(work(lwork))
allocate(zwork(lwork))
! allocate real SHT matrices for lmaxapw
if (allocated(rbshtapw)) deallocate(rbshtapw)
allocate(rbshtapw(lmmaxapw,lmmaxapw))
if (allocated(rfshtapw)) deallocate(rfshtapw)
allocate(rfshtapw(lmmaxapw,lmmaxapw))
! allocate real SHT matrices for lmaxvr
if (allocated(rbshtvr)) deallocate(rbshtvr)
allocate(rbshtvr(lmmaxvr,lmmaxvr))
if (allocated(rfshtvr)) deallocate(rfshtvr)
allocate(rfshtvr(lmmaxvr,lmmaxvr))
! allocate complex SHT matrices for lmaxapw
if (allocated(zbshtapw)) deallocate(zbshtapw)
allocate(zbshtapw(lmmaxapw,lmmaxapw))
if (allocated(zfshtapw)) deallocate(zfshtapw)
allocate(zfshtapw(lmmaxapw,lmmaxapw))
! allocate complex SHT matrices for lmaxvr
if (allocated(zbshtvr)) deallocate(zbshtvr)
allocate(zbshtvr(lmmaxvr,lmmaxvr))
if (allocated(zfshtvr)) deallocate(zfshtvr)
allocate(zfshtvr(lmmaxvr,lmmaxvr))
! generate spherical covering set for lmaxapw
call sphcover(lmmaxapw,tp)
! generate real and complex spherical harmonics and set the backward SHT arrays
do itp=1,lmmaxapw
  call genrlm(lmaxapw,tp(:,itp),rlm)
  rbshtapw(itp,1:lmmaxapw)=rlm(1:lmmaxapw)
  call genylm(lmaxapw,tp(:,itp),ylm)
  zbshtapw(itp,1:lmmaxapw)=ylm(1:lmmaxapw)
end do
! generate spherical covering set for lmaxvr
call sphcover(lmmaxvr,tp)
! generate real and complex spherical harmonics and set the backward SHT arrays
do itp=1,lmmaxvr
  call genrlm(lmaxvr,tp(:,itp),rlm)
  rbshtvr(itp,1:lmmaxvr)=rlm(1:lmmaxvr)
  call genylm(lmaxvr,tp(:,itp),ylm)
  zbshtvr(itp,1:lmmaxvr)=ylm(1:lmmaxvr)
end do
! find the forward SHT arrays
! real, lmaxapw
rfshtapw(:,:)=rbshtapw(:,:)
call dgetrf(lmmaxapw,lmmaxapw,rfshtapw,lmmaxapw,ipiv,info)
if (info.ne.0) goto 10
call dgetri(lmmaxapw,rfshtapw,lmmaxapw,ipiv,work,lwork,info)
if (info.ne.0) goto 10
! real, lmaxvr
rfshtvr(:,:)=rbshtvr(:,:)
call dgetrf(lmmaxvr,lmmaxvr,rfshtvr,lmmaxvr,ipiv,info)
if (info.ne.0) goto 10
call dgetri(lmmaxvr,rfshtvr,lmmaxvr,ipiv,work,lwork,info)
if (info.ne.0) goto 10
! complex, lmaxapw
zfshtapw(:,:)=zbshtapw(:,:)
call zgetrf(lmmaxapw,lmmaxapw,zfshtapw,lmmaxapw,ipiv,info)
if (info.ne.0) goto 10
call zgetri(lmmaxapw,zfshtapw,lmmaxapw,ipiv,zwork,lwork,info)
if (info.ne.0) goto 10
! complex, lmaxvr
zfshtvr(:,:)=zbshtvr(:,:)
call zgetrf(lmmaxvr,lmmaxvr,zfshtvr,lmmaxvr,ipiv,info)
if (info.ne.0) goto 10
call zgetri(lmmaxvr,zfshtvr,lmmaxvr,ipiv,zwork,lwork,info)
if (info.ne.0) goto 10
deallocate(tp,rlm,ylm,ipiv,work,zwork)
return
10 continue
write(*,*)
write(*,'("Error(genshtmat): unable to find inverse spherical harmonic &
 &transform")')
write(*,'(" => improper spherical covering")')
write(*,*)
stop
end subroutine
!EOC

