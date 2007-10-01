
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
! allocate global SHT arrays
if (allocated(rbshtapw)) deallocate(rbshtapw)
allocate(rbshtapw(lmmaxapw,lmmaxapw))
if (allocated(rfshtvr)) deallocate(rfshtvr)
allocate(rfshtvr(lmmaxvr,lmmaxvr))
if (allocated(zbshtapw)) deallocate(zbshtapw)
allocate(zbshtapw(lmmaxapw,lmmaxapw))
if (allocated(zfshtapw)) deallocate(zfshtapw)
allocate(zfshtapw(lmmaxapw,lmmaxapw))
if (allocated(zfshtvr)) deallocate(zfshtvr)
allocate(zfshtvr(lmmaxvr,lmmaxvr))
! generate spherical covering set
call sphcover(lmmaxapw,tp)
! generate real and complex spherical harmonics and set the backward SHT arrays
do itp=1,lmmaxapw
  call genrlm(lmaxapw,tp(1,itp),rlm)
  rbshtapw(itp,1:lmmaxapw)=rlm(1:lmmaxapw)
  call genylm(lmaxapw,tp(1,itp),ylm)
  zbshtapw(itp,1:lmmaxapw)=ylm(1:lmmaxapw)
end do
! find the forward SHT arrays
! real, lmaxvr
rfshtvr(1:lmmaxvr,1:lmmaxvr)=rbshtapw(1:lmmaxvr,1:lmmaxvr)
call dgetrf(lmmaxvr,lmmaxvr,rfshtvr,lmmaxvr,ipiv,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(genshtmat): dgetrf returned non-zero info : ",I8)') info
  write(*,'(" => improper spherical covering for lmaxvr")')
  write(*,*)
  stop
end if
call dgetri(lmmaxvr,rfshtvr,lmmaxvr,ipiv,work,lwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(genshtmat): dgetri returned non-zero info : ",I8)') info
  write(*,'(" => improper spherical covering for lmaxvr")')
  write(*,*)
  stop
end if
! complex, lmaxapw
zfshtapw(:,:)=zbshtapw(:,:)
call zgetrf(lmmaxapw,lmmaxapw,zfshtapw,lmmaxapw,ipiv,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(genshtmat): zgetrf returned non-zero info : ",I8)') info
  write(*,'(" => improper spherical covering for lmaxapw")')
  write(*,*)
  stop
end if
call zgetri(lmmaxapw,zfshtapw,lmmaxapw,ipiv,zwork,lwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(genshtmat): zgetri returned non-zero info : ",I8)') info
  write(*,'(" => improper spherical covering for lmaxapw")')
  write(*,*)
  stop
end if
! complex, lmaxvr
zfshtvr(1:lmmaxvr,1:lmmaxvr)=zbshtapw(1:lmmaxvr,1:lmmaxvr)
call zgetrf(lmmaxvr,lmmaxvr,zfshtvr,lmmaxvr,ipiv,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(genshtmat): zgetrf returned non-zero info : ",I8)') info
  write(*,'(" => improper spherical covering for lmaxvr")')
  write(*,*)
  stop
end if
call zgetri(lmmaxvr,zfshtvr,lmmaxvr,ipiv,zwork,lwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(genshtmat): zgetri returned non-zero info : ",I8)') info
  write(*,'(" => improper spherical covering for lmaxvr")')
  write(*,*)
  stop
end if
deallocate(tp,rlm,ylm,ipiv,work,zwork)
return
end subroutine
!EOC

