
! Copyright (C) 2002-2005 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine moke
use modmain
implicit none
! local variables
integer iw,iostat
complex(8) zt1,zt2,zt3
! allocatable arrays
real(8), allocatable :: w(:)
real(8), allocatable :: sig1(:,:)
real(8), allocatable :: sig2(:,:)
complex(8), allocatable :: kerr(:)
! calculate dielectric function for the 11 and 12 components
noptcomp=2
optcomp(1,1)=1
optcomp(2,1)=1
optcomp(1,2)=1
optcomp(2,2)=2
call dielectric
! allocate local arrays
allocate(w(nwdos))
allocate(sig1(nwdos,2),sig2(nwdos,2))
allocate(kerr(nwdos))
! read diagonal contribution to optical conductivity
open(50,file='SIGMA_11.OUT',action='READ',status='OLD', &
 form='FORMATTED',iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(moke): error opening SIGMA_11.OUT")')
  write(*,*)
  stop
end if
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sig1(iw,1)
end do
read(50,*)
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sig2(iw,1)
end do
close(50)
! read off-diagonal contribution to optical conductivity
open(50,file='SIGMA_12.OUT',action='READ',status='OLD', &
 form='FORMATTED',iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(moke): error opening SIGMA_12.OUT")')
  write(*,*)
  stop
end if
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sig1(iw,2)
end do
read(50,*)
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sig2(iw,2)
end do
close(50)
! calculate the complex Kerr angle
do iw=1,nwdos
  if (w(iw).gt.0.d0) then
    zt1=cmplx(sig1(iw,1),sig2(iw,1),8)
    zt2=cmplx(sig1(iw,2),sig2(iw,2),8)
    zt3=zt1*sqrt(1.d0+fourpi*zi*zt1/w(iw))
    if (abs(zt3).gt.1.d-8) then
      kerr(iw)=-zt2/zt3
    else
      kerr(iw)=0.d0
    end if
  else
    kerr(iw)=0.d0
  end if
end do
open(50,file='KERR.OUT',action='WRITE',form='FORMATTED')
do iw=1,nwdos
  write(50,'(2G18.10)') w(iw),dble(kerr(iw))*180.d0/pi
end do
write(50,'("     ")')
do iw=1,nwdos
  write(50,'(2G18.10)') w(iw),aimag(kerr(iw))*180.d0/pi
end do
close(50)
write(*,*)
write(*,'("Info(moke):")')
write(*,'(" complex Kerr angle in degrees written to KERR.OUT")')
write(*,*)
deallocate(w,sig1,sig2,kerr)
return
end subroutine

