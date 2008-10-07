
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengroup(ngen,srgen,stgen,ngrp,srgrp,stgrp)
implicit none
! arguments
integer, intent(in) :: ngen
real(8), intent(in) :: srgen(3,3,ngen)
real(8), intent(in) :: stgen(3,ngen)
integer, intent(out) :: ngrp
real(8), intent(out) :: srgrp(3,3,192)
real(8), intent(out) :: stgrp(3,192)
! local variables
integer i,j,k
real(8), parameter :: eps=1.d-6
real(8) sr(3,3),st(3)
! external functions
logical seitzeq
external seitzeq
! store the identity
ngrp=1
srgrp(1,1,1)=1.d0; srgrp(1,2,1)=0.d0; srgrp(1,3,1)=0.d0
srgrp(2,1,1)=0.d0; srgrp(2,2,1)=1.d0; srgrp(2,3,1)=0.d0
srgrp(3,1,1)=0.d0; srgrp(3,2,1)=0.d0; srgrp(3,3,1)=1.d0
stgrp(:,1)=0.d0
10 continue
! right multiply by the generators
do i=1,ngen
  do j=1,ngrp
    call seitzmul(eps,srgrp(:,:,j),stgrp(:,j),srgen(:,:,i),stgen(:,i),sr,st)
! check if the new element already exists
    do k=1,ngrp
      if (seitzeq(eps,srgrp(:,:,k),stgrp(:,k),sr,st)) goto 20
    end do
    goto 40
20 continue
  end do
end do
! left multiply by the generators
do i=1,ngen
  do j=1,ngrp
    call seitzmul(eps,srgen(:,:,i),stgen(:,i),srgrp(:,:,j),stgrp(:,j),sr,st)
! check if the new element already exists
    do k=1,ngrp
      if (seitzeq(eps,srgrp(:,:,k),stgrp(:,k),sr,st)) goto 30
    end do
    goto 40
30 continue
  end do
end do
! all elements accounted for
return
40 continue
! add new element
ngrp=ngrp+1
if (ngrp.gt.192) then
  write(*,*)
  write(*,'("Error(gengroup): more than 192 group elements")')
  write(*,*)
  stop
end if
srgrp(:,:,ngrp)=sr(:,:)
stgrp(:,ngrp)=st(:)
goto 10
return
end subroutine

