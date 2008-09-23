
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phveff
use modmain
implicit none
! local variables
integer is,ia,ias,jas,i
integer i1,i2,i3,ir
real(8) v1(3),v2(3)
! external functions
real(8) rfirvec
external rfirvec
! muffin-tin density
ias=0
jas=0
do is=1,nspecies
  do ia=1,natoms0(is)
    ias=ias+1
    do i=1,nphcell
      jas=jas+1
      veffmt(:,:,jas)=veffmt0(:,:,ias)
    end do
  end do
end do
! interstitial density
ir=0
do i3=0,ngrid(3)-1
  v1(3)=dble(i3)/dble(ngrid(3))
  do i2=0,ngrid(2)-1
    v1(2)=dble(i2)/dble(ngrid(2))
    do i1=0,ngrid(1)-1
      v1(1)=dble(i1)/dble(ngrid(1))
      ir=ir+1
      call r3mv(avec,v1,v2)
      veffir(ir)=rfirvec(ngrid0,ainv0,v2,veffir0)
    end do
  end do
end do
call genveffig
return
end subroutine

