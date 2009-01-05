
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

real(8) function rfirvec(ngrid,ainv,vc,rfir)
implicit none
! arguments
integer, intent(in) :: ngrid(3)
real(8), intent(in) :: ainv(3,3)
real(8), intent(in) :: vc(3)
real(8), intent(in) :: rfir(ngrid(1),ngrid(2),ngrid(3))
! local variables
integer i1,i2,i3,j1,j2,j3
real(8), parameter :: eps=1.d-6
real(8) v1,v2,v3,p1,p2,p3,q1,q2,q3
real(8) f00,f01,f10,f11,f0,f1
! input vector in lattice coordinates
v1=ainv(1,1)*vc(1)+ainv(1,2)*vc(2)+ainv(1,3)*vc(3)
v2=ainv(2,1)*vc(1)+ainv(2,2)*vc(2)+ainv(2,3)*vc(3)
v3=ainv(3,1)*vc(1)+ainv(3,2)*vc(2)+ainv(3,3)*vc(3)
! map lattice coordinates to [0,1) interval
i1=int(v1)
i2=int(v2)
i3=int(v3)
v1=v1-dble(i1)
v2=v2-dble(i2)
v3=v3-dble(i3)
if (v1.lt.0.d0) v1=v1+1.d0
if (v2.lt.0.d0) v2=v2+1.d0
if (v3.lt.0.d0) v3=v3+1.d0
if (1.d0-v1.lt.eps) v1=0.d0
if (1.d0-v2.lt.eps) v2=0.d0
if (1.d0-v3.lt.eps) v3=0.d0
if (v1.lt.eps) v1=0.d0
if (v2.lt.eps) v2=0.d0
if (v3.lt.eps) v3=0.d0
! determine coordinates on grid
v1=dble(ngrid(1))*v1
v2=dble(ngrid(2))*v2
v3=dble(ngrid(3))*v3
i1=int(v1)
i2=int(v2)
i3=int(v3)
! use trilinear interpolation with neighbouring points
p1=v1-dble(i1)
p2=v2-dble(i2)
p3=v3-dble(i3)
q1=1.d0-p1
q2=1.d0-p2
q3=1.d0-p3
i1=modulo(i1,ngrid(1))+1
i2=modulo(i2,ngrid(2))+1
i3=modulo(i3,ngrid(3))+1
j1=modulo(i1,ngrid(1))+1
j2=modulo(i2,ngrid(2))+1
j3=modulo(i3,ngrid(3))+1
f00=rfir(i1,i2,i3)*q1+rfir(j1,i2,i3)*p1
f01=rfir(i1,i2,j3)*q1+rfir(j1,i2,j3)*p1
f10=rfir(i1,j2,i3)*q1+rfir(j1,j2,i3)*p1
f11=rfir(i1,j2,j3)*q1+rfir(j1,j2,j3)*p1
f0=f00*q2+f10*p2
f1=f01*q2+f11*p2
rfirvec=f0*q3+f1*p3
return
end function

