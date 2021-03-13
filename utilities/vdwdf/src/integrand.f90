! Copyright (C) 2007-2010 D. Nabok, P. Puschnig and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!============================================================
! This subroutine gives the integrand function for
! using in Cuba program
!
! Integrand = n(r1)*phi(r1,r2)*n(r2)
!============================================================
integer function integrand(ndim, x, ncomp, func)
!-----------------------------------------
use param
!-----------------------------------------

implicit none
! input
integer          :: ndim         ! number of dimensions of the integral
integer          :: ncomp        ! number of components of the integral
double precision :: x(ndim)      ! point in ndim space
! output
double precision :: func(ncomp)  ! value of integrand in point x
!-----------------------------------------

real*8  :: xsc(6), xpr(6), r1(3), r2(3), dr(3)
real*8  :: dis
real*8  :: dens1, dens2
real*8  :: graddens1, graddens2
integer :: i

!-------------------------------------------------|
!------------------ INTEGRAND --------------------|
!-------------------------------------------------|

! NOTE THAT: 
!   point x() is in fractional supercell units, i.e, x(i) in [0,1]
!   point r belongs to the primitive unit cell ( r  = [xsc(1),xsc(2),xsc(3)] )
!   point r' runs over full space ( r' = [xsc(4),xsc(5),xsc(6)] )

! unscale r' vector according to the supercell dimensions
   xsc(1)=x(1);    xsc(4)=(2*dble(nrx)+1.0)*x(4)-dble(nrx)
   xsc(2)=x(2);    xsc(5)=(2*dble(nry)+1.0)*x(5)-dble(nry)
   xsc(3)=x(3);    xsc(6)=(2*dble(nrz)+1.0)*x(6)-dble(nrz)
   !write(6,*) 'x=', x
   !write(6,*) 'xsc=', xsc

!  shifting to unitcell
   xpr(1)=xsc(1)
   xpr(2)=xsc(2)
   xpr(3)=xsc(3)
   do i=4,6
     if ( xsc(i).lt.0.0d0 ) then
          xpr(i)=xsc(i)-dble(int(xsc(i)+eps))+1.0d0
     else
          xpr(i)=xsc(i)-dble(int(xsc(i)+eps))
     end if
   end do
   !write(*,*) 'xpr=', xpr
   
!-------------------------------------------
!  density and gradient interpolation
   call blend103(xpr(1:3),dens1,graddens1)
   !write(6,*) 'dens1=', dens1
   !write(6,*) 'graddens1=', graddens1

   call blend103(xpr(4:6),dens2,graddens2)
   !write(6,*) 'dens2=', dens2
   !write(6,*) 'graddens2=', graddens2
 
!-------------------------------------------
!  calculation real space radius vectors
   r1(:)=xsc(1)*vectors(:,1)+xsc(2)*vectors(:,2)+xsc(3)*vectors(:,3)
   r2(:)=xsc(4)*vectors(:,1)+xsc(5)*vectors(:,2)+xsc(6)*vectors(:,3)

   dr   = r2 - r1
   dis  = dsqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
   !write(6,*) 'dis=', dis

!-------------------------------------------
!  Choice of vdW-DF type
  
   if ( (vdWDF_version=='vdW-DF') .or.(vdWDF_version=='vdW-DF2') ) then
   
       call  vdWDF_LL(dis,dens1,graddens1,dens2,graddens2,func)
   
   else if (vdWDF_version=='VV09') then
   
       call  vdWDF_VV(dis,dens1,graddens1,dens2,graddens2,func)
   
   else
   
       write(*,*) 'ERROR(integrand.f90): Unknown vdW-DF version!'
       stop
   
   end if
   
   integrand = 0
   
return
end function integrand

