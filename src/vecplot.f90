
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: vecplot
! !INTERFACE:
subroutine vecplot
! !DESCRIPTION:
!   Outputs a 2D or 3D vector field for plotting. The vector field can be the
!   magnetisation vector field, ${\bf m}$; the exchange-correlation magnetic
!   field, ${\bf B}_{\rm xc}$; or the electric field
!   ${\bf E}\equiv-\nabla V_{\rm C}$. The magnetisation is obtained from the
!   spin density matrix, $\rho_{\alpha\beta}$, by solving
!   $$ \rho_{\alpha\beta}({\bf r})=\frac{1}{2}\left(n({\bf r})
!    \delta_{\alpha\beta}+\sigma\cdot {\bf m}({\bf r})\right), $$
!   where $n\equiv\tr\rho_{\alpha\beta}$ is the total density. In the case of 2D
!   plots, the magnetisation vectors are still 3D, but are in the coordinate
!   system of the plane.
!
! !REVISION HISTORY:
!   Created August 2004 (JKD)
!   Included electric field plots, August 2006 (JKD)
!EOP
!BOC
use modmain
implicit none
! local variables
integer is,ia,ias,ir,lm
real(8) vl1(3),vl2(3),vc1(3),vc2(3),vc3(3),vc4(3),t1
! allocatable arrays
real(8), allocatable :: rvfmt(:,:,:,:)
real(8), allocatable :: rvfir(:,:)
! external functions
real(8) r3dot
external r3dot
if ((task.eq.72).or.(task.eq.73).or.(task.eq.82).or.(task.eq.83)) then
  if (.not.spinpol) then
    write(*,*)
    write(*,'("Error(vecplot): spin-unpolarised magnetisation/field is zero")')
    write(*,*)
    stop
  end if
end if
! initialise universal variables
call init0
! read magnetisation from file
call readstate
allocate(rvfmt(lmmaxvr,nrmtmax,natmtot,3))
allocate(rvfir(ngrtot,3))
select case(task)
case(72,73)
! magnetisation
  if (ndmag.eq.3) then
! non-collinear
    rvfmt(:,:,:,:)=magmt(:,:,:,:)
    rvfir(:,:)=magir(:,:)
  else
! collinear
    rvfmt(:,:,:,1:2)=0.d0
    rvfir(:,1:2)=0.d0
    rvfmt(:,:,:,3)=magmt(:,:,:,1)
    rvfir(:,3)=magir(:,1)
  end if
case(82,83)
! effective magnetic field
  if (ndmag.eq.3) then
! non-collinear
    rvfmt(:,:,:,:)=bxcmt(:,:,:,:)
    rvfir(:,:)=bxcir(:,:)
  else
! collinear
    rvfmt(:,:,:,1:2)=0.d0
    rvfir(:,1:2)=0.d0
    rvfmt(:,:,:,3)=bxcmt(:,:,:,1)
    rvfir(:,3)=bxcir(:,1)
  end if
case(142,143)
! electric field
  call gradrf(vclmt,vclir,rvfmt,rvfir)
! use the negative of the gradient
  rvfmt(:,:,:,:)=-rvfmt(:,:,:,:)
  rvfir(:,:)=-rvfir(:,:)
case(152,153)
  if (ndmag.lt.3) then
    write(*,*)
    write(*,'("Error(vecplot): collinear m(r)xB(r) is zero")')
    write(*,*)
    stop
  end if
  call rvfcross(magmt,bxcmt,magir,bxcir,rvfmt,rvfir)
end select
select case(task)
case(72,82,142,152)
! determine the projection of the magnetisation/field onto the plotting plane
  vl1(:)=vclp2d(:,2)-vclp2d(:,1)
  vl2(:)=vclp2d(:,3)-vclp2d(:,1)
  call r3mv(avec,vl1,vc1)
  call r3mv(avec,vl2,vc2)
  t1=sqrt(vc1(1)**2+vc1(2)**2+vc1(3)**2)
  vc1(:)=vc1(:)/t1
  t1=sqrt(vc2(1)**2+vc2(2)**2+vc2(3)**2)
  vc2(:)=vc2(:)/t1
  call r3cross(vc1,vc2,vc3)
  t1=sqrt(vc3(1)**2+vc3(2)**2+vc3(3)**2)
  vc3(:)=vc3(:)/t1
  call r3cross(vc3,vc1,vc2)
  t1=sqrt(vc2(1)**2+vc2(2)**2+vc2(3)**2)
  vc2(:)=vc2(:)/t1
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nrmt(is)
        do lm=1,lmmaxvr
          vc4(:)=rvfmt(lm,ir,ias,:)
          rvfmt(lm,ir,ias,1)=r3dot(vc4,vc1)
          rvfmt(lm,ir,ias,2)=r3dot(vc4,vc2)
          rvfmt(lm,ir,ias,3)=r3dot(vc4,vc3)
        end do
      end do
    end do
  end do
  do ir=1,ngrtot
    vc4(:)=rvfir(ir,:)
    rvfir(ir,1)=r3dot(vc4,vc1)
    rvfir(ir,2)=r3dot(vc4,vc2)
    rvfir(ir,3)=r3dot(vc4,vc3)
  end do
  if (task.eq.72) then
    open(50,file='MAG2D.OUT',action='WRITE',form='FORMATTED')
  else if (task.eq.82) then
    open(50,file='BXC2D.OUT',action='WRITE',form='FORMATTED')
  else if (task.eq.142) then
    open(50,file='EF2D.OUT',action='WRITE',form='FORMATTED')
  else
    open(50,file='MCBXC2D.OUT',action='WRITE',form='FORMATTED')
  end if
  call plot2d(50,3,lmaxvr,lmmaxvr,rvfmt,rvfir)
  close(50)
  write(*,*)
  write(*,'("Info(vecplot):")')
  if (task.eq.72) then
    write(*,'(" 2D magnetisation density written to MAG2D.OUT")')
  else if (task.eq.82) then
    write(*,'(" 2D exchange-correlation field written to BXC2D.OUT")')
  else if (task.eq.142) then
    write(*,'(" 2D electric field written to EF2D.OUT")')
  else
    write(*,'(" 2D m(r) x B_xc(r) written to MCBXC2D.OUT")')
  end if
  write(*,*)
case(73,83,143,153)
  if (task.eq.73) then
    open(50,file='MAG3D.OUT',action='WRITE',form='FORMATTED')
  else if (task.eq.83) then
    open(50,file='BXC3D.OUT',action='WRITE',form='FORMATTED')
  else if (task.eq.143) then
    open(50,file='EF3D.OUT',action='WRITE',form='FORMATTED')
  else
    open(50,file='MCBXC3D.OUT',action='WRITE',form='FORMATTED')
  end if
  call plot3d(50,3,lmaxvr,lmmaxvr,rvfmt,rvfir)
  close(50)
  write(*,*)
  write(*,'("Info(vecplot):")')
  if (task.eq.73) then
    write(*,'(" 3D magnetisation density written to MAG3D.OUT")')
  else if (task.eq.83) then
    write(*,'(" 3D exchange-correlation field written to BXC3D.OUT")')
  else if (task.eq.143) then
    write(*,'(" 3D electric field written to EF3D.OUT")')
  else
    write(*,'(" 3D m(r) x B_xc(r) written to MCBXC3D.OUT")')
  end if
  write(*,*)
end select
deallocate(rvfmt,rvfir)
return
end subroutine
!EOC

