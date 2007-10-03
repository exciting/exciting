
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rschroddme
! !INTERFACE:
subroutine rschroddme(m,l,k,e,np,nr,r,vr,nn,p0,p1,q0,q1)
! !INPUT/OUTPUT PARAMETERS:
!   m   : order of energy derivative (in,integer)
!   l   : angular momentum quantum number (in,integer)
!   k   : quantum number k, zero if Dirac eqn. is not to be used (in,integer)
!   e   : energy (in,real)
!   np  : order of predictor-corrector polynomial (in,integer)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   nn  : number of nodes (out,integer)
!   p0  : m th energy derivative of P (out,real(nr))
!   p1  : radial derivative of p0 (out,real(nr))
!   q0  : m th energy derivative of Q (out,real(nr))
!   q1  : radial derivative of q0 (out,real(nr))
! !DESCRIPTION:
!   Finds the solution to the $m$th energy derivative of the scalar relativistic
!   radial Schr\"{o}dinger equation using the routine {\tt rschrodint}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
integer, intent(in) :: m
integer, intent(in) :: l
integer, intent(in) :: k
real(8), intent(in) :: e
integer, intent(in) :: np
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: vr(nr)
integer, intent(out) :: nn
real(8), intent(out) :: p0(nr)
real(8), intent(out) :: p1(nr)
real(8), intent(out) :: q0(nr)
real(8), intent(out) :: q1(nr)
! local variables
integer im,kpa,ir
! fine-structure constant
real(8), parameter :: alpha=1.d0/137.03599911d0
real(8) rm
! allocatable arrays
real(8), allocatable :: p0p(:)
real(8), allocatable :: g0(:),g1(:)
real(8), allocatable :: f0(:),f1(:)
real(8), allocatable :: cf(:,:)
if (nr.le.0) then
  write(*,*)
  write(*,'("Error(rschroddme): invalid nr : ",I8)') nr
  write(*,*)
  stop
end if
if ((m.lt.0).or.(m.gt.6)) then
  write(*,*)
  write(*,'("Error(rschroddme): m out of range : ",I8)') m
  write(*,*)
  stop
end if
if (k.eq.0) then
! use the scalar relativistic Schrodinger equation
  allocate(p0p(nr))
  if (m.eq.0) then
    call rschrodint(m,l,e,np,nr,r,vr,nn,p0p,p0,p1,q0,q1)
  else
    do im=0,m
      call rschrodint(im,l,e,np,nr,r,vr,nn,p0p,p0,p1,q0,q1)
      p0p(:)=p0(:)
    end do
  end if
  deallocate(p0p)
else
! use the Dirac equation
  allocate(g0(nr),g1(nr))
  allocate(f0(nr),f1(nr))
  allocate(cf(3,nr))
  if (k.eq.l) then
    kpa=k
  else if (k.eq.l+1) then
    kpa=-k
  else
    write(*,*)
    write(*,'("Error(rschroddme): incompatible l and k : ",2I8)') l,k
    write(*,*)
    stop
  end if
  call rdiracdme(m,kpa,e,np,nr,r,vr,nn,g0,g1,f0,f1)
! determine equivalent scalar relativistic functions
  do ir=1,nr
    rm=1.d0-0.5d0*(alpha**2)*vr(ir)
    p0(ir)=g0(ir)
    p1(ir)=g1(ir)
    q0(ir)=(p1(ir)-p0(ir)/r(ir))/(2.d0*rm)
  end do
  call fderiv(1,nr,r,q0,q1,cf)
  deallocate(g0,g1,f0,f1,cf)
end if
return
end subroutine
!EOC

