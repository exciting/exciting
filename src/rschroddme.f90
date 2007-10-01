
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rschroddme
! !INTERFACE:
subroutine rschroddme(m,l,e,np,nr,r,vr,nn,p0,p1,q0,q1)
! !INPUT/OUTPUT PARAMETERS:
!   m   : order of energy derivative (in,integer)
!   l   : angular momentum quantum number (in,integer)
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
integer im
! automatic arrays
real(8) p0p(nr)
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
if (m.eq.0) then
  call rschrodint(m,l,e,np,nr,r,vr,nn,p0p,p0,p1,q0,q1)
else
  do im=0,m
    call rschrodint(im,l,e,np,nr,r,vr,nn,p0p,p0,p1,q0,q1)
    p0p(:)=p0(:)
  end do
end if
return
end subroutine
!EOC
