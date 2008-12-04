
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rdiracdme
! !INTERFACE:
subroutine rdiracdme(m,kpa,e,np,nr,r,vr,nn,g0,g1,f0,f1)
! !INPUT/OUTPUT PARAMETERS:
!   m   : order of energy derivative (in,integer)
!   kpa : quantum number kappa (in,integer)
!   e   : energy (in,real)
!   np  : order of predictor-corrector polynomial (in,integer)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   nn  : number of nodes (out,integer)
!   g0  : m th energy derivative of the major component multiplied by r
!         (out,real(nr))
!   g1  : radial derivative of g0 (out,real(nr))
!   f0  : m th energy derivative of the minor component multiplied by r
!         (out,real(nr))
!   f1  : radial derivative of f0 (out,real(nr))
! !DESCRIPTION:
!   Finds the solution to the $m$th energy derivative of the radial Dirac
!   equation using the routine {\tt rdiracint}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m
integer, intent(in) :: kpa
real(8), intent(in) :: e
integer, intent(in) :: np
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: vr(nr)
integer, intent(out) :: nn
real(8), intent(out) :: g0(nr)
real(8), intent(out) :: g1(nr)
real(8), intent(out) :: f0(nr)
real(8), intent(out) :: f1(nr)
! local variables
integer im
! automatic arrays
real(8) g0p(nr),f0p(nr)
if (nr.le.0) then
  write(*,*)
  write(*,'("Error(rdiracdme): invalid nr : ",I8)') nr
  write(*,*)
  stop
end if
if ((m.lt.0).or.(m.gt.6)) then
  write(*,*)
  write(*,'("Error(rdiracdme): m out of range : ",I8)') m
  write(*,*)
  stop
end if
if (m.eq.0) then
  call rdiracint(m,kpa,e,np,nr,r,vr,nn,g0p,f0p,g0,g1,f0,f1)
else
  do im=0,m
    call rdiracint(im,kpa,e,np,nr,r,vr,nn,g0p,f0p,g0,g1,f0,f1)
    g0p(:)=g0(:)
    f0p(:)=f0(:)
  end do
end if
return
end subroutine
!EOC
