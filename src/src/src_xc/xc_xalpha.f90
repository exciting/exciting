
! Copyright (C) 1998-2006 ABINIT group (DCA, XG, GMR).
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: xc_xalpha
! !INTERFACE:
subroutine xc_xalpha(n,rho,exc,vxc)
! !INPUT/OUTPUT PARAMETERS:
!   n   : number of density points (in,integer)
!   rho : charge density (in,real(n))
!   exc : exchange-correlation energy density (out,real(n))
!   vxc : exchange-correlation potential (out,real(n))
! !DESCRIPTION:
!   $X_{\alpha}$ approximation to the exchange-correlation potential and energy
!   density. See J. C. Slater, {\it Phys. Rev.} {\bf 81}, 385 (1951).
!
! !REVISION HISTORY:
!   Modified an ABINIT routine, September 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: rho(n)
real(8), intent(out) :: exc(n)
real(8), intent(out) :: vxc(n)
! local variables
integer i
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: alpha=1.d0
real(8) r,efac,rs,rsm1,vfac
vfac=(1.5d0/pi)**(2.d0/3.d0)
efac=0.75d0*vfac
! loop over density points
do i=1,n
  r=rho(i)
  if (r.gt.1.d-20) then
    rs=(3.d0/(4.d0*pi*r))**(1.d0/3.d0)
    rsm1=1.0d0/rs
! compute energy density
    exc(i)=-alpha*efac*rsm1
! compute potential
    vxc(i)=-alpha*vfac*rsm1
  end if
end do
end subroutine
!EOC

