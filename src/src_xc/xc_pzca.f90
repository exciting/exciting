
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: xc_pzca
! !INTERFACE:
subroutine xc_pzca(n,rho,ex,ec,vx,vc)
! !INPUT/OUTPUT PARAMETERS:
!   n   : number of density points (in,integer)
!   rho : charge density (in,real(n))
!   ex  : exchange energy density (out,real(n))
!   ec  : correlation energy density (out,real(n))
!   vx  : exchange potential (out,real(n))
!   vc  : correlation potential (out,real(n))
! !DESCRIPTION:
!   Spin-unpolarised exchange-correlation potential and energy of the
!   Perdew-Zunger parameterisation of Ceperley-Alder electron gas: {\it Phys.
!   Rev. B} {\bf 23}, 5048 (1981) and {\it Phys. Rev. Lett.} {\bf 45}, 566
!   (1980).
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: rho(n)
real(8), intent(out) :: ex(n)
real(8), intent(out) :: ec(n)
real(8), intent(out) :: vx(n)
real(8), intent(out) :: vc(n)
! local variables
integer i
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: thrd=1.d0/3.d0
real(8), parameter :: thrd2=2.d0/3.d0
real(8), parameter :: thrd4=4.d0/3.d0
real(8), parameter :: g=-0.1423d0,b1=1.0529d0,b2=0.3334d0
real(8), parameter :: a=0.0311d0,b=-0.048d0,c=0.0020d0,d=-0.0116d0
real(8) r,rs,srs,lrs
if (n.le.0) then
  write(*,*)
  write(*,'("Error(xc_pzca): invalid n : ",I8)') n
  write(*,*)
  stop
end if
do i=1,n
  r=rho(i)
  if (r.gt.1.d-20) then
    rs=(3.d0/(4.d0*pi*r))**thrd
! exchange energy and potential
    ex(i)=(-3.d0/(4.d0*pi))*(3.d0*r*(pi**2))**thrd
    vx(i)=thrd4*ex(i)
! correlation energy and potential
    if (rs.ge.1.d0) then
      srs=sqrt(rs)
      ec(i)=g/(1.d0+b1*srs+b2*rs)
      vc(i)=ec(i)*(1.d0+(7.d0/6.d0)*b1*srs+thrd4*b2*rs)/(1.d0+b1*srs+b2*rs)
    else
      lrs=dlog(rs)
      ec(i)=a*lrs+b+c*rs*lrs+d*rs
      vc(i)=a*lrs+(b-thrd*a)+thrd2*c*rs*lrs+thrd*(2.d0*d-c)*rs
    end if
  else
    ex(i)=0.d0
    ec(i)=0.d0
    vx(i)=0.d0
    vc(i)=0.d0
  end if
end do
return
end subroutine
!EOC
