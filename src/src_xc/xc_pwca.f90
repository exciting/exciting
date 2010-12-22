
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: xc_pwca
! !INTERFACE:
subroutine xc_pwca(n,rhoup,rhodn,ex,ec,vxup,vxdn,vcup,vcdn)
! !INPUT/OUTPUT PARAMETERS:
!   n     : number of density points (in,integer)
!   rhoup : spin-up charge density (in,real(n))
!   rhodn : spin-down charge density (in,real(n))
!   ex    : exchange energy density (out,real(n))
!   ec    : correlation energy density (out,real(n))
!   vxup  : spin-up exchange potential (out,real(n))
!   vxdn  : spin-down exchange potential (out,real(n))
!   vcup  : spin-up correlation potential (out,real(n))
!   vcdn  : spin-down correlation potential (out,real(n))
! !DESCRIPTION:
!   Spin-polarised exchange-correlation potential and energy of the Perdew-Wang
!   parameterisation of the Ceperley-Alder electron gas: {\it Phys. Rev. B}
!   {\bf 45}, 13244 (1992) and {\it Phys. Rev. Lett.} {\bf 45}, 566 (1980).
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: rhoup(n)
real(8), intent(in) :: rhodn(n)
real(8), intent(out) :: ex(n)
real(8), intent(out) :: ec(n)
real(8), intent(out) :: vxup(n)
real(8), intent(out) :: vxdn(n)
real(8), intent(out) :: vcup(n)
real(8), intent(out) :: vcdn(n)
! local variables
integer i
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: thrd=1.d0/3.d0
real(8), parameter :: thrd2=2.d0/3.d0
real(8), parameter :: thrd4=4.d0/3.d0
real(8), parameter :: f20=1.709921d0
! beyond RPA
real(8), parameter :: p=1.d0
real(8) a(3),a1(3),b1(3),b2(3),b3(3),b4(3)
data a  / 0.0310907d0, 0.01554535d0, 0.0168869d0 /
data a1 / 0.21370d0,   0.20548d0,    0.11125d0   /
data b1 / 7.5957d0,   14.1189d0,    10.357d0     /
data b2 / 3.5876d0,    6.1977d0,     3.6231d0    /
data b3 / 1.6382d0,    3.3662d0,     0.88026d0   /
data b4 / 0.49294d0,   0.62517d0,    0.49671d0   /
real(8) rup,rdn,r,rs,srs
real(8) z,z4,ec0,ec1,ac,fz
real(8) t1,t2,t3,dzf,dzec
real(8) drec0,drec1,drac,drec
real(8) q0(3),q1(3),q1p
if (n.le.0) then
  write(*,*)
  write(*,'("Error(xc_pwca): invalid n : ",I8)') n
  write(*,*)
  stop
end if
do i=1,n
  rup=rhoup(i); rdn=rhodn(i)
! total density
  r=rup+rdn
  if ((rup.ge.0.d0).and.(rdn.ge.0.d0).and.(r.gt.1.d-12)) then
    rs=(3.d0/(4.d0*pi*r))**thrd
    srs=sqrt(rs)
    z=(rup-rdn)/r
    z4=z**4
! exchange energy density
    ex(i)=-(3.d0/(4.d0*pi*rs))*((9.d0*pi/4.d0)**thrd) &
     *((1.d0+z)**thrd4+(1.d0-z)**thrd4)/2.d0
! correlation energy density
    q0(1)=-2.d0*a(1)*(1.d0+a1(1)*rs)
    q1(1)=2.d0*a(1)*(b1(1)*srs+b2(1)*rs+b3(1)*(srs**3)+b4(1)*rs**(p+1.d0))
    ec0=q0(1)*log(1.d0+1.d0/q1(1))
    q0(2)=-2.d0*a(2)*(1.d0+a1(2)*rs)
    q1(2)=2.d0*a(2)*(b1(2)*srs+b2(2)*rs+b3(2)*(srs**3)+b4(2)*rs**(p+1.d0))
    ec1=q0(2)*log(1.d0+1.d0/q1(2))
    q0(3)=-2.d0*a(3)*(1.d0+a1(3)*rs)
    q1(3)=2.d0*a(3)*(b1(3)*srs+b2(3)*rs+b3(3)*(srs**3)+b4(3)*rs**(p+1.d0))
    ac=-q0(3)*log(1.d0+1.d0/q1(3))
    fz=((1.d0+z)**thrd4+(1.d0-z)**thrd4-2.d0)/(2.d0**thrd4-2.d0)
    ec(i)=ec0+ac*(fz/f20)*(1.d0-z4)+(ec1-ec0)*fz*z4
! exchange potential
    t1=((3.d0/(2.d0*pi))**thrd2)/(2.d0*rs)
    t2=-t1*((1.d0-z)**thrd4+(1.d0+z)**thrd4)
    t3=t1*((1.d0-z)**thrd-(1.d0+z)**thrd)
    vxup(i)=t2-(z-1.d0)*t3
    vxdn(i)=t2-(z+1.d0)*t3
! correlation potential
    dzf=thrd4*((1.d0+z)**thrd-(1.d0-z)**thrd)/(2.d0**thrd4-2.d0)
    dzec=4.d0*(z**3)*fz*(ec1-ec0-ac/f20)+dzf*((z4)*ec1-(z4)*ec0 &
     +(1.d0-z4)*ac/f20)
    q1p=a(1)*(b1(1)/srs+2.d0*b2(1)+3.d0*b3(1)*srs+2.d0*(p+1.d0)*b4(1)*rs**p)
    drec0=-2.d0*a(1)*a1(1)*log(1.d0+1.d0/q1(1))-q0(1)*q1p/(q1(1)**2+q1(1))
    q1p=a(2)*(b1(2)/srs+2.d0*b2(2)+3.d0*b3(2)*srs+2.d0*(p+1.d0)*b4(2)*rs**p)
    drec1=-2.d0*a(2)*a1(2)*log(1.d0+1.d0/q1(2))-q0(2)*q1p/(q1(2)**2+q1(2))
    q1p=a(3)*(b1(3)/srs+2.d0*b2(3)+3.d0*b3(3)*srs+2.d0*(p+1.d0)*b4(3)*rs**p)
    drac=2.d0*a(3)*a1(3)*log(1.d0+1.d0/q1(3))+q0(3)*q1p/(q1(3)**2+q1(3))
    drec=drec0*(1.d0-fz*z4)+drec1*fz*z4+drac*(fz/f20)*(1.d0-z4)
    t1=ec(i)-(rs/3.d0)*drec
    vcup(i)=t1-(z-1.d0)*dzec
    vcdn(i)=t1-(z+1.d0)*dzec
  else
    ex(i)=0.d0
    ec(i)=0.d0
    vxup(i)=0.d0
    vxdn(i)=0.d0
    vcup(i)=0.d0
    vcdn(i)=0.d0
  end if
end do
return
end subroutine
!EOC
