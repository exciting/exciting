
! This routine is based on code written by K. Burke.

!BOP
! !ROUTINE: xc_pbe
! !INTERFACE:
subroutine xc_pbe(n,kappa,mu,beta,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho, &
 g3up,g3dn,ex,ec,vxup,vxdn,vcup,vcdn)
! !INPUT/OUTPUT PARAMETERS:
!   n     : number of density points (in,integer)
!   kappa : parameter for large-gradient limit (in,real)
!   mu    : gradient expansion coefficient (in,real)
!   beta  : gradient expansion coefficient (in,real)
!   rhoup : spin-up charge density (in,real(n))
!   rhodn : spin-down charge density (in,real(n))
!   grho  : |grad rho| (in,real(n))
!   gup   : |grad rhoup| (in,real(n))
!   gdn   : |grad rhodn| (in,real(n))
!   g2up  : grad^2 rhoup (in,real(n))
!   g2dn  : grad^2 rhodn (in,real(n))
!   g3rho : (grad rho).(grad |grad rho|) (in,real(n))
!   g3up  : (grad rhoup).(grad |grad rhoup|) (in,real(n))
!   g3dn  : (grad rhodn).(grad |grad rhodn|) (in,real(n))
!   ex    : exchange energy density (out,real(n))
!   ec    : correlation energy density (out,real(n))
!   vxup  : spin-up exchange potential (out,real(n))
!   vxdn  : spin-down exchange potential (out,real(n))
!   vcup  : spin-up correlation potential (out,real(n))
!   vcdn  : spin-down correlation potential (out,real(n))
! !DESCRIPTION:
!   Spin-polarised exchange-correlation potential and energy of the generalised
!   gradient approximation functional of J. P. Perdew, K. Burke and M. Ernzerhof
!   {\it Phys. Rev. Lett.} {\bf 77}, 3865 (1996) and {\bf 78}, 1396(E) (1997).
!   The parameter $\kappa$, which controls the large-gradient limit, can be set
!   to $0.804$ or $1.245$ corresponding to the value in the original article or
!   the revised version of Y. Zhang and W. Yang, {\it Phys. Rev. Lett.}
!   {\bf 80}, 890 (1998).
!
! !REVISION HISTORY:
!   Modified routines written by K. Burke, October 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: kappa
real(8), intent(in) :: mu
real(8), intent(in) :: beta
real(8), intent(in) :: rhoup(n)
real(8), intent(in) :: rhodn(n)
real(8), intent(in) :: grho(n)
real(8), intent(in) :: gup(n)
real(8), intent(in) :: gdn(n)
real(8), intent(in) :: g2up(n)
real(8), intent(in) :: g2dn(n)
real(8), intent(in) :: g3rho(n)
real(8), intent(in) :: g3up(n)
real(8), intent(in) :: g3dn(n)
real(8), intent(out) :: ex(n)
real(8), intent(out) :: ec(n)
real(8), intent(out) :: vxup(n)
real(8), intent(out) :: vxdn(n)
real(8), intent(out) :: vcup(n)
real(8), intent(out) :: vcdn(n)
! local variables
integer i
real(8), parameter :: thrd=1.d0/3.d0
real(8), parameter :: thrd2=2.d0/3.d0
real(8), parameter :: pi=3.1415926535897932385d0
real(8) r,r2,kf,s,u,v
real(8) rs,z,g,ks,ksg
real(8) t,uu,vv,ww
real(8) g2rho,exup,exdn
do i=1,n
  if ((rhoup(i).gt.1.d-12).and.(rhodn(i).gt.1.d-12)) then
! exchange energy density and potential
! spin-up
    r=rhoup(i)
    r2=2.d0*r
    kf=(r2*3.d0*pi**2)**thrd
    s=gup(i)/(2.d0*kf*r)
    u=g3up(i)/((r**2)*(2.d0*kf)**3)
    v=g2up(i)/(r*(2.d0*kf)**2)
    call x_pbe(kappa,mu,r2,s,u,v,exup,vxup(i))
! spin-down
    r=rhodn(i)
    r2=2.d0*r
    kf=(r2*3.d0*pi**2)**thrd
    s=gdn(i)/(2.d0*kf*r)
    u=g3dn(i)/((r**2)*(2.d0*kf)**3)
    v=g2dn(i)/(r*(2.d0*kf)**2)
    call x_pbe(kappa,mu,r2,s,u,v,exdn,vxdn(i))
! total density
    r=rhoup(i)+rhodn(i)
! average exchange energy density
    ex(i)=(exup*rhoup(i)+exdn*rhodn(i))/r
! correlation
    rs=(3.d0/(4.d0*pi*r))**thrd
    z=(rhoup(i)-rhodn(i))/r
    g=((1.d0+z)**thrd2+(1.d0-z)**thrd2)/2.d0
    kf=(r*3.d0*pi**2)**thrd
    ks=sqrt(4.d0*kf/pi)
    ksg=2.d0*ks*g
    t=grho(i)/(ksg*r)
    uu=g3rho(i)/((r**2)*ksg**3)
    g2rho=g2up(i)+g2dn(i)
    vv=g2rho/(r*ksg**2)
    ww=(gup(i)**2-gdn(i)**2-z*grho(i)**2)/(r*r*ksg**2)
    call c_pbe(beta,rs,z,t,uu,vv,ww,ec(i),vcup(i),vcdn(i))
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

