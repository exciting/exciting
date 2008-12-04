
! This routine is based on code written by K. Burke.

subroutine c_pbe(beta,rs,z,t,uu,vv,ww,ec,vcup,vcdn)
implicit none
! arguments
real(8), intent(in) :: beta
real(8), intent(in) :: rs
real(8), intent(in) :: z
real(8), intent(in) :: t
real(8), intent(in) :: uu
real(8), intent(in) :: vv
real(8), intent(in) :: ww
real(8), intent(out) :: ec
real(8), intent(out) :: vcup
real(8), intent(out) :: vcdn
! local variables
real(8), parameter :: thrd=1.d0/3.d0
real(8), parameter :: thrdm=-thrd
real(8), parameter :: thrd2=2.d0*thrd
real(8), parameter :: thrd4=4.d0*thrd
real(8), parameter :: sixthm=thrdm/2.d0
real(8), parameter :: gam=0.5198420997897463295d0
real(8), parameter :: fzz=8.d0/(9.d0*gam)
real(8), parameter :: gamma=0.0310906908696548950d0
real(8), parameter :: eta=1.d-12
real(8) rtrs,eu,eurs,ep,eprs,alfm,alfrsm,z4,f
real(8) ecrs,fz,ecz,comm,g,g3,pon,b,b2,t2,t4
real(8) q4,q5,g4,t6,rsthrd,gz,fac,bg,bec,q8,q9
real(8) hb,hrs,fact0,fact1,hbt,hrst,hz,ht,hzt
real(8) fact2,fact3,htt,pref,fact5,h,dvcup,dvcdn
real(8) delt
delt=beta/gamma
rtrs=sqrt(rs)
call c_pbe_gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0,0.49294d0, &
 rtrs,eu,eurs)
call c_pbe_gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0,0.62517d0, &
 rtrs,ep,eprs)
call c_pbe_gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,0.49671d0, &
 rtrs,alfm,alfrsm)
z4=z**4
f=((1.d0+z)**thrd4+(1.d0-z)**thrd4-2.d0)/gam
! local contribution to correlation energy density
ec=eu*(1.d0-f*z4)+ep*f*z4-alfm*f*(1.d0-z4)/fzz
ecrs=eurs*(1.d0-f*z4)+eprs*f*z4-alfrsm*f*(1.d0-z4)/fzz
fz=thrd4*((1.d0+z)**thrd-(1.d0-z)**thrd)/gam
ecz=4.d0*(z**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu-(1.d0-z4)*alfm/fzz)
comm=ec-rs*ecrs/3.d0-z*ecz
! local contribution to correlation potential
vcup=comm+ecz
vcdn=comm-ecz
g=((1.d0+z)**thrd2+(1.d0-z)**thrd2)/2.d0
g3=g**3
pon=-ec/(g3*gamma)
b=delt/(exp(pon)-1.d0)
b2=b*b
t2=t*t
t4=t2*t2
q4=1.d0+b*t2
q5=1.d0+b*t2+b2*t4
! gradient correction to energy density
h=g3*(beta/delt)*log(1.d0+delt*q4*t2/q5)
g4=g3*g
t6=t4*t2
rsthrd=rs/3.d0
gz=(((1.d0+z)**2+eta)**sixthm-((1.d0-z)**2+eta)**sixthm)/3.d0
fac=delt/b+1.d0
bg=-3.d0*b2*ec*fac/(beta*g4)
bec=b2*fac/(beta*g3)
q8=q5*q5+delt*q4*q5*t2
q9=1.d0+2.d0*b*t2
hb=-beta*g3*b*t6*(2.d0+b*t2)/q8
hrs=-rsthrd*hb*bec*ecrs
fact0=2.d0*delt-6.d0*b
fact1=q5*q9+q4*q9*q9
hbt=2.d0*beta*g3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
hrst=rsthrd*t2*hbt*bec*ecrs
hz=3.d0*gz*h/g+hb*(bg*gz+bec*ecz)
ht=2.d0*beta*g3*q9/q8
hzt=3.d0*gz*ht/g+hbt*(bg*gz+bec*ecz)
fact2=q4*q5+b*t2*(q4*q9+q5)
fact3=2.d0*b*q5*q9+delt*fact2
htt=4.d0*beta*g3*t*(2.d0*b/q8-(q9*fact3/q8)/q8)
comm=h+hrs+hrst+t2*ht/6.d0+7.d0*t2*t*htt/6.d0
pref=hz-gz*t2*ht/g
fact5=gz*(2.d0*ht+t*htt)/g
comm=comm-pref*z-uu*htt-vv*ht-ww*(hzt-fact5)
! gradient correction to potential
dvcup=comm+pref
dvcdn=comm-pref
! add gradient corrections
ec=ec+h
vcup=vcup+dvcup
vcdn=vcdn+dvcdn
return
end subroutine

