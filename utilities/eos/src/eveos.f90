real(8) function eveos(etype,param,v)
implicit none
! arguments
integer, intent(in) :: etype
real(8), intent(in) :: param(*)
real(8), intent(in) :: v
! local variables
real(8) v0,e0,b0,b0p,b0pp
real(8) t1,t2,t3,t4,t5,t6,t7
eveos=0.d0
select case(etype)
case(1)
! Universal equation of state
  v0=param(1)
  e0=param(2)
  b0=param(3)
  b0p=param(4)
  if (v0.lt.1.d-5) v0=1.d-5
  if (abs(b0p-1.d0).lt.1.d-5) b0p=b0p+1.d-5
  t1=b0*v0
  t2=b0p-1.d0
  t3=(v/v0)**(1.d0/3.d0)
  t4=exp(-3.d0/2.d0*t2*(-1.d0+t3))
  t5=t2**2
  t6=1.d0/t5
  eveos=-2.d0*t1*t4*(3.d0*t3*b0p-3.d0*t3+5.d0-3.d0*b0p)*t6+4.d0*t1*t6+e0
case(2)
! Murnaghan equation of state
  v0=param(1)
  e0=param(2)
  b0=param(3)
  b0p=param(4)
  if (v0.lt.1.d-5) v0=1.d-5
  if (abs(b0p).lt.1.d-5) b0p=1.d-5
  if (abs(b0p-1.d0).lt.1.d-5) b0p=b0p+1.d-5
  t1=(v0/v)**b0p
  t2=1.d0/(b0p-1.d0)
  eveos=b0*(b0p-1.d0+t1)/b0p*t2*v-b0*v0*t2+e0
case(3)
! Birch-Murnaghan third-order equation of state
  v0=param(1)
  e0=param(2)
  b0=param(3)
  b0p=param(4)
  if (v0.lt.1.d-5) v0=1.d-5
  t1=(v0/v)**(1.d0/3.d0)
  t2=t1**2
  t3=t2-1.d0
  eveos=9.d0/8.d0*b0*v0*t3**2*(b0p*t3/2.d0-2.d0*t2+3.d0)+e0
case(4)
! Birch-Murnaghan fourth-order equation of state
  v0=param(1)
  e0=param(2)
  b0=param(3)
  b0p=param(4)
  b0pp=param(5)
  if (v0.lt.1.d-5) v0=1.d-5
  t1=(v0/v)**(1.d0/3.d0)
  t2=t1**2
  t3=t2-1.d0
  t4=t3**2/4.d0
  t5=b0p**2
  eveos=3.d0/8.d0*b0*v0*t4*(9.d0*t4*b0*b0pp+9.d0*t4*t5-63.d0*t4*b0p+143.d0*t4 &
   +6.d0*b0p*t3-24.d0*t2+36.d0)+e0
case(5)
! Natural strain third-order equation of state
  v0=param(1)
  e0=param(2)
  b0=param(3)
  b0p=param(4)
  if (v0.lt.1.d-5) v0=1.d-5
  t1=b0*v0
  t2=log(v0/v)
  t3=t2**2
  t4=t3*t2
  eveos=t1*t3/2.d0+t1*t4*b0p/6.d0-t1*t4/3.d0+e0
case(6)
! Natural strain fourth-order equation of state
  v0=param(1)
  e0=param(2)
  b0=param(3)
  b0p=param(4)
  b0pp=param(5)
  if (v0.lt.1.d-5) v0=1.d-5
  t1=b0*v0
  t2=log(v0/v)
  t3=t2**2
  t4=t3**2
  t5=b0**2
  t6=b0p**2
  t7=t3*t2
  eveos=t1*t4/8.d0+t5*v0*t4*b0pp/24.d0-t1*t4*b0p/8.d0+t1*t4*t6/24.d0 &
   +t1*t7*b0p/6.d0-t1*t7/3.d0+t1*t3/2.d0+e0
case(7)
! cubic polynomial
  v0=param(1)
  e0=param(2)
  b0=param(3)
  b0p=param(4)
  if (v0.lt.1.d-5) v0=1.d-5
  t1=v0**2
  t2=v0-v
  t3=t2**2
  eveos=(1.d0+b0p)*b0/t1*t3*t2/6.d0+b0/v0*t3/2.d0+e0
case default
  write(*,*)
  write(*,'("Error(eveos): etype not defined : ",I4)') etype
  write(*,*)
  stop
end select
return
end function
