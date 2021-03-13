! Copyright (C) 2007-2010 D. Nabok, P. Puschnig and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!============================================================
subroutine vdWDF_VV(dis,dens1,graddens1,dens2,graddens2,func)
!-----------------------------------------
use param
!-----------------------------------------

implicit none
real*8  :: dis, func
real*8  :: dens1, dens2
real*8  :: graddens1, graddens2

!-----------------------------------------
real*8  :: kf1, kf2
real*8  :: omegapsq1,omegapsq2,omegagsq1,omegagsq2
real*8  :: omega01, omega02
real*8  :: Aa,Bb,Dd,Kappa,kappa1,kappa2
real*8  :: dcut=1.0d-06

!-----------------------------------------   
!-----------------------------------------   

   func = 0.0d0
  
   kf1 = (3.0d0*pi*pi*dens1)**(1.0d0/3.0d0)
   kf2 = (3.0d0*pi*pi*dens2)**(1.0d0/3.0d0)
   
   omegaPsq1=4.0d0*pi*dens1
   omegaPsq2=4.0d0*pi*dens2
   
   omegagsq1=Cfac*(graddens1/dens1)**4.0d0
   omegagsq2=Cfac*(graddens2/dens2)**4.0d0
   
   omega01=dsqrt(omegagsq1+omegapsq1/3.0d0)
   omega02=dsqrt(omegagsq2+omegapsq2/3.0d0)
   
   if (dis > dcut) then
       
       kappa1=(4.0d0*kf1/pi)
       kappa2=(4.0d0*kf2/pi)
       Kappa=dis/2.0d0*dsqrt(kappa1*kappa2/(kappa1+kappa2))
       
       Aa=2.0d0*Kappa/dsqrt(pi)*dexp(-Kappa*Kappa)
       Bb=derf(Kappa)-Aa
       Dd=4.0d0/3.0d0*Kappa*Kappa*Aa*Bb-Bb*Bb
       
       func = const*3.0d0/64.0d0/pi/pi*omegaPsq1*omegaPsq2*Dd*dis**(-6.0d0)/ &
              (omega01*omega02*(omega01+omega02))
   else
   
       ! Use assymptotic formula (Vydrov2009, formula (17))
       kappa1=(4.0d0*kf1/pi)          
       func = const*3.0d0/64.0d0/pi/pi*omegaPsq1*omegaPsq2*kappa1**3.0d0/(288.0d0*pi) / &
              (omega01*omega02*(omega01+omega02))
   
   end if

return
end
