!
!
!
! Copyright (C) 2013 exciting team
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rdiracint
! !INTERFACE:
!
!
Subroutine rdiracint (m, kpa, e, nr, r, vr, nn, g0p, f0p, g0, g1, f0, f1,sloppy)
! !INPUT/OUTPUT PARAMETERS:
!   m   : order of energy derivative (in,integer)
!   kpa : quantum number kappa (in,integer)
!   e   : energy (in,real)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   nn  : number of nodes (out,integer)
!   g0p : m-1 th energy derivative of the major component multiplied by r
!         (in,real(nr))
!   f0p : m-1 th energy derivative of the minor component multiplied by r
!         (in,real(nr))
!   g0  : m th energy derivative of the major component multiplied by r
!         (out,real(nr))
!   g1  : radial derivative of g0 (out,real(nr))
!   f0  : m th energy derivative of the minor component multiplied by r
!         (out,real(nr))
!   f1  : radial derivative of f0 (out,real(nr))
!   sloppy : whether a quick-and-dirty integration algorithm should be used.
! !DESCRIPTION:
!   Integrates the $m$th energy derivative of the radial Dirac equation from
!   $r=0$ outwards. 
!
! !REVISION HISTORY:
!   Created February 2013 (Andris)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: m
      Integer, Intent (In) :: kpa
      Real (8), Intent (In) :: e
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: vr (nr)
      Integer, Intent (Out) :: nn
      Real (8), Intent (In) :: g0p (nr)
      Real (8), Intent (In) :: f0p (nr)
      Real (8), Intent (Out) :: g0 (nr)
      Real (8), Intent (Out) :: g1 (nr)
      Real (8), Intent (Out) :: f0 (nr)
      Real (8), Intent (Out) :: f1 (nr)
      Logical, Intent (In) :: sloppy
! local variables
      Integer :: ir, ir0, iter,itmax,step
      parameter (itmax=32)
      real (8) gest(itmax),fest(itmax),errorf,error,errorbound
      parameter (errorbound=1d-10)
! fine-structure constant
      Real (8), Parameter :: alpha = 1d0 / 137.03599911d0
! rest mass of electron
      Real (8), Parameter :: e0 = 1.d0 / (alpha**2) 
      Real (8) :: rkpa, ri,A,B1,B2,CC,t1,t2,t3,det,detg,detf,rmult,rimult
      real (8) :: cf(4,nr),rvr(nr),r1,r2,vr2,f0old,g0old,f1old,g1old,logBA
      real (8) :: g0pcf(4,nr),f0pcf(4,nr),g00p,f00p,mc
! automatic arrays
!      Real (8) :: ttt(itmax,itmax)

     do iter=1,nr
        rvr(iter)=vr(iter)*r(iter)
     enddo
     call spline4(nr,r,1,rvr,cf)
     if (m.ne.0) then
       mc=dble(m)*alpha
       call spline4(nr,r,1,g0p,g0pcf)
       call spline4(nr,r,1,f0p,f0pcf)
     endif

      If (nr .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(rdiracint): invalid nr : ", I8)') nr
         Write (*,*)
         Stop
      End If
      If ((m .Lt. 0) .Or. (m .Gt. 6)) Then
         Write (*,*)
         Write (*, '("Error(rdiracint): m out of range : ", I8)') m
         Write (*,*)
         Stop
      End If
      rkpa = dble (kpa)
! estimate r -> 0 boundary values
!!      f0 (1) = 0.d0
!!      g1 (1) = 1.d0
!+ alpha * dble (m) * f0p (1)
!!      f1 (1) = 0
!1d-4
      ri = 1.d0 / r (1)
      A=rkpa*ri
      B1=alpha*(e+2.d0*e0-vr(1))
      B2=alpha*(e-vr(1))
!      det =A*A-B1*B2
!      detg=-g1(1)*A +f1(1)*B1
!      detf=-g1(1)*B2+f1(1)*A
!      g0(1)=detg/det
!      f0(1)=detf/det

!      write(*,*) vr(1)*r(1)
!      f0(1)=1d0/(10d0*alpha)*(kpa+sqrt(kpa*kpa-(10d0*alpha)**2))
!      t1=sqrt(kpa*kpa-(vr(1)*r(1)*alpha)**2)
!      if (m.ne.0) then
        g0(1)=1d0
        f0(1)=1d0/(-vr(1)*r(1)*alpha)*(kpa+sqrt(kpa*kpa-(vr(1)*r(1)*alpha)**2))
        g1(1) =  B1*f0(1) - A*g0(1)
        f1(1) = -B2*g0(1) + A*f0(1)
!      endif
        if (m.ne.0) then
          g1(1) =g1(1)+mc*f0p(1)
          f1(1) =f1(1)-mc*g0p(1)
        endif
!        write(*,*) r(1),g0(1),f0(1)

!      If (m .Ne. 0) Then
!         g1 (1) = g1 (1) + alpha * dble (m) * f0p (1)
!         f1 (1) = f1 (1) - alpha * dble (m) * g0p (1)
!      End If

      nn = 0
      error=1d0
      Do ir = 2, nr
         ri = 1.d0 / r (ir)
! predictor-corrector order
         logBA=log(r(ir)/r(ir-1))

    if (.not.sloppy) then 
!Plan A
      gest(:)=1d100
      fest(:)=1d100
      step=1
      error=1d0
      gest(itmax)=0d0
      do while ((step.le.itmax).and.(error.gt.errorbound))
         g0old=g0(ir-1)
         f0old=f0(ir-1)
         g1old =g1(ir-1)
         f1old =f1(ir-1)
         r2=r(ir-1)  
         ri=1d0/r2
         rmult=exp(logBA/dble(step*2))
         rimult=1d0/rmult
         do iter=1,step*2
            r1=r2
            r2=r2*rmult
!            ri=1d0/r2
            ri=ri*rimult
            vr2=rvr(ir-1)+(r2-r(ir-1))*(cf(1,ir-1)+(r2-r(ir-1))*(cf(2,ir-1)+(r2-r(ir-1))*cf(3,ir-1)))
            vr2=vr2*ri
!*ri
            A=-(r2-r1)*0.5d0
            t1=g0old-g1old*A
! *(r2-r1)*0.5d0
            t2=f0old-f1old*A
          if (m.ne.0) then
            g00p=g0p(ir-1)+(r2-r(ir-1))*(g0pcf(1,ir-1)+(r2-r(ir-1))*(g0pcf(2,ir-1)+(r2-r(ir-1))*g0pcf(3,ir-1)))       
            f00p=f0p(ir-1)+(r2-r(ir-1))*(f0pcf(1,ir-1)+(r2-r(ir-1))*(f0pcf(2,ir-1)+(r2-r(ir-1))*f0pcf(3,ir-1)))
            t1=t1-A*f00p*mc
            t2=t2+A*g00p*mc
          endif

          B1=-alpha*(e+2.d0*e0-vr2)
          B2=alpha*(e-vr2)
          CC=rkpa*ri

          det=1d0/(1d0-A*A*(CC*CC+B1*B2))
          detg=t1*(1d0+A*CC)+t2*A*B1
          detf=t1*A*B2+t2*(1d0-A*CC)
          g0old=detg*det
          f0old=detf*det
          g1old = -B1*f0old - CC*g0old
          f1old = -B2*g0old + CC*f0old
          if (m.ne.0) then
            g1old =g1old+mc*f00p
            f1old =f1old-mc*g00p
          endif

         enddo
         gest(itmax-step+1)=g0old
         fest(itmax-step+1)=f0old
         do iter=1,step-1
           t1=1d0/((dble(step)/dble(step-iter))**2-1d0)
!           t1=ttt(step,iter)
           gest(itmax-step+1+iter)=gest(itmax-step+iter)+(gest(itmax-step+iter)-gest(itmax-step+1+iter))*t1
           fest(itmax-step+1+iter)=fest(itmax-step+iter)+(fest(itmax-step+iter)-fest(itmax-step+1+iter))*t1
         enddo
         
         error=abs((gest(itmax)-gest(itmax-1))/gest(itmax))

         step=step+1
        enddo

         g0(ir)=gest(itmax)
         f0(ir)=fest(itmax)
         g1(ir) = -B1*f0(ir) - CC*g0(ir)
         f1(ir) = -B2*g0(ir) + CC*f0(ir)
         if (m.ne.0) then
            g1(ir) =g1(ir)+mc*f0p(ir)
            f1(ir) =f1(ir)-mc*g0p(ir)
         endif
    endif         
!Plan B
     if ((error.gt.errorbound).or.sloppy) then
!      if (error.gt.errorbound) then
!       write(*,*) ir
!      if (.true.) then
       if (sloppy) then
         step=1 
       else
         step=512
       endif

         g0old=g0(ir-1)
         f0old=f0(ir-1)
         g1old =g1(ir-1)
         f1old =f1(ir-1)
         r2=r(ir-1)
         
         rmult=exp(logBA/dble(step*2))
         do iter=1,step*2
            r1=r2
            r2=r2*rmult
            ri=1d0/r2
            vr2=rvr(ir-1)+(r2-r(ir-1))*(cf(1,ir-1)+(r2-r(ir-1))*(cf(2,ir-1)+(r2-r(ir-1))*cf(3,ir-1)))
            vr2=vr2*ri
            A=-(r2-r1)*0.5d0
            t1=g0old-g1old*A
            t2=f0old-f1old*A
          if (m.ne.0) then
            g00p=g0p(ir-1)+(r2-r(ir-1))*(g0pcf(1,ir-1)+(r2-r(ir-1))*(g0pcf(2,ir-1)+(r2-r(ir-1))*g0pcf(3,ir-1)))
            f00p=f0p(ir-1)+(r2-r(ir-1))*(f0pcf(1,ir-1)+(r2-r(ir-1))*(f0pcf(2,ir-1)+(r2-r(ir-1))*f0pcf(3,ir-1)))
            t1=t1-A*f00p*mc
            t2=t2+A*g00p*mc
          endif


          B1=-alpha*(e+2.d0*e0-vr2)
          B2=alpha*(e-vr2)
          CC=rkpa*ri

          det=1d0/(1d0-A*A*(CC*CC+B1*B2))
          detg=t1*(1d0+A*CC)+t2*A*B1
          detf=t1*A*B2+t2*(1d0-A*CC)
          g0old=detg*det
          f0old=detf*det
          g1old = -B1*f0old - CC*g0old
          f1old = -B2*g0old + CC*f0old
          if (m.ne.0) then
            g1old =g1old+mc*f00p
            f1old =f1old-mc*g00p
          endif
         enddo
         g0(ir)=g0old
         f0(ir)=f0old
         g1(ir)=g1old
         f1(ir)=f1old

     endif
!     write(*,*) r(ir), g0(ir) 

!, f0(ir)
!        write(*,*) g0old,f0old
!        stop
! check for overflow
         If (Abs(g0(ir)) .Gt. 1.d100) Then
!            if ((e.lt.-4d2).and.(e.gt.-5d2)) write(*,*) ir,e
            g0 (ir:nr) = g0 (ir)
            g1 (ir:nr) = g1 (ir)
            f0 (ir:nr) = f0 (ir)
            f1 (ir:nr) = f1 (ir)
!    call timesec(ts1)
!    write(*,*) ts1-ts0

!            read(*,*)
!            stop
            Return
         End If
! check for node
         If (g0(ir-1)*g0(ir) .Lt. 0.d0) nn = nn + 1
      End Do
!    call timesec(ts1)
!    write(*,*) ts1-ts0

!      stop
!      read(*,*)
      Return
End Subroutine
!EOC
