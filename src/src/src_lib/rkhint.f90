!
!
!
! Copyright (C) 2013 exciting team
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rkhint
! !INTERFACE:
!
!
Subroutine rkhint (m, l, e, nr, r, vr, nn, p0p, q0p, p0, p1, q0, q1)
! !INPUT/OUTPUT PARAMETERS:
!   m   : order of energy derivative (in,integer)
!   l   : angular momentum quantum number (in,integer)
!   e   : energy (in,real)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   nn  : number of nodes (out,integer)
!   p0p : m-1 th energy derivative of P (in,real(nr))
!   p0  : m th energy derivative of P (out,real(nr))
!   p1  : radial derivative of p0 (out,real(nr))
!   q0  : m th energy derivative of Q (out,real(nr))
!   q1  : radial derivative of q0 (out,real(nr))
! !DESCRIPTION:
!   Integrates the $m$th energy derivative of the scalar relativistic radial
!   Schr\"{o}dinger equation from $r=0$ outwards.  
!
! !REVISION HISTORY:
!   Created July 2013 (Andris)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: m
      Integer, Intent (In) :: l
      Real (8), Intent (In) :: e
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: vr (nr)
      Integer, Intent (Out) :: nn
      Real (8), Intent (In) :: p0p (nr)
      Real (8), Intent (In) :: q0p (nr)
      Real (8), Intent (Out) :: p0 (nr)
      Real (8), Intent (Out) :: p1 (nr)
      Real (8), Intent (Out) :: q0 (nr)
      Real (8), Intent (Out) :: q1 (nr)
! local variables
      Integer :: ir, ir0, iter,itmax,step
      parameter (itmax=32)
      real (8) pest(itmax),qest(itmax),errorq,error,errorbound
      parameter (errorbound=1d-10)
! fine-structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
      Real (8) :: rm, ri, tmp1, tmp2,t1,t2
! automatic arrays
! temporary stuff
      Real (8) :: A,B1,B2,CC,t3,det,detp,detq,rmult
      real (8) :: cf(4,nr),rvr(nr),r1,r2,vr2,p0old,q0old,p1old,q1old,logBA
      real (8) :: pcf(4,nr),qcf(4,nr),p00p,q00p,ttt(itmax,itmax)
! estimate r -> 0 boundary values
       p00p=0d0
       q00p=0d0
!      itmax=32
      
!      write(*,*) 
!      write(*,*) nr,nr
!      read(*,*) 
!      write(*,*) nr
!     do ir=1,itmax
!      do iter=1,itmax
!       ttt(ir,iter)=1d0/((dble(ir)/dble(ir-iter))**2-1d0)
!      enddo
!     enddo

      do iter=1,nr
!        if (abs(vr(1)*r(1)+7d0).lt.1d-1) then
!         write(*,*) r(iter),vr(iter),vr(iter)*r(iter)
!        endif
!        rvr(iter)=vr(iter)*(r(iter)**2)
         rvr(iter)=vr(iter)*r(iter)
!         rvr(iter)=vr(iter)
      enddo
      
      call spline4(nr,r,1,rvr,cf)
      If (m .Ne. 0) Then
         call spline4(nr,r,1,p0p,pcf)
         call spline4(nr,r,1,q0p,qcf)
      endif



      p1 (1) = 1.d0
      rm = 1.d0 + 0.5d0 * (alpha**2) * (e - vr(1))
      If (m .Ne. 0) Then      
! m=1
       q1 (1) = - p0p (1)*(1d0+(alpha**2)*dble (l*(l+1))) / (4.d0*(rm*r(1))**2)
      else
       q1 (1) = 0d0
      endif
!      rm = 1.d0 + 0.5d0 * (alpha**2) * (e - vr(1))
      t1 = dble (l*(l+1)) / (2.d0*rm*r(1)**2)

      A=1d0/r(1)
      B1=2.d0*rm
      B2=-(t1+vr(1)-e)
      det =A*A-B1*B2
      detp= p1(1)*A +q1(1)*B1
      detq=-p1(1)*B2-q1(1)*A
       p0 (1) =detp/det
       q0 (1) =detq/det
!      write(*,*) detp/det,detq/det


!      If (m .Ne. 0) Then
!         q1 (1) = q1 (1) - dble (m) * p0p (1)
!      End If
      nn = 0
      Do ir = 2, nr
!         rm = 1.d0 + 0.5d0 * (alpha**2) * (energyref - vr (ir))*rmfactor
!         ri = 1.d0 / r (ir)
!         t1 = dble (l*(l+1)) / (2.d0*rm*r(ir)**2)
!         t2 = t1 + vr (ir) - e
! predictor-corrector order


      logBA=log(r(ir)/r(ir-1))
      pest(:)=1d100
      qest(:)=1d100
      step=1
      error=1d0
      do while ((step.le.itmax).and.(error.gt.errorbound))
         p0old=p0(ir-1)
         q0old=q0(ir-1)
         p1old =p1(ir-1)
         q1old =q1(ir-1)

         r2=r(ir-1)
         rmult=exp(logBA/dble(step*2))
         do iter=1,step*2
            r1=r2
            r2=r2*rmult
            vr2=rvr(ir-1)+(r2-r(ir-1))*(cf(1,ir-1)+(r2-r(ir-1))*(cf(2,ir-1)+(r2-r(ir-1))*cf(3,ir-1)))
!            vr2=vr2/(r2*r2) 
             vr2=vr2/r2
!            vr2=(vr(ir)*r(ir)*(r2-r(ir-1))+vr(ir-1)*r(ir-1)*(r(ir)-r2))/(r2*(r(ir)-r(ir-1)))
            rm = 1.d0 + 0.5d0 * (alpha**2) * (e - vr2) 
            t1 = dble (l*(l+1)) / (2.d0*rm*r2**2)
            t2 = t1 + vr2 - e
            A=-(r2-r1)*0.5d0
            tmp1=p0old+p1old*(r2-r1)*0.5d0
            tmp2=q0old+q1old*(r2-r1)*0.5d0
            If (m .Ne. 0) Then
              p00p=p0p(ir-1)+(r2-r(ir-1))*(pcf(1,ir-1)+(r2-r(ir-1))*(pcf(2,ir-1)+(r2-r(ir-1))*pcf(3,ir-1)))
              q00p=q0p(ir-1)+(r2-r(ir-1))*(qcf(1,ir-1)+(r2-r(ir-1))*(qcf(2,ir-1)+(r2-r(ir-1))*qcf(3,ir-1)))
! Assuming m=1
              tmp1=tmp1 - A*q00p*(alpha**2)
              tmp2=tmp2 + A*p00p*(1d0+(alpha**2)*dble (l*(l+1)) / (4.d0*(rm*r2)**2))
              
            End If

          B1=-2.d0*rm
          B2=-t2

          CC=-1d0/r2

          det=1d0-A*A*(CC*CC+B1*B2)
          detp=tmp1*(1d0+A*CC)+tmp2*A*B1
          detq=tmp1*A*B2+tmp2*(1d0-A*CC)
          p0old=detp/det
          q0old=detq/det
          p1old = -B1*q0old - CC*p0old+ q00p*(alpha**2)
!m=1
          q1old = -B2*p0old + CC*q0old- p00p*(1d0+(alpha**2)*dble (l*(l+1)) / (4.d0*(rm*r2)**2))
         enddo
         pest(itmax-step+1)=p0old
         qest(itmax-step+1)=q0old
         do iter=1,step-1
           t1=1d0/((dble(step)/dble(step-iter))**2-1d0)
!           t1=ttt(step,iter)
           pest(itmax-step+1+iter)=pest(itmax-step+iter)+(pest(itmax-step+iter)-pest(itmax-step+1+iter))*t1
           qest(itmax-step+1+iter)=qest(itmax-step+iter)+(qest(itmax-step+iter)-qest(itmax-step+1+iter))*t1
         enddo

         error =abs((pest(itmax)-pest(itmax-1))/pest(itmax))
!         write(*,*) step,pest(itmax)
         step=step+1
        enddo
        rm = 1.d0 + 0.5d0 * (alpha**2) * (e - vr(ir))
        p0(ir)=pest(itmax)
        q0(ir)=qest(itmax)
        p1(ir) = -B1*q0(ir) - CC*p0(ir)+ q00p*(alpha**2)
        q1(ir) = -B2*p0(ir) + CC*q0(ir)- p00p*(1d0+(alpha**2)*dble (l*(l+1)) / (4.d0*(rm*r(ir))**2))
!          p1old = -B1*q0old - CC*p0old
!          q1old = -B2*p0old + CC*q0old-dble (m) * p00p


!         g0(ir)=gest(itmax)
!         f0(ir)=fest(itmax)
!         g1(ir) = -B1*f0(ir) - CC*g0(ir)
!         f1(ir) = -B2*g0(ir) + CC*f0(ir)

!         stop
         
!plan B
       if (error.gt.errorbound) then
!         write(*,*) 'not converged',vr(1)*r(1),ir
!         if (.true.) then
         step=512
         p0old=p0(ir-1)
         q0old=q0(ir-1)
         p1old =p1(ir-1)
         q1old =q1(ir-1)

         r2=r(ir-1)
         rmult=exp(logBA/dble(step))

         do iter=1,step
            r1=r2
            r2=r2*rmult

!            r1=dble(iter-1)/dble(itmax)*(r(ir)-r(ir-1))+r(ir-1)
!            r2=dble(iter)/dble(itmax)*(r(ir)-r(ir-1))+r(ir-1)
            vr2=rvr(ir-1)+(r2-r(ir-1))*(cf(1,ir-1)+(r2-r(ir-1))*(cf(2,ir-1)+(r2-r(ir-1))*cf(3,ir-1)))
            vr2=vr2/r2
!            vr2=vr2/(r2*r2)
!            vr2=(vr(ir)*r(ir)*(r2-r(ir-1))+vr(ir-1)*r(ir-1)*(r(ir)-r2))/(r2*(r(ir)-r(ir-1)))
            rm = 1.d0 + 0.5d0 * (alpha**2) * (e-vr2)
            t1 = dble (l*(l+1)) / (2.d0*rm*r2**2)
            t2 = t1 + vr2 - e

            A=-(r2-r1)*0.5d0
            tmp1=p0old+p1old*(r2-r1)*0.5d0
            tmp2=q0old+q1old*(r2-r1)*0.5d0
            If (m .Ne. 0) Then
              p00p=p0p(ir-1)+(r2-r(ir-1))*(pcf(1,ir-1)+(r2-r(ir-1))*(pcf(2,ir-1)+(r2-r(ir-1))*pcf(3,ir-1)))
              q00p=q0p(ir-1)+(r2-r(ir-1))*(qcf(1,ir-1)+(r2-r(ir-1))*(qcf(2,ir-1)+(r2-r(ir-1))*qcf(3,ir-1)))
! Assuming m=1
              tmp1=tmp1 - A*q00p*(alpha**2)
              tmp2=tmp2 + A*p00p*(1d0+(alpha**2)*dble (l*(l+1)) / (4.d0*(rm*r2)**2))
            End If

          B1=-2.d0*rm
          B2=-t2

          CC=-1d0/r2

          det=1d0-A*A*(CC*CC+B1*B2)
          detp=tmp1*(1d0+A*CC)+tmp2*A*B1
          detq=tmp1*A*B2+tmp2*(1d0-A*CC)
          p0old=detp/det
          q0old=detq/det
          p1old = -B1*q0old - CC*p0old+q00p*(alpha**2)
          q1old = -B2*p0old + CC*q0old-p00p*(1d0+(alpha**2)*dble (l*(l+1)) / (4.d0*(rm*r2)**2))

         enddo
!         write(*,*)
!         write(*,*) p0old
!         stop
!!         if (m.eq.1) then
!           write(*,*) r(ir),q1(ir),q1old
!           write(*,*) r(ir),p1(ir),p1old
!           stop
!!         endif
         p0(ir)=p0old
         q0(ir)=q0old
         p1(ir) = p1old
         q1(ir) = q1old
       endif
!         If (m .Ne. 0) Then
!            q1 (ir) = q1 (ir) - dble (m) * p0p (ir)
!         End If
!         write(*,*) p0(ir),q0(ir)
!         stop
         

! check for overflow
         If (Abs(p0(ir)) .Gt. 1.d100) Then
            p0 (ir:nr) = p0 (ir)
            p1 (ir:nr) = p1 (ir)
            q0 (ir:nr) = q0 (ir)
            q1 (ir:nr) = q1 (ir)
!            write(*,*) 'rschrodint:',e
            Return
         End If
! check for node
         If (p0(ir-1)*p0(ir) .Lt. 0.d0) nn = nn + 1
      End Do
!      write(*,*) 'rschrodint:',e
      Return
End Subroutine
!EOC
