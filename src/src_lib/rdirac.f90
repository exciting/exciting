!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rdirac
! !INTERFACE:
!
!
Subroutine rdirac (m, n, l, k, np, nr, r, vr, eval, g0, f0,dirac_eq)
! !INPUT/OUTPUT PARAMETERS:
!   m    : energy-derivative order (in,integer)
!   n    : principal quantum number (in,integer)
!   l    : quantum number l (in,integer)
!   k    : quantum number k (l or l+1) (in,integer)
!   np   : order of predictor-corrector polynomial (in,integer)
!   nr   : number of radial mesh points (in,integer)
!   r    : radial mesh (in,real(nr))
!   vr   : potential on radial mesh (in,real(nr))
!   eval : eigenvalue without rest-mass energy (inout,real)
!   g0   : major component of the radial wavefunction (out,real(nr))
!   f0   : minor component of the radial wavefunction (out,real(nr))
! !DESCRIPTION:
!   Finds the solution to the radial Dirac equation for a given potential $v(r)$
!   and quantum numbers $n$, $k$ and $l$. The method involves integrating the
!   equation using the predictor-corrector method and adjusting $E$ until the
!   number of nodes in the wavefunction equals $n-l-1$. The calling routine must
!   provide an initial estimate for the eigenvalue. Note that the arrays
!   {\tt g0} and {\tt f0} represent the radial functions multiplied by $r$.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Logical dirac_eq
      Integer, Intent (In) :: m
      Integer, Intent (In) :: n
      Integer, Intent (In) :: l
      Integer, Intent (In) :: k
      Integer, Intent (In) :: np
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: vr (nr)
      Real (8), Intent (Inout) :: eval
      Real (8), Intent (Out) :: g0 (nr)
      Real (8), Intent (Out) :: f0 (nr)
! local variables
      Integer, Parameter :: maxit = 200
      Integer :: kpa, it, nn, ir, irm, nnd, nndp,count
      real(8) ::e_hi,e_lo,f_lo,f_hi,e_guess,f_mi,e_mi
      real(8) :: xa,xb,xc,ya,yb,yc,int_c,int_b,int_a
      logical :: lo_found,hi_found
! energy convergence tolerance
      Real (8), Parameter :: eps = 1.d-11
      Real (8) :: t1, de
! automatic arrays
      Real (8) :: g1 (nr), f1 (nr), fr (nr), gr (nr), cf (3, nr)
  
      If (k .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(rdirac): k <= 0 : ", I8)') k
         Write (*,*)
         Stop
      End If
      If (k .Gt. n) Then
         Write (*,*)
         Write (*, '("Error(rdirac): incompatible n and k : ", 2I8)') &
        & n, k
         Write (*,*)
         Stop
      End If
      If ((k .Eq. n) .And. (l .Ne. k-1)) Then
         Write (*,*)
         Write (*, '("Error(rdirac): incompatible n, k and l : ", 3I8)') n, k, l
         Write (*,*)
         Stop
      End If
      If (k .Eq. l) Then
         kpa = k
      Else If (k .Eq. l+1) Then
         kpa = - k
      Else
         Write (*,*)
         Write (*, '("Error(rdirac): incompatible l and k : ", 2I8)') &
        & l, k
         Write (*,*)
         Stop
      End If
      de = 0.01d0
      nndp = 0
! count=1
      lo_found=.false.
      hi_found=.false.
!         write(*,*) 'bracketting'
!write(*,*) 'initial guess', eval
! integrate the Dirac equation
       if (dirac_eq) then
         Call rdiracdme (m, kpa, eval, np, nr, r, vr, nn, g0, g1, f0, f1)
       else
         Call rschroddme (m, l, 0, eval, np, nr, r, vr, nn, g0, g1, f0, f1)
       endif
!      if (nn.eq.n-l)   f_hi=g0(ir)

         if (nn.le.n-l-1) then
             
           do while (nn.le.n-l-1)
               e_lo=eval
               lo_found=(nn.eq.n-l-1) 
               de=de*2d0
               eval = eval + de
       if (dirac_eq) then
         Call rdiracdme (m, kpa, eval, np, nr, r, vr, nn, g0, g1, f0, f1)
       else
         Call rschroddme (m, l, 0, eval, np, nr, r, vr, nn, g0, g1, f0, f1)
       endif
!               count=count+1
           enddo
           e_hi=eval
           hi_found=(nn.eq.n-l)
           
         else
           do while (nn.gt.n-l-1)
               e_hi=eval
               hi_found=(nn.eq.n-l)
               de=de*2d0
               eval = eval - de
       if (dirac_eq) then
         Call rdiracdme (m, kpa, eval, np, nr, r, vr, nn, g0, g1, f0, f1)
       else
         Call rschroddme (m, l, 0, eval, np, nr, r, vr, nn, g0, g1, f0, f1)
       endif
!               count=count+1
           enddo
           e_lo=eval
           lo_found=(nn.eq.n-l)
         endif
         do while ((.not.lo_found).or.(.not.hi_found))
!           count=count+1
           eval=0.5d0*(e_lo+e_hi)
       if (dirac_eq) then
         Call rdiracdme (m, kpa, eval, np, nr, r, vr, nn, g0, g1, f0, f1)
       else
         Call rschroddme (m, l, 0, eval, np, nr, r, vr, nn, g0, g1, f0, f1)
       endif
           if (nn.gt.n-l-1) then
             e_hi=eval
             hi_found=(nn.eq.n-l)
           else
             e_lo=eval
             lo_found=(nn.eq.n-l-1)
           endif
         enddo
         
 
! searching for a boundary
      ir=1
      do while ((abs(g0(ir)).lt.1d50).and.(ir.lt.nr))
         ir=ir+1
!         write(*,*) ir,g0(ir),f0(ir)
      enddo
!      write(*,*) ir
!! Now we try to turn ir-th point to 0
!       stop  

if (eval.eq.e_lo) then
 f_lo=g0(ir)
       if (dirac_eq) then
         Call rdiracdme (m, kpa, e_hi, np, nr, r, vr, nn, g0, g1, f0, f1)
       else
         Call rschroddme (m, l, 0, e_hi, np, nr, r, vr, nn, g0, g1, f0, f1)
       endif
 f_hi=g0(ir)
else
 f_hi=g0(ir)
       if (dirac_eq) then
         Call rdiracdme (m, kpa, e_lo, np, nr, r, vr, nn, g0, g1, f0, f1)
       else
         Call rschroddme (m, l, 0, e_lo, np, nr, r, vr, nn, g0, g1, f0, f1)
       endif
 f_lo=g0(ir)
endif
! count=count+1
!f_lo=g0(ir)
!Call rdiracdme (0, kpa, e_hi, np, nr, r, vr, nn, g0, g1, f0, f1)
!f_hi=g0(ir)

!write(*,*) nn
!write(*,*) f_lo,f_hi      
!write(*,*)

! do it=0,100
!   eval=e_lo+dble(it)*0.01*(e_hi-e_lo)
!   Call rdiracdme (0, kpa, eval, np, nr, r, vr, nn, g0, g1, f0, f1)
!   write(*,*) g0(ir)
! enddo
! write(*,*)
! stop
      Do it = 1, maxit
! finding parameters A, B and C for an interpolating function y=A+B*exp(-C*x) 
!that passes through the three points: (e_lo,f_lo), (e_mi,f_mi) and (e_lo,f_lo)

e_mi=0.5d0*(e_hi+e_lo)
       if (dirac_eq) then
         Call rdiracdme (m, kpa, e_mi, np, nr, r, vr, nn, g0, g1, f0, f1)
       else
         Call rschroddme (m, l, 0, e_mi, np, nr, r, vr, nn, g0, g1, f0, f1)
       endif
f_mi=g0(ir)
! count=count+1
! step 1. searching for C
!      write(*,*) e_lo,f_lo
!      write(*,*) e_mi,f_mi
!      write(*,*) e_hi,f_hi
!      write(*,*) "diff",e_hi-e_mi,e_mi-e_lo

      
      if (abs(f_hi-f_mi)/(e_hi-e_mi).lt.abs(f_mi-f_lo)/(e_mi-e_lo)) then
        xa=1d-8
        ya=((f_hi-f_lo)*(e_mi-e_lo)/(e_hi-e_lo)-(f_mi-f_lo))*xa        
        xb=1d0
        yb=(f_hi-f_mi)-(f_hi-f_lo)*exp(-xb*(e_mi-e_lo)/((e_hi-e_lo)))+(f_mi-f_lo)*exp(-xb)
        do while ((yb*ya.ge.0).and.(abs((yb-ya)/ya).gt.1d-10))

           xa=xb
           ya=yb
           xb=xb*2d0
           yb=(f_hi-f_mi)-(f_hi-f_lo)*exp(-xb*(e_mi-e_lo)/((e_hi-e_lo)))+(f_mi-f_lo)*exp(-xb)
        enddo
       else
        xb=-1d-8
        yb=((f_hi-f_lo)*(e_mi-e_lo)/(e_hi-e_lo)-(f_mi-f_lo))*xb
        xa=-1d0
        ya=(f_hi-f_mi)-(f_hi-f_lo)*exp(-xa*(e_mi-e_lo)/((e_hi-e_lo)))+(f_mi-f_lo)*exp(-xa)
        
        do while ((yb*ya.ge.0).and.(abs(yb).lt.1d100))
           xb=xa
           yb=ya
           xa=xa*2d0
           ya=(f_hi-f_mi)-(f_hi-f_lo)*exp(-xa*(e_mi-e_lo)/((e_hi-e_lo)))+(f_mi-f_lo)*exp(-xa)
        enddo
       endif

      if ((abs((yb-ya)/ya).gt.1d-10).and.(abs(yb).lt.1d100)) then
        do while ((xb-xa)/xb.gt.1d-12)
           xc=0.5d0*(xa+xb)
           yc=(f_hi-f_mi)-(f_hi-f_lo)*exp(-xc*(e_mi-e_lo)/((e_hi-e_lo)))+(f_mi-f_lo)*exp(-xc)
           if (ya*yc.lt.0) then
             xb=xc
             yb=yc
           else
             xa=xc
             ya=yc
           endif
        enddo
        int_c=0.5d0*(xa+xb)
! step 2. calculating A and B
        int_b=(f_mi-f_lo)/(exp(-int_c*(e_mi-e_lo)/(e_hi-e_lo))-1d0)        
        int_a=f_lo-int_b
! Now we use y=A+B*exp(-C*x), to estimate the root/eigenvalue
        if (-int_a/int_b.gt.0d0) then
         e_guess=e_lo+(-log(-int_a/int_b)/int_c)*(e_hi-e_lo)
        else
         e_guess=e_hi+abs(e_hi)
        endif
        
       else
        if (f_mi*f_lo.lt.0d0) then
         e_guess=0.5d0*(e_mi+e_lo)
        else
         e_guess=0.5d0*(e_mi+e_hi)
        endif
       endif
!        write(*,*)
!        write(*,*) e_lo,e_hi
!        write(*,*) e_mi,e_guess
 
        if ((e_guess.gt.e_hi).or.(e_guess.lt.e_lo)) then
         if (f_mi*f_lo.lt.0d0) then
          e_guess=0.5d0*(e_mi+e_lo)
         else
          e_guess=0.5d0*(e_mi+e_hi)
         endif
        endif
         
        if (IsNaN(e_guess)) then
          write(*,*) 'watch out for NaN'
          stop
        endif
       if (dirac_eq) then
         Call rdiracdme (m, kpa, e_guess, np, nr, r, vr, nn, g0, g1, f0, f1)
       else
         Call rschroddme (m, l, 0, e_guess, np, nr, r, vr, nn, g0, g1, f0, f1)
       endif
!        count=count+1 
!        write(*,*) e_lo,e_hi,e_mi,e_guess
!        write(*,*) f_lo,f_hi,f_mi

        if (e_hi-e_guess.lt.1d-11*(abs(e_guess)+0d0)) then
           e_lo=e_guess
           Go To 20
        elseif (e_guess-e_lo.lt.1d-11*(abs(e_guess)+0d0)) then
           Go To 20
        endif
       

        if (e_guess.lt.e_mi) then
           e_hi=e_mi
           e_mi=e_guess
           f_hi=f_mi
           f_mi=g0(ir)
        else
           e_lo=e_mi
           e_mi=e_guess
           f_lo=f_mi
           f_mi=g0(ir)
        endif

        if (f_mi*f_lo.lt.0d0) then
          f_hi=f_mi
          e_hi=e_mi
        else
          f_lo=f_mi
          e_lo=e_mi
        endif


!        if (e_hi-e_lo.lt.1d-12*(abs(e_guess)+1d0)) Go To 20
!       write(*,*) e_lo,e_hi
!       write(*,*) f_lo,f_hi
!        write(*,*) ':'
!        read(*,*)

!!        if (e_hi-e_lo.lt.1d-12*(abs(e_lo)+1d0)) Go To 20       

!       write(*,*) e_lo,e_hi,e_guess

      End Do

!      stop
      Write (*,*)
      Write (*, '("Error(rdirac): maximum iterations exceeded")')
      Write (*,*)
      Stop
20    Continue
!      stop
      eval=e_lo
!      write(*,*) count
!      write(*,*) e_lo,it
!      write(*,*) 'it=',it
!      if (n.eq.1) then
!        write(*,*) ':'
!        read(*,*)
!      endif
! find effective infinity and set wavefunction to zero after that point
! major component
      irm = nr
      Do ir = 2, nr
         If ((g0(ir-1)*g0(ir) .Lt. 0.d0) .Or. (g1(ir-1)*g1(ir) .Lt. &
        & 0.d0)) irm = ir
      End Do
      g0 (irm:nr) = 0.d0
! minor component
      irm = nr
      Do ir = 2, nr
         If ((f0(ir-1)*f0(ir) .Lt. 0.d0) .Or. (f1(ir-1)*f1(ir) .Lt. &
        & 0.d0)) irm = ir
      End Do
      f0 (irm:nr) = 0.d0
! normalise
      if (dirac_eq) then
        Do ir = 1, nr
          fr (ir) = g0 (ir) ** 2 + f0 (ir) ** 2
        End Do
      else
        Do ir = 1, nr
          fr (ir) = g0 (ir) ** 2 
        End Do
      endif
      Call fderiv (-1, nr, r, fr, gr, cf)
      t1 = Sqrt (Abs(gr(nr)))
      If (t1 .Gt. 0.d0) Then
         t1 = 1.d0 / t1
      Else
         Write (*,*)
         Write (*, '("Error(rdirac): zero wavefunction")')
         Write (*,*)
         Stop
      End If
      g0 (:) = t1 * g0 (:)
      f0 (:) = t1 * f0 (:)
      Return
End Subroutine
!EOC
