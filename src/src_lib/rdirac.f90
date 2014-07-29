!
!
!
! Copyright (C) 2013 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rdirac
! !INTERFACE:
!
!
Subroutine rdirac (m, n, l, k, nr, r, vr, eval, g0, f0,dirac_eq,sloppy)
use modinput
! !INPUT/OUTPUT PARAMETERS:
!   m    : energy-derivative order (in,integer)
!   n    : principal quantum number (in,integer)
!   l    : quantum number l (in,integer)
!   k    : quantum number k (l or l+1) (in,integer)
!   nr   : number of radial mesh points (in,integer)
!   r    : radial mesh (in,real(nr))
!   vr   : potential on radial mesh (in,real(nr))
!   eval : eigenvalue without rest-mass energy (inout,real)
!   g0   : major component of the radial wavefunction (out,real(nr))
!   f0   : minor component of the radial wavefunction (out,real(nr))
!   dirac_eq : flag to pick the equation (Dirac or Schroedinger) 
!   sloppy   : flag to pick a quick-and-dirty algorithm for integrating the Dirac equation
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
!   Rewritten 2013 (Andris)
!    
!EOP
!BOC
      Implicit None
! arguments
      Logical dirac_eq,sloppy
      Integer, Intent (In) :: m
      Integer, Intent (In) :: n
      Integer, Intent (In) :: l
      Integer, Intent (In) :: k
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: vr (nr)
      Real (8), Intent (Inout) :: eval
      Real (8), Intent (Out) :: g0 (nr)
      Real (8), Intent (Out) :: f0 (nr)
      
! local variables
      Integer, Parameter :: maxit = 500
      Integer :: kpa, it, nn, ir, maxr,nnp
      real(8) ::e_hi,e_lo
! energy convergence tolerance
      Real (8), Parameter :: eps = 1.d-10
      Real (8) :: t1, de, step_e,rm,rm0
      Real (8) :: large
      parameter (large=1d100)
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0

! automatic arrays
      Real (8) :: g1 (nr), f1 (nr), fr (nr), gr (nr), cf (3, nr),g0p(nr),f0p(nr),g1p(nr),f1p(nr)
  
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

      de = 1d0
      e_lo=-large
      e_hi= large 
      step_e=1d1
      it=0
!         write(*,*) 'bracketting'
!write(*,*) 'initial guess', eval
! integrate the Dirac equation
      do while ((it.lt.maxit).and.(abs(de).gt.eps*(1d0+abs(eval))))
        it=it+1
!        read(*,*) 
!        write(*,*) eval
        if (dirac_eq) then
          Call rdiracdme (m, kpa, eval, nr, r, vr, nn, g0, g1, f0, f1,sloppy)
        else
          Call rschroddme (m, l, 0, eval, nr, r, vr, nn, g0, g1, f0, f1)
        endif
        if (nn.lt.n-l-1) then
!           write(*,*) e_lo,e_hi,nn
!          if (e_lo.lt.eval) then
            e_lo=eval
!          endif
          if (e_hi.eq.large) then
            eval=eval+step_e
            step_e=step_e*2d0
          else
            eval=(e_lo+e_hi)*0.5d0
          endif
        elseif (nn.gt.n-l) then
!            write(*,*) e_lo,e_hi,nn
!          if (e_hi.gt.eval) then
            e_hi=eval
!          endif
          if (e_lo.eq.-large) then
            eval=eval-step_e
            step_e=step_e*2d0
          else
            eval=(e_lo+e_hi)*0.5d0
          endif
        else
!          write(*,*) e_lo,e_hi,nn
          if (dirac_eq) then
!            Call rdiracdme (m+1, kpa, eval, np, nr, r, vr, nnp, g0p, g1p, f0p, f1p,sloppy)
             Call rdiracint (m+1, kpa, eval, nr, r, vr, nnp, g0, f0, g0p, g1p, f0p, f1p,sloppy)
          else
             Call rschroddme (m+1, l, 0, eval, nr, r, vr, nnp , g0p, g1p, f0p, f1p)
             
          endif          
          ir=nr
          do while ((abs(g0(ir)).gt.1d100).or.(abs(g0p(ir)).gt.1d100))
            ir=ir-1
          enddo
          maxr=ir
          de=-g0(maxr)/g0p(maxr)
 !         write(*,*) (nn.eq.n-l-1),de
          if (nn.eq.n-l-1) then
             e_lo=eval
             if ((de.gt.0d0).and.(abs(de).lt.1d-1)) then
!            if (de.gt.0d0) then
!             if ((de.gt.0d0).and.(.not.sloppy)) then
               eval=eval+de
             elseif (e_hi.eq.large) then
               eval=eval+step_e
               step_e=step_e*2d0
             else
               eval=(e_lo+e_hi)*0.5d0
             endif
                 
          else   ! case if (nn.eq.n-l) then
             e_hi=eval 
             if ((de.lt.0d0).and.(abs(de).lt.1d-1)) then
!             if ((de.lt.0d0).and.(.not.sloppy)) then
               eval=eval+de
             elseif (e_lo.eq.-large) then
               eval=eval-step_e
               step_e=step_e*2d0
             else
               eval=(e_lo+e_hi)*0.5d0
             endif
            
          endif

        endif
!        read(*,*)
      enddo
!      write(*,*) it
!     read(*,*)

      if (abs(de).ge.eps*(1d0+abs(eval))) then
        Write (*,*)
        Write (*, '("Error(rdirac): maximum iterations exceeded")')
        Write (*,*)
        Stop
      endif

!      if (.not.sloppy) then
       g0=g0+g0p*de
       g1=g1+g1p*de
!      endif
      ir = 2
      nn=0
      do while ((ir.le.nr).and.(nn.lt.n-l-1))
        if (g0(ir-1)*g0(ir) .Lt. 0.d0) then
         nn=nn+1
        endif
        ir=ir+1
      enddo
      do while ((ir.le.nr).and.(g1(ir)*g1(ir-1).gt.0d0))
        ir=ir+1
      enddo
      ir=ir+1
      do while ((ir.lt.nr).and.(g1(ir)*g1(ir-1).gt.0d0).and.(g0(ir)*g0(ir-1).gt.0d0))
        ir=ir+1
      enddo
!      write(*,*) ir
      
      if (ir.lt.nr-1) then
!           write(*,*) ir
           if (dirac_eq) then
             Call rdiracdme (m  , kpa, eval, ir, r, vr, nn , g0 , g1, f0 , f1,sloppy)
             Call rdiracint (m+1, kpa, eval, ir, r, vr, nnp, g0, f0, g0p, g1p, f0p, f1p,sloppy)
           else
             Call rschroddme (m  , l, 0, eval, ir, r, vr, nn , g0 , g1, f0 , f1)
             Call rschroddme (m+1, l, 0, eval, ir, r, vr, nnp, g0p, g1p, f0p, f1p)
           endif
           de=-g0(ir)/g0p(ir)
           g0=g0+g0p*de
           g1=g1+g1p*de
           f0=f0+f0p*de
           f1=f1+f1p*de
!      write(*,*) ir
!      do maxr=1,nr
!        write(*,*) maxr,g0(maxr)
!      enddo
!      read(*,*)
      g0(ir+1:nr)=0d0
      f0(ir+1:nr)=0d0
      g1(ir+1:nr)=0d0
      f1(ir+1:nr)=0d0
     else
        f0=f0+f0p*de
!        f1=f1+f1p*de
     endif
    

! normalise
      if (dirac_eq) then
        Do ir = 1, nr
          fr (ir) = g0 (ir) ** 2 + f0 (ir) ** 2
        End Do
      elseif (input%groundstate%ValenceRelativity.eq."kh") then
        
        Do ir = 1, nr
!          rm=1d0/(1d0-0.5d0*alpha**2*vr(ir))
          fr (ir) = g0 (ir) ** 2 + (f0 (ir)*alpha) ** 2
        End Do
      elseif (input%groundstate%ValenceRelativity.eq."iora") then        
        Do ir = 1, nr
          rm0 = 1.d0 - 0.5d0 * (alpha**2) *  vr(ir)
          rm = rm0/(1d0 - 0.5d0*eval*(alpha**2)/rm0)
          fr (ir) = g0 (ir) ** 2 + (f0 (ir)*alpha*rm/rm0) ** 2
!          fr (ir) = g0 (ir) ** 2 + (f0 (ir)*alpha) ** 2

!          rm=1d0/(1d0-0.5d0*alpha**2*vr(ir))
!          fr (ir) = g0 (ir) ** 2 + (0.5d0*g1 (ir)*rm*alpha) ** 2
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
      if (dirac_eq) then
        f0 (:) = t1 * f0 (:)
      elseif (input%groundstate%ValenceRelativity.eq."kh") then
        f0 (:) = t1 * f0 (:)*alpha
      elseif (input%groundstate%ValenceRelativity.eq."iora") then
        Do ir = 1, nr
!          rm=1d0/(1d0-0.5d0*alpha**2*vr(ir)) 
!          f0 (ir) = t1 * 0.5d0 * g1 (ir) * rm * alpha 
          rm0 = 1.d0 - 0.5d0 * (alpha**2) *  vr(ir)
          rm = rm0/(1d0 - 0.5d0*eval*(alpha**2)/rm0)
          f0 (ir) = t1 * f0 (ir) * rm * alpha/rm0
!          f0 (ir) = t1 * f0 (ir) * alpha
!          fr (ir) = g0 (ir) ** 2 + (f0 (ir)*alpha*rm/rm0) ** 2

        enddo
      endif
      Return
End Subroutine
!EOC
