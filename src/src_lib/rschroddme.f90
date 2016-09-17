!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rschroddme
! !INTERFACE:
!
!
Subroutine rschroddme (m, l, k, e, nr, r, vr, nn, p0, p1, q0, q1)
use mod_timing
use modinput
! !INPUT/OUTPUT PARAMETERS:
!   m   : order of energy derivative (in,integer)
!   l   : angular momentum quantum number (in,integer)
!   k   : quantum number k, zero if Dirac eqn. is not to be used (in,integer)
!   e   : energy (in,real)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   nn  : number of nodes (out,integer)
!   p0  : m th energy derivative of P (out,real(nr))
!   p1  : radial derivative of p0 (out,real(nr))
!   q0  : m th energy derivative of Q (out,real(nr))
!   q1  : radial derivative of q0 (out,real(nr))
! !DESCRIPTION:
!   Finds the solution to the $m$th energy derivative of the scalar relativistic
!   radial Schr\"{o}dinger equation using the routine {\tt rschrodint}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
      Integer, Intent (In) :: m
      Integer, Intent (In) :: l
      Integer, Intent (In) :: k
      Real (8), Intent (In) :: e
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: vr (nr)
      Integer, Intent (Out) :: nn
      Real (8), Intent (Out) :: p0 (nr)
      Real (8), Intent (Out) :: p1 (nr)
      Real (8), Intent (Out) :: q0 (nr)
      Real (8), Intent (Out) :: q1 (nr)
! local variables
      Integer :: im, kpa, ir
! fine-structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
      Real (8) :: rm,rm0,rmfactor
      character(16) :: relativity
! allocatable arrays
      Real (8), Allocatable :: p0p (:),q0p(:),pe(:,:),qe(:,:)
      Real (8), Allocatable :: g0 (:), g1 (:)
      Real (8), Allocatable :: f0 (:), f1 (:)
      Real (8), Allocatable :: cf (:, :)
! timer
      Real (8) ts0, ts1

      call timesec(ts0)
      If (nr .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(rschroddme): invalid nr : ", I8)') nr
         Write (*,*)
         Stop
      End If
      If ((m .Lt. 0) .Or. (m .Gt. 6)) Then
         Write (*,*)
         Write (*, '("Error(rschroddme): m out of range : ", I8)') m
         Write (*,*)
         Stop
      End If

#ifdef SPECIES
      relativity = "zora"
#else
      relativity = input%groundstate%ValenceRelativity
#endif

      If (k .Eq. 0) Then
         if ((relativity.eq."zora").or.(relativity.eq."none")) then
           if (relativity.eq."zora") then
             rmfactor=1d0
           else
             rmfactor=0d0
           endif
           Allocate (p0p(nr))
           If (m .Eq. 0) Then
             Call rschrodint (m, l, e, nr, r, vr, nn, rmfactor, p0p, p0, p1, q0, q1)
           Else
             Do im = 0, m
               Call rschrodint (im, l, e, nr, r, vr, nn, rmfactor, p0p, p0, p1, q0, q1)
               p0p (:) = p0 (:)
             End Do
           End If
           Deallocate (p0p)
         elseif ((relativity.eq."kh*").or.(relativity.eq."kh")) then
! Koelling-Harmon with or without small component
           Allocate (p0p(nr))
           Allocate (q0p(nr))
           if (m.eq.0) then 
             call rkhint (0, l, e, nr, r, vr, nn, p0p, q0p, p0, p1, q0, q1)
           elseif (m.eq.1) then
             call rkhint (0, l, e, nr, r, vr, nn, p0p, q0p, p0, p1, q0, q1)
             p0p (:) = p0 (:)
             q0p (:) = q0 (:)
             call rkhint (1, l, e, nr, r, vr, nn, p0p, q0p, p0, p1, q0, q1)
           else
             write(*,*) 'Error(rschroddme): energy derivative',m,'not implemented.'
             stop
           endif
           Deallocate (p0p,q0p)
         elseif ((relativity.eq."iora").or.(relativity.eq."iora*")) then
! IORA with or without small component
           Allocate (p0p(nr),q0p(nr))
          
           if (m.le.3) then
             call rlkhint (0, l, e, nr, r, vr, nn, p0p, q0p, p0, p1, q0, q1)
             if (m.gt.0) then
               Allocate(pe(nr,m),qe(nr,m))
               pe(:,1)=p0(:)
               qe(:,1)=q0(:)
               do ir=1,nr
                 rm0 = 1.d0 - 0.5d0 * (alpha**2) *  vr(ir)
                 rm = rm0/(1d0 - 0.5d0*e*(alpha**2)/rm0)
                 q0p(ir)=(alpha*rm/rm0)**2*qe(ir,1)
                 p0p(ir)=(1d0+dble(l*(l+1))*(alpha/(2d0*rm0*r(ir)))**2)*pe(ir,1)
               enddo
               call rlkhint (1, l, e, nr, r, vr, nn, p0p, q0p, p0, p1, q0, q1)
               if (m.gt.1) then
                 pe(:,2)=p0(:)
                 qe(:,2)=q0(:)
                 do ir=1,nr
                   rm0 = 1.d0 - 0.5d0 * (alpha**2) *  vr(ir)
                   rm = rm0/(1d0 - 0.5d0*e*(alpha**2)/rm0)
                   q0p(ir)=2d0*(alpha*rm/rm0)**2*qe(ir,2)+(alpha*rm/rm0)**4/rm*qe(ir,1)
                   p0p(ir)=2d0*(1d0+dble(l*(l+1))*(alpha/(2d0*rm0*r(ir)))**2)*pe(ir,2)
                 enddo
                 call rlkhint (2, l, e, nr, r, vr, nn, p0p, q0p, p0, p1, q0, q1)
                 if (m.gt.2) then
                   pe(:,3)=p0(:)
                   qe(:,3)=q0(:)
                   do ir=1,nr
                     rm0 = 1.d0 - 0.5d0 * (alpha**2) *  vr(ir)
                     rm = rm0/(1d0 - 0.5d0*e*(alpha**2)/rm0)
                     q0p(ir)=3d0*(alpha*rm/rm0)**2*qe(ir,3)+2d0*(alpha*rm/rm0)**4/rm*qe(ir,2)+1.5d0*(alpha*rm/rm0)**6/(rm**2)*qe(ir,1)
                     p0p(ir)=3d0*(1d0+dble(l*(l+1))*(alpha/(2d0*rm0*r(ir)))**2)*pe(ir,3)
                   enddo
                   call rlkhint (3, l, e, nr, r, vr, nn, p0p, q0p, p0, p1, q0, q1)
                 endif 
               endif
               deallocate(pe,qe)
             endif
           else
             write(*,*) 'Error(rschroddme): energy derivative',m,'not implemented.'
             stop
           endif
           Deallocate (p0p,q0p)
         else  
           write(*,*) 'Error(rschroddme):',trim(relativity),' not implemented.'
           write(*,*) 'case sensitivity issue?'
           stop
         endif
      Else
! use the Dirac equation
         Allocate (g0(nr), g1(nr))
         Allocate (f0(nr), f1(nr))
         Allocate (cf(3, nr))
         If (k .Eq. l) Then
            kpa = k
         Else If (k .Eq. l+1) Then
            kpa = - k
         Else
            Write (*,*)
            Write (*, '("Error(rschroddme): incompatible l and k : ", 2&
           &I8)') l, k
            Write (*,*)
            Stop
         End If
         Call rdiracdme (m, kpa, e, nr, r, vr, nn, g0, g1, f0, f1)
! determine equivalent scalar relativistic functions
         Do ir = 1, nr
            rm = 1.d0 - 0.5d0 * (alpha**2) * vr (ir)
            p0 (ir) = g0 (ir)
            p1 (ir) = g1 (ir)
            q0 (ir) = (p1(ir)-p0(ir)/r(ir)) / (2.d0*rm)
         End Do
         Call fderiv (1, nr, r, q0, q1, cf)
         Deallocate (g0, g1, f0, f1, cf)
      End If
      call timesec(ts1)
      time_rschrod=ts1-ts0+time_rschrod
      Return
End Subroutine
!EOC
