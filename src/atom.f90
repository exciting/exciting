!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: atom
! !INTERFACE:
!
!
Subroutine atom (ptnucl, zn, nst, n, l, k, occ, xctype, xcgrad, np, nr, &
& r, eval, rho, vr, rwf)
! !USES:
      Use modxcifc
! !INPUT/OUTPUT PARAMETERS:
!   ptnucl : .true. if the nucleus is a point particle (in,logical)
!   zn     : nuclear charge (in,real)
!   nst    : number of states to solve for (in,integer)
!   n      : priciple quantum number of each state (in,integer(nst))
!   l      : quantum number l of each state (in,integer(nst))
!   k      : quantum number k (l or l+1) of each state (in,integer(nst))
!   occ    : occupancy of each state (inout,real(nst))
!   xctype : exchange-correlation type (in,integer)
!   xcgrad : 1 for GGA functional, 0 otherwise (in,integer)
!   np     : order of predictor-corrector polynomial (in,integer)
!   nr     : number of radial mesh points (in,integer)
!   r      : radial mesh (in,real(nr))
!   eval   : eigenvalue without rest-mass energy for each state (out,real(nst))
!   rho    : charge density (out,real(nr))
!   vr     : self-constistent potential (out,real(nr))
!   rwf    : major and minor components of radial wavefunctions for each state
!            (out,real(nr,2,nst))
! !DESCRIPTION:
!   Solves the Dirac-Kohn-Sham equations for an atom using the
!   exchange-correlation functional {\tt xctype} and returns the self-consistent
!   radial wavefunctions, eigenvalues, charge densities and potentials. The
!   variable {\tt np} defines the order of polynomial used for performing
!   numerical integration. Requires the exchange-correlation interface routine
!   {\tt xcifc}.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!   Fixed s.c. convergence problem, October 2003 (JKD)
!   Added support for GGA functionals, June 2006 (JKD)
!
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: ptnucl
      Real (8), Intent (In) :: zn
      Integer, Intent (In) :: nst
      Integer, Intent (In) :: n (nst)
      Integer, Intent (In) :: l (nst)
      Integer, Intent (In) :: k (nst)
      Real (8), Intent (Inout) :: occ (nst)
      Integer, Intent (In) :: xctype(3)
      Integer, Intent (In) :: xcgrad
      Integer, Intent (In) :: np
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (Out) :: eval (nst)
      Real (8), Intent (Out) :: rho (nr)
      Real (8), Intent (Out) :: vr (nr)
      Real (8), Intent (Out) :: rwf (nr, 2, nst)
      Integer, Parameter :: maxscl = 200
      Integer :: ir, ist, iscl
      Real (8), Parameter :: fourpi = 12.566370614359172954d0
! fine-structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
! potential convergence tolerance
      Real (8), Parameter :: eps = 1.d-6
      Real (8) :: sum, dv, dvp, ze, beta, t1
! allocatable arrays
      Real (8), Allocatable :: vn (:), vh (:), ex (:), ec (:), vx (:), &
     & vc (:), vrp (:)
      Real (8), Allocatable :: ri (:), fr1 (:), fr2 (:), gr1 (:), gr2 &
     & (:), cf (:, :)
      Real (8), Allocatable :: grho (:), g2rho (:), g3rho (:)
      If (nst .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(atom): invalid nst : ", I8)') nst
         Write (*,*)
         Stop
      End If
      If (np .Lt. 2) Then
         Write (*,*)
         Write (*, '("Error(atom): np < 2 : ", I8)') np
         Write (*,*)
         Stop
      End If
      If (nr .Lt. np) Then
         Write (*,*)
         Write (*, '("Error(atom): nr < np : ", 2I8)') nr, np
         Write (*,*)
         Stop
      End If
! allocate local arrays
      Allocate (vn(nr), vh(nr), ex(nr), ec(nr), vx(nr), vc(nr), &
     & vrp(nr))
      Allocate (ri(nr), fr1(nr), fr2(nr), gr1(nr), gr2(nr), cf(3, nr))
      If (xcgrad .Eq. 1) Then
         Allocate (grho(nr), g2rho(nr), g3rho(nr))
      End If
! find total electronic charge
      ze = 0.d0
      Do ist = 1, nst
         ze = ze + occ (ist)
      End Do
! set up nuclear potential
      Call potnucl (ptnucl, nr, r, zn, vn)
      Do ir = 1, nr
         ri (ir) = 1.d0 / r (ir)
! initialise the effective potential to the nuclear potential
         vr (ir) = vn (ir)
      End Do
      dvp = 0.d0
      vrp (:) = 0.d0
! initialise mixing parameter
      beta = 0.5d0
! initialise eigenvalues to relativistic values (minus the rest mass energy)
      Do ist = 1, nst
         t1 = Sqrt (dble(k(ist)**2)-(zn*alpha)**2)
         t1 = (dble(n(ist)-Abs(k(ist)))+t1) ** 2
         t1 = 1.d0 + ((zn*alpha)**2) / t1
         eval (ist) = (1.d0/alpha**2) / Sqrt (t1) - 1.d0 / alpha ** 2
      End Do
! start of self-consistent loop
      Do iscl = 1, maxscl
! solve the Dirac equation for each state
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
         Do ist = 1, nst
            Call rdirac (n(ist), l(ist), k(ist), np, nr, r, vr, &
           & eval(ist), rwf(:, 1, ist), rwf(:, 2, ist))
         End Do
!$OMP END DO
!$OMP END PARALLEL
! compute the charge density
         Do ir = 1, nr
            sum = 0.d0
            Do ist = 1, nst
               sum = sum + occ (ist) * (rwf(ir, 1, ist)**2+rwf(ir, 2, &
              & ist)**2)
            End Do
            fr1 (ir) = sum
            fr2 (ir) = sum * ri (ir)
            rho (ir) = (1.d0/fourpi) * sum * ri (ir) ** 2
         End Do
         Call fderiv (-1, nr, r, fr1, gr1, cf)
         Call fderiv (-1, nr, r, fr2, gr2, cf)
! find the Hartree potential
         t1 = gr2 (nr)
         Do ir = 1, nr
            vh (ir) = gr1 (ir) * ri (ir) + t1 - gr2 (ir)
         End Do
! normalise charge density and potential
         t1 = ze / gr1 (nr)
         rho (:) = t1 * rho (:)
         vh (:) = t1 * vh (:)
! compute the exchange-correlation energy and potential
         If (xcgrad .Eq. 1) Then
! GGA functional
! |grad rho|
            Call fderiv (1, nr, r, rho, grho, cf)
! grad^2 rho
            Call fderiv (2, nr, r, rho, g2rho, cf)
            Do ir = 1, nr
               g2rho (ir) = g2rho (ir) + 2.d0 * ri (ir) * grho (ir)
            End Do
! approximate (grad rho).(grad |grad rho|)
            Do ir = 1, nr
               g3rho (ir) = grho (ir) * g2rho (ir)
            End Do
            Call xcifc (xctype, n=nr, rho=rho, grho=grho, g2rho=g2rho, &
           & g3rho=g3rho, ex=ex, ec=ec, vx=vx, vc=vc)
         Else
! LDA functional
            Call xcifc (xctype, n=nr, rho=rho, ex=ex, ec=ec, vx=vx, &
           & vc=vc)
         End If
! self-consistent potential
         vr (:) = vh (:) + vx (:) + vc (:)
! determine change in potential
         sum = 0.d0
         Do ir = 1, nr
            sum = sum + (vr(ir)-vrp(ir)) ** 2
         End Do
         dv = Sqrt (sum) / dble (nr)
         If (iscl .Gt. 2) Then
! reduce beta if change in potential is diverging
            If (dv .Gt. dvp) beta = beta * 0.8d0
            beta = Max (beta, 0.01d0)
         End If
         dvp = dv
         Do ir = 1, nr
! mix old and new potentials
            vr (ir) = (1.d0-beta) * vrp (ir) + beta * vr (ir)
            vrp (ir) = vr (ir)
! add nuclear potential
            vr (ir) = vr (ir) + vn (ir)
         End Do
! check for convergence
         If ((iscl .Gt. 2) .And. (dv .Lt. eps)) Go To 10
! end self-consistent loop
      End Do
      Write (*,*)
      Write (*, '("Warning(atom): maximum iterations exceeded")')
10    Continue
      Deallocate (vn, vh, ex, ec, vx, vc, vrp)
      Deallocate (ri, fr1, fr2, gr1, gr2, cf)
      If (xcgrad .Eq. 1) Then
         Deallocate (grho, g2rho, g3rho)
      End If
      Return
End Subroutine
!EOC
