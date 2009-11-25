!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module modxcifc
Contains
!
!BOP
! !ROUTINE: xcifc
! !INTERFACE:
!
!
      Subroutine xcifc (xctype, n, rho, rhoup, rhodn, grho, gup, gdn, &
     & g2rho, g2up, g2dn, g3rho, g3up, g3dn, ex, ec, vx, vc, vxup, &
     & vxdn, vcup, vcdn)
! !INPUT/OUTPUT PARAMETERS:
!   xctype : type of exchange-correlation functional (in,integer)
!   n      : number of density points (in,integer,optional)
!   rho    : spin-unpolarised charge density (in,real(n),optional)
!   rhoup  : spin-up charge density (in,real(n),optional)
!   rhodn  : spin-down charge density (in,real(n),optional)
!   grho   : |grad rho| (in,real(n),optional)
!   gup    : |grad rhoup| (in,real(n),optional)
!   gdn    : |grad rhodn| (in,real(n),optional)
!   g2rho  : grad^2 rho (in,real(n),optional)
!   g2up   : grad^2 rhoup (in,real(n),optional)
!   g2dn   : grad^2 rhodn (in,real(n),optional)
!   g3rho  : (grad rho).(grad |grad rho|) (in,real(n),optional)
!   g3up   : (grad rhoup).(grad |grad rhoup|) (in,real(n),optional)
!   g3dn   : (grad rhodn).(grad |grad rhodn|) (in,real(n),optional)
!   ex     : exchange energy density (out,real(n),optional)
!   ec     : correlation energy density (out,real(n),optional)
!   vx     : spin-unpolarised exchange potential (out,real(n),optional)
!   vc     : spin-unpolarised correlation potential (out,real(n),optional)
!   vxup   : spin-up exchange potential (out,real(n),optional)
!   vxdn   : spin-down exchange potential (out,real(n),optional)
!   vcup   : spin-up correlation potential (out,real(n),optional)
!   vcdn   : spin-down correlation potential (out,real(n),optional)
! !DESCRIPTION:
!   Interface to the exchange-correlation routines. This makes it relatively
!   simple to add new functionals which do not necessarily depend only on
!   $\rho$.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
         Implicit None
! mandatory arguments
         Integer, Intent (In) :: xctype
! optional arguments
         Integer, Optional, Intent (In) :: n
         Real (8), Optional, Intent (In) :: rho (*)
         Real (8), Optional, Intent (In) :: rhoup (*)
         Real (8), Optional, Intent (In) :: rhodn (*)
         Real (8), Optional, Intent (In) :: grho (*)
         Real (8), Optional, Intent (In) :: gup (*)
         Real (8), Optional, Intent (In) :: gdn (*)
         Real (8), Optional, Intent (In) :: g2rho (*)
         Real (8), Optional, Intent (In) :: g2up (*)
         Real (8), Optional, Intent (In) :: g2dn (*)
         Real (8), Optional, Intent (In) :: g3rho (*)
         Real (8), Optional, Intent (In) :: g3up (*)
         Real (8), Optional, Intent (In) :: g3dn (*)
         Real (8), Optional, Intent (Out) :: ex (*)
         Real (8), Optional, Intent (Out) :: ec (*)
         Real (8), Optional, Intent (Out) :: vx (*)
         Real (8), Optional, Intent (Out) :: vc (*)
         Real (8), Optional, Intent (Out) :: vxup (*)
         Real (8), Optional, Intent (Out) :: vxdn (*)
         Real (8), Optional, Intent (Out) :: vcup (*)
         Real (8), Optional, Intent (Out) :: vcdn (*)
! local variables
         Real (8) :: kappa, mu, beta
! automatic arrays
         Real (8), Allocatable :: ra (:, :)
         Select Case (Abs(xctype))
         Case (1)
! No density-derived exchange-correlation energy or potential
            If ( .Not. (present(n))) Go To 10
            If (n .Le. 0) Go To 20
            If (present(ex)) ex (1:n) = 0.d0
            If (present(ec)) ec (1:n) = 0.d0
            If (present(vx)) vx (1:n) = 0.d0
            If (present(vc)) vc (1:n) = 0.d0
            If (present(vxup)) vxup (1:n) = 0.d0
            If (present(vxdn)) vxdn (1:n) = 0.d0
            If (present(vcup)) vcup (1:n) = 0.d0
            If (present(vcdn)) vcdn (1:n) = 0.d0
         Case (2)
! Perdew-Zunger parameterisation of Ceperley-Alder electron gas
! J. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981)
! D.M. Ceperly and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980)
            If (present(n) .And. present(rho) .And. present(ex) .And. &
           & present(ec) .And. present(vx) .And. present(vc)) Then
               If (n .Le. 0) Go To 20
               Call xc_pzca (n, rho, ex, ec, vx, vc)
            Else
               Go To 10
            End If
         Case (3)
! Perdew-Wang parameterisation of the spin-polarised Ceperley-Alder electron gas
! J. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)
! D.M. Ceperly and B.J. Alder, Phys. Rev. Lett. 45, 566 (1980)
            If (present(n) .And. present(rhoup) .And. present(rhodn) &
           & .And. present(ex) .And. present(ec) .And. present(vxup) &
           & .And. present(vxdn) .And. present(vcup) .And. &
           & present(vcdn)) Then
! spin-polarised density
               If (n .Le. 0) Go To 20
               Call xc_pwca (n, rhoup, rhodn, ex, ec, vxup, vxdn, vcup, &
              & vcdn)
            Else If (present(n) .And. present(rho) .And. present(ex) &
           & .And. present(ec) .And. present(vx) .And. present(vc)) &
           & Then
! divide spin-unpolarised density into up and down
               If (n .Le. 0) Go To 20
               Allocate (ra(n, 1))
               ra (1:n, 1) = 0.5d0 * rho (1:n)
               Call xc_pwca (n, ra(:, 1), ra(:, 1), ex, ec, vx, vx, vc, &
              & vc)
               Deallocate (ra)
            Else
               Go To 10
            End If
         Case (4)
! X-alpha approximation
! J. C. Slater, Phys. Rev. 81, 385 (1951)
            If (present(n) .And. present(rho) .And. present(ex) .And. &
           & present(ec) .And. present(vx) .And. present(vc)) Then
               If (n .Le. 0) Go To 20
               Call xc_xalpha (n, rho, ex, vx)
! set correlation energy and potential to zero
               ec (1:n) = 0.d0
               vc (1:n) = 0.d0
            Else
               Go To 10
            End If
         Case (5)
! U. von Barth and L. Hedin parameterisation of LSDA
! J. Phys. C, 5, 1629 (1972)
            If (present(n) .And. present(rhoup) .And. present(rhodn) &
           & .And. present(ex) .And. present(ec) .And. present(vxup) &
           & .And. present(vxdn) .And. present(vcup) .And. &
           & present(vcdn)) Then
! spin-polarised density
               If (n .Le. 0) Go To 20
               Call xc_vbh (n, rhoup, rhodn, ex, ec, vxup, vxdn, vcup, &
              & vcdn)
            Else If (present(n) .And. present(rho) .And. present(ex) &
           & .And. present(ec) .And. present(vx) .And. present(vc)) &
           & Then
! divide spin-unpolarised density into up and down
               If (n .Le. 0) Go To 20
               Allocate (ra(n, 1))
               ra (1:n, 1) = 0.5d0 * rho (1:n)
               Call xc_vbh (n, ra(:, 1), ra(:, 1), ex, ec, vx, vx, vc, &
              & vc)
               Deallocate (ra)
            Else
               Go To 10
            End If
         Case (20, 21, 22)
! original PBE kappa
            kappa = 0.804d0
            If (xctype .Eq. 21) Then
! Zhang-Yang kappa
               kappa = 1.245d0
            End If
! original PBE mu and beta
            mu = 0.2195149727645171d0
            beta = 0.06672455060314922d0
            If (xctype .Eq. 22) Then
! PBEsol parameters
               mu = 10.d0 / 81.d0
               beta = 0.046d0
            End If
! Perdew-Burke-Ernzerhof generalised gradient approximation
! Phys. Rev. Lett. 77, 3865 (1996); 78, 1396(E) (1997)
! Revised PBE, Zhang-Yang, Phys. Rev. Lett. 80, 890 (1998)
            If (present(n) .And. present(rhoup) .And. present(rhodn) &
           & .And. present(grho) .And. present(gup) .And. present(gdn) &
           & .And. present(g2up) .And. present(g2dn) .And. &
           & present(g3rho) .And. present(g3up) .And. present(g3dn) &
           & .And. present(ex) .And. present(ec) .And. present(vxup) &
           & .And. present(vxdn) .And. present(vcup) .And. &
           & present(vcdn)) Then
               Call xc_pbe (n, kappa, mu, beta, rhoup, rhodn, grho, &
              & gup, gdn, g2up, g2dn, g3rho, g3up, g3dn, ex, ec, vxup, &
              & vxdn, vcup, vcdn)
            Else If (present(n) .And. present(rho) .And. present(grho) &
           & .And. present(g2rho) .And. present(g3rho) .And. &
           & present(ex) .And. present(ec) .And. present(vx) .And. &
           & present(vc)) Then
               Allocate (ra(n, 4))
               ra (1:n, 1) = 0.5d0 * rho (1:n)
               ra (1:n, 2) = 0.5d0 * grho (1:n)
               ra (1:n, 3) = 0.5d0 * g2rho (1:n)
               ra (1:n, 4) = 0.25d0 * g3rho (1:n)
               Call xc_pbe (n, kappa, mu, beta, ra(:, 1), ra(:, 1), &
              & grho, ra(:, 2), ra(:, 2), ra(:, 3), ra(:, 3), g3rho, &
              & ra(:, 4), ra(:, 4), ex, ec, vx, vx, vc, vc)
               Deallocate (ra)
            Else
               Go To 10
            End If
         Case (26)
! Wu-Cohen exchange with PBE correlation generalised gradient functional
! Zhigang Wu and R. E. Cohen, Phys. Rev. B 73, 235116 (2006)
            If (present(n) .And. present(rho) .And. present(grho) .And. &
           & present(g2rho) .And. present(g3rho) .And. present(ex) &
           & .And. present(ec) .And. present(vx) .And. present(vc)) &
           & Then
               Call xc_wc06 (n, rho, grho, g2rho, g3rho, ex, ec, vx, &
              & vc)
            Else
               Go To 10
            End If
         Case (30)
! Armiento-Mattsson generalised gradient functional
! R. Armiento and A. E. Mattsson, Phys. Rev. B 72, 085108 (2005)
            If (present(n) .And. present(rho) .And. present(grho) .And. &
           & present(g2rho) .And. present(g3rho) .And. present(ex) &
           & .And. present(ec) .And. present(vx) .And. present(vc)) &
           & Then
               Call xc_am05 (n, rho, grho, g2rho, g3rho, ex, ec, vx, &
              & vc)
            Else
               Go To 10
            End If
         Case Default
            Write (*,*)
            Write (*, '("Error(xcifc): xctype not defined : ", I8)') &
           & xctype
            Write (*,*)
            Stop
         End Select
! set exchange potential to zero for EXX
         If (xctype .Le.-2) Then
            If (present(vx)) vx (1:n) = 0.d0
            If (present(vxup)) vxup (1:n) = 0.d0
            If (present(vxdn)) vxdn (1:n) = 0.d0
         End If
         Return
10       Continue
         Write (*,*)
         Write (*, '("Error(xcifc): missing arguments for exchange-corr&
        &elation type ", I5)') xctype
         Write (*,*)
         Stop
20       Continue
         Write (*,*)
         Write (*, '("Error(xcifc): n <= 0 : ", I8)') n
         Write (*,*)
         Stop
      End Subroutine
!EOC
!
!BOP
! !ROUTINE: getxcdata
! !INTERFACE:
!
!
      Subroutine getxcdata (xctype, xcdescr, xcspin, xcgrad)
! !INPUT/OUTPUT PARAMETERS:
!   xctype  : type of exchange-correlation functional (in,integer)
!   xcdescr : description of functional (out,character(256))
!   xcspin  : spin treatment (out,integer)
!   xcgrad  : gradient treatment (out,integer)
! !DESCRIPTION:
!   Returns data on the exchange-correlation functional labeled by {\tt xctype}.
!   The character array {\tt xctype} contains a short description of the
!   functional including journal references. The variable {\tt xcspin} is set to
!   1 or 0 for spin-polarised or -unpolarised functionals, respectively. For
!   functionals which require the gradients of the density {\tt xcgrad} is set
!   to 1, otherwise it is set to 0.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
         Implicit None
         Integer, Intent (In) :: xctype
         Character (256), Intent (Out) :: xcdescr
         Integer, Intent (Out) :: xcspin
         Integer, Intent (Out) :: xcgrad
         Select Case (Abs(xctype))
         Case (1)
            xcdescr = 'No density-derived exchange-correlation energy o&
           &r potential'
! spin-polarisation or gradient status not required
            xcspin = - 1
            xcgrad = - 1
            Return
         Case (2)
            xcdescr = 'Perdew-Zunger/Ceperley-Alder, Phys. Rev. B 23, 5&
           &048 (1981)'
            xcspin = 0
            xcgrad = 0
            Return
         Case (3)
            xcdescr = 'Perdew-Wang/Ceperley-Alder, Phys. Rev. B 45, 132&
           &44 (1992)'
            xcspin = 1
            xcgrad = 0
         Case (4)
            xcdescr = 'X-alpha approximation, J. C. Slater, Phys. Rev. &
           &81, 385 (1951)'
            xcspin = 0
            xcgrad = 0
         Case (5)
            xcdescr = 'von Barth-Hedin, J. Phys. C 5, 1629 (1972)'
            xcspin = 1
            xcgrad = 0
         Case (20)
            xcdescr = 'Perdew-Burke-Ernzerhof, Phys. Rev. Lett. 77, 386&
           &5 (1996)'
            xcspin = 1
            xcgrad = 1
         Case (21)
            xcdescr = 'Revised PBE, Zhang-Yang, Phys. Rev. Lett. 80, 89&
           &0 (1998)'
            xcspin = 1
            xcgrad = 1
         Case (22)
            xcdescr = 'PBEsol, arXiv:0711.0156 (2007)'
            xcspin = 1
            xcgrad = 1
         Case (26)
            xcdescr = 'Wu-Cohen exchange + PBE correlation, Phys. Rev. &
           &B 73, 235116 (2006)'
            xcspin = 0
            xcgrad = 1
         Case (30)
            xcdescr = 'Armiento-Mattsson functional, Phys. Rev. B 72, 8&
           &5108 (2005)'
            xcspin = 0
            xcgrad = 1
         Case Default
            Write (*,*)
            Write (*, '("Error(getxcdata): xctype not defined : ", I8)') xctype
            Write (*,*)
            Stop
         End Select
         Return
      End Subroutine
!EOC
!
End Module
