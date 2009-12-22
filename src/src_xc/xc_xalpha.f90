!
!
!
! Copyright (C) 1998-2006 ABINIT group (DCA, XG, GMR).
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: xc_xalpha
! !INTERFACE:
!
!
Subroutine xc_xalpha (n, rho, exc, vxc)
! !INPUT/OUTPUT PARAMETERS:
!   n   : number of density points (in,integer)
!   rho : charge density (in,real(n))
!   exc : exchange-correlation energy density (out,real(n))
!   vxc : exchange-correlation potential (out,real(n))
! !DESCRIPTION:
!   $X_{\alpha}$ approximation to the exchange-correlation potential and energy
!   density. See J. C. Slater, {\it Phys. Rev.} {\bf 81}, 385 (1951).
!
! !REVISION HISTORY:
!   Modified an ABINIT routine, September 2006 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Real (8), Intent (In) :: rho (n)
      Real (8), Intent (Out) :: exc (n)
      Real (8), Intent (Out) :: vxc (n)
! local variables
      Integer :: i
      Real (8), Parameter :: pi = 3.1415926535897932385d0
      Real (8), Parameter :: alpha = 1.d0
      Real (8) :: r, efac, rs, rsm1, vfac
      vfac = (1.5d0/pi) ** (2.d0/3.d0)
      efac = 0.75d0 * vfac
! loop over density points
      Do i = 1, n
         r = rho (i)
         If (r .Gt. 1.d-20) Then
            rs = (3.d0/(4.d0*pi*r)) ** (1.d0/3.d0)
            rsm1 = 1.0d0 / rs
! compute energy density
            exc (i) = - alpha * efac * rsm1
! compute potential
            vxc (i) = - alpha * vfac * rsm1
         End If
      End Do
End Subroutine
!EOC
