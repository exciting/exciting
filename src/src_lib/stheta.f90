!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: stheta
! !INTERFACE:
Real (8) Function stheta (stype, x)
! !INPUT/OUTPUT PARAMETERS:
!   stype : smearing type (in,integer)
!   x     : real argument (in,real)
! !DESCRIPTION:
!   Returns the Heaviside step function corresponding to the smooth approximation
!   to the Dirac delta function:
!   $$ \tilde\Theta(x)=\int_{-\infty}^x dt\,\tilde\delta(t). $$
!   See function {\tt sdelta} for details.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: stype
      Real (8), Intent (In) :: x
! external functions
      Real (8) stheta_mp, stheta_fd, stheta_sq
      External stheta_mp, stheta_fd, stheta_sq
      stheta = 0.d0
      Select Case (stype)
      Case (0)
         stheta = stheta_mp (0, x)
         Return
      Case (1)
         stheta = stheta_mp (1, x)
         Return
      Case (2)
         stheta = stheta_mp (2, x)
         Return
      Case (3)
         stheta = stheta_fd (x)
         Return
      Case (4)
         stheta = stheta_sq (x)
      Case Default
         Write (*,*)
         Write (*, '("Error(stheta): sytpe not defined : ",I8)') stype
         Write (*,*)
         Stop
      End Select
End Function
!EOC
