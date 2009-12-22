!
!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: sdelta
! !INTERFACE:
Function sdelta (stype, x)
! !INPUT/OUTPUT PARAMETERS:
!   stype : smearing type (in,integer)
!   x     : real argument (in,real)
! !DESCRIPTION:
!   Returns a normalised smooth approximation to the Dirac delta function. These
!   functions are defined such that
!   $$ \int\tilde{\delta}(x)dx=1. $$
!   The effective width, $w$, of the delta function may be varied by using the
!   normalising transformation
!   $$ \tilde{\delta}_w(x)\equiv\frac{\tilde{\delta}(x/w)}{w}. $$
!   Currently implimented are:
!   \begin{list}{}{\itemsep -2pt}
!    \item[0.] Gaussian
!    \item[1.] Methfessel-Paxton order 1
!    \item[2.] Methfessel-Paxton order 2
!    \item[3.] Fermi-Dirac
!    \item[4.] Square-wave impulse
!   \end{list}
!   See routines {\tt stheta}, {\tt sdelta\_mp}, {\tt sdelta\_fd} and
!   {\tt sdelta\_sq}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
      Real (8) :: sdelta
! arguments
      Integer, Intent (In) :: stype
      Real (8), Intent (In) :: x
! external functions
      Real (8) :: sdelta_mp, sdelta_fd, sdelta_sq
      External sdelta_mp, sdelta_fd, sdelta_sq
      sdelta = 0.d0
      Select Case (stype)
      Case (0)
         sdelta = sdelta_mp (0, x)
         Return
      Case (1)
         sdelta = sdelta_mp (1, x)
         Return
      Case (2)
         sdelta = sdelta_mp (2, x)
         Return
      Case (3)
         sdelta = sdelta_fd (x)
         Return
      Case (4)
         sdelta = sdelta_sq (x)
      Case Default
         Write (*,*)
         Write (*, '("Error(sdelta): sytpe not defined : ", I8)') stype
         Write (*,*)
         Stop
      End Select
End Function
!EOC
!
!BOP
! !ROUTINE: getsdata
! !INTERFACE:
!
!
Subroutine getsdata (stype, sdescr)
! !INPUT/OUTPUT PARAMETERS:
!   stype  : smearing type (in,integer)
!   sdescr : smearing scheme description (out,character(256))
! !DESCRIPTION:
!   Returns a description of the smearing scheme as string {\tt sdescr} up
!   to 256 characters long.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: stype
      Character (256), Intent (Out) :: sdescr
      Select Case (stype)
      Case (0)
         sdescr = 'Gaussian'
         Return
      Case (1)
         sdescr = 'Methfessel-Paxton order 1, Phys. Rev. B 40, 3616 (19&
        &89)'
         Return
      Case (2)
         sdescr = 'Methfessel-Paxton order 2, Phys. Rev. B 40, 3616 (19&
        &89)'
         Return
      Case (3)
         sdescr = 'Fermi-Dirac'
         Return
      Case (4)
         sdescr = 'Square-wave impulse'
      Case Default
         Write (*,*)
         Write (*, '("Error(getsdata): sytpe not defined : ", I8)') &
        & stype
         Write (*,*)
         Stop
      End Select
End Subroutine
!EOC
