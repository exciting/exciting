!
!
!
! Copyright (C) 2003-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rfmtinp
! !INTERFACE:
Real (8) Function rfmtinp (lrstp, lmax, nr, r, ld, rfmt1, rfmt2)
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   lmax  : maximum angular momentum (in,integer)
!   nr    : number of radial mesh points (in,integer)
!   r     : radial mesh (in,real(nr))
!   ld    : the leading dimension (in,integer)
!   rfmt1 : first real function inside muffin-tin (in,real(ld,nr))
!   rfmt2 : second real function inside muffin-tin (in,real(ld,nr))
! !DESCRIPTION:
!   Calculates the inner product of two real fuctions in the muffin-tin. So
!   given two real functions of the form
!   $$ f({\bf r})=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}f_{lm}(r)R_{lm}
!    (\hat{\bf r}) $$
!   where $R_{lm}$ are the real spherical harmonics, the function returns
!   $$ I=\int\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}f_{lm}^1(r)f_{lm}^2(r)r^2
!    dr\;. $$
!   The radial integral is performed using a high accuracy cubic spline method.
!   See routines {\tt genrlm} and {\tt fderiv}.
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lrstp
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Integer, Intent (In) :: ld
      Real (8), Intent (In) :: rfmt1 (ld, nr)
      Real (8), Intent (In) :: rfmt2 (ld, nr)
! local variables
      Integer :: lmmax, ir, irc
! automatic arrays
      Real (8) :: rc (nr), fr (nr), gr (nr), cf (3, nr)
! external functions
      Real (8) :: ddot
      External ddot
      If (lrstp .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(rfmtinp): lrstp <= 0 : ", I8)') lrstp
         Write (*,*)
         Stop
      End If
      If (lmax .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(rfmtinp): lmax < 0 : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      If (nr .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(rfmtinp): nr <= 0 : ", I8)') nr
         Write (*,*)
         Stop
      End If
      lmmax = (lmax+1) ** 2
      irc = 0
      Do ir = 1, nr, lrstp
         irc = irc + 1
         rc (irc) = r (ir)
         fr (irc) = ddot (lmmax, rfmt1(:, ir), 1, rfmt2(:, ir), 1) * &
        & (r(ir)**2)
      End Do
      Call fderiv (-1, irc, rc, fr, gr, cf)
      rfmtinp = gr (irc)
      Return
End Function
!EOC
