!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: zfmtinp
! !INTERFACE:
Complex (8) Function zfmtinp (tsh, lmax, nr, r, ld, zfmt1, zfmt2)
! !INPUT/OUTPUT PARAMETERS:
!   tsh   : .true. if the functions are in spherical harmonics (in,logical)
!   lmax  : maximum angular momentum
!   nr    : number of radial mesh points (in,integer)
!   r     : radial mesh (in,real(nr))
!   ld    : leading dimension (in,integer)
!   zfmt1 : first complex muffin-tin function in spherical harmonics/
!           coordinates (in,complex(ld,nr))
!   zfmt2 : second complex muffin-tin function in spherical harmonics/
!           coordinates (in,complex(ld,nr))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions in the muffin-tin. In
!   other words, given two complex functions of the form
!   $$ f({\bf r})=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}f_{lm}(r)Y_{lm}
!    (\hat{\bf r}), $$
!   the function returns
!   $$ I=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}\int f_{lm}^{1*}(r)
!    f_{lm}^2(r)r^2\,dr\;. $$
!   Note that if {\tt tsh} is {\tt .false.} the functions are in spherical
!   coordinates rather than spherical harmonics. In this case $I$ is multiplied
!   by $4\pi/(l_{\rm max}+1)^2$.
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: tsh
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Integer, Intent (In) :: ld
      Complex (8), Intent (In) :: zfmt1 (ld, nr)
      Complex (8), Intent (In) :: zfmt2 (ld, nr)
! local variables
      Integer :: lmmax, ir
      Real (8), Parameter :: fourpi = 12.566370614359172954d0
      Real (8) :: t1, t2
      Complex (8) zt1
! automatic arrays
      Real (8) :: fr1 (nr), fr2 (nr), gr (nr), cf (3, nr)
! external functions
      Complex (8) zdotc
      External zdotc
      lmmax = (lmax+1) ** 2
      Do ir = 1, nr
         t1 = r (ir) ** 2
         zt1 = zdotc (lmmax, zfmt1(:, ir), 1, zfmt2(:, ir), 1)
         fr1 (ir) = t1 * dble (zt1)
         fr2 (ir) = t1 * aimag (zt1)
      End Do
      Call fderiv (-1, nr, r, fr1, gr, cf)
      t1 = gr (nr)
      Call fderiv (-1, nr, r, fr2, gr, cf)
      t2 = gr (nr)
      zfmtinp = cmplx (t1, t2, 8)
      If ( .Not. tsh) zfmtinp = zfmtinp * fourpi / dble (lmmax)
      Return
End Function
!EOC
