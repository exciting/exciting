!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: gradrfmt
! !INTERFACE:
!
!
Subroutine gradrfmt (lmax, nr, r, ld1, ld2, rfmt, grfmt)
! !INPUT/OUTPUT PARAMETERS:
!   lmax  : maximum angular momentum (in,integer)
!   nr    : number of radial mesh points (in,integer)
!   r     : radial mesh (in,real(nr))
!   ld1   : leading dimension 1 (in,integer)
!   ld2   : leading dimension 2 (in,integer)
!   rfmt  : real muffin-tin function (in,real(ld1,nr))
!   grfmt : gradient of rfmt (out,real(ld1,ld2,3))
! !DESCRIPTION:
!   Calculates the gradient of a real muffin-tin function. In other words, given
!   the real spherical harmonic expansion coefficients, $f_{lm}(r)$, of a
!   function $f({\bf r})$, the routine returns ${\bf F}_{lm}$ where
!   $$ \sum_{lm}{\bf F}_{lm}(r)R_{lm}(\hat{\bf r})=\nabla f({\bf r}), $$
!   and $R_{lm}$ is a real spherical harmonic function. This is done by first
!   converting the function to a complex spherical harmonic expansion and then
!   using the routine {\tt gradzfmt}. See routine {\tt genrlm}.
!
! !REVISION HISTORY:
!   Created August 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Integer, Intent (In) :: ld1
      Integer, Intent (In) :: ld2
      Real (8), Intent (In) :: rfmt (ld1, nr)
      Real (8), Intent (Out) :: grfmt (ld1, ld2, 3)
! local variables
      Integer :: lmmax, ir, i
! allocatable arrays
      Complex (8), Allocatable :: zfmt (:, :)
      Complex (8), Allocatable :: gzfmt (:, :, :)
      lmmax = (lmax+1) ** 2
      Allocate (zfmt(lmmax, nr))
      Allocate (gzfmt(lmmax, nr, 3))
! convert real to complex spherical harmonic expansion
      Do ir = 1, nr
         Call rtozflm (lmax, rfmt(:, ir), zfmt(:, ir))
      End Do
! compute the gradient
      Call gradzfmt (lmax, nr, r, lmmax, nr, zfmt, gzfmt)
! convert complex to real spherical harmonic expansion
      Do i = 1, 3
         Do ir = 1, nr
            Call ztorflm (lmax, gzfmt(:, ir, i), grfmt(:, ir, i))
         End Do
      End Do
      Deallocate (zfmt, gzfmt)
      Return
End Subroutine
!EOC
