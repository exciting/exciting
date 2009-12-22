!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_fxc_lrcd
      Implicit None
!
Contains
!
!BOP
! !ROUTINE: fxc_lrcd
! !INTERFACE:
!
!
      Subroutine fxc_lrcd (msiz, sw, alpha, beta, w, fxc)
! !USES:
         Use mod_constants, Only: fourpi
         Use modxs, Only: unitout
! !INPUT/OUTPUT PARAMETERS:
!   msiz  : matrix size of local field effects (in,integer)
!   sw    : true for inclusion of local field effects (in,logical)
!   alpha : real constant (in,real)
!   w     : frequency grid (in,complex(:))
!   fxc   : xc-kernel Fourier coefficients (out,complex(:,:))
! !DESCRIPTION:
!   Dynamical long range xc-kernel; S. Botti, PRB 72, 125203 (2005).
!   Calculates the symmetrized xc-kernel for the static long range model.
!   According to the switch {\tt sw} either $$
!         f_{\rm xc}({\bf G},{\bf G'}) = - \frac{1}{4\pi} (\alpha+\beta\omega^2)
!         \delta({\bf G},{\bf G'}), $$
!   if the switch is true, or $$
!         f_{\rm xc}({\bf G},{\bf G'}) = - \frac{1}{4\pi} (\alpha+\beta\omega^2)
!         \delta({\bf G},{\bf G'})\delta({\bf G},{\bf 0}), $$
!   otherwise.
!
! !REVISION HISTORY:
!   Created March 2006 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Intent (In) :: msiz
    ! true if all G-components of fxc are to be considered
         Logical, Intent (In) :: sw
         Real (8), Intent (In) :: alpha, beta
         Complex (8), Intent (In) :: w
         Complex (8), Intent (Out) :: fxc (:, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'fxc_lrcd'
         Complex (8) :: zt1
         Integer :: sh (2), ig
         sh = shape (fxc)
         If ((sh(1) .Lt. msiz) .Or. (sh(2) .Lt. msiz)) Then
            Write (unitout, '(a, 2i9, a, i9, a)') 'Error(' // trim &
           & (thisnam) // '): size of fxc is to small (required)', sh, &
           & '(', msiz, ')'
            Call terminate
         End If
         fxc (:, :) = (0.d0, 0.d0)
         zt1 = - (alpha+beta*w**2) / fourpi
         If ( .Not. sw) Then
            fxc (1, 1) = zt1
         Else
            Do ig = 1, msiz
               fxc (ig, ig) = zt1
            End Do
         End If
      End Subroutine fxc_lrcd
!EOC
!
End Module m_fxc_lrcd
