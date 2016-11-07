!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_fxc_lrc
      use modmpi
      Implicit None
!
Contains
!
!BOP
! !ROUTINE: fxc_lrc
! !INTERFACE:
!
!
      Subroutine fxc_lrc (msiz, sw, alpha, fxc)
! !USES:
         Use mod_constants, Only: fourpi
         Use modxs, Only: unitout
! !INPUT/OUTPUT PARAMETERS:
!   msiz  : matrix size of local field effects (in,integer)
!   sw    : true for inclusion of local field effects (in,logical)
!   alpha : real constant (in,real)
!   fxc   : xc-kernel Fourier coefficients (out,complex(:,:))
! !DESCRIPTION:
!   Static long range xc-kernel; S. Botti, PRB 70, 045301 (2004).
!   Calculates the symmetrized xc-kernel for the static long range model.
!   According to the switch {\tt sw} either $$
!         f_{\rm xc}({\bf G},{\bf G'}) = - \frac{\alpha}{4\pi}
!         \delta({\bf G},{\bf G'}), $$
!   if the switch is true, or $$
!         f_{\rm xc}({\bf G},{\bf G'}) = - \frac{\alpha}{4\pi}
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
         Real (8), Intent (In) :: alpha
         Complex (8), Intent (Out) :: fxc (:, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'fxc_lrc'
         Real (8) :: t1
         Integer :: sh (2), ig
         sh = shape (fxc)
         If ((sh(1) .Lt. msiz) .Or. (sh(2) .Lt. msiz)) Then
            Write (unitout, '(a, 2i9, a, i9, a)') 'Error(' // trim &
           & (thisnam) // '): size of fxc is to small (required)', sh, &
           & '(', msiz, ')'
            Call terminate
         End If
         fxc (:, :) = (0.d0, 0.d0)
         t1 = - alpha / fourpi
         If ( .Not. sw) Then
            fxc (1, 1) = t1
         Else
            Do ig = 1, msiz
               fxc (ig, ig) = t1
            End Do
         End If
      End Subroutine fxc_lrc
End Module m_fxc_lrc
!EOC
