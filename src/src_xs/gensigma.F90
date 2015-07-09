!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_gensigma
      Implicit None
Contains
!
!
      Subroutine gensigma (w, eps, oc, sigma)
         Use mod_constants, Only: pi, zi
         Use modxs
         Implicit None
    ! arguments
         Real (8), Intent (In) :: w (:)
         Complex (8), Intent (In) :: eps (:)
         Integer, Intent (In) :: oc (2)
         Complex (8), Intent (Out) :: sigma (:)
    ! local variables
         Character (*), Parameter :: thisnam = 'gensigma'
         Real (8) :: delt
         If (any(shape(eps) .Ne. shape(sigma))) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): input and&
           & output arrays have diffenrent shape'
            Call terminate
         End If
    ! optical conductivity
         delt = 0.d0
         If (oc(1) .Eq. oc(2)) delt = 1.d0
!         sigma (:) = aimag (eps(:)) * w (:) / (4.d0*pi)
!         sigma (:) = sigma (:) + zi * &
!        & (-(dble(eps(:))-delt)*w(:)/(4.d0*pi))

        sigma (:) = - zi * (eps(:) - delt) * w (:) / (4.d0*pi)

      End Subroutine gensigma
!
End Module m_gensigma
