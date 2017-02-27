!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: kramkron
! !INTERFACE:
!
!
Subroutine kramkron (i1, i2, eps, n, w, im, re)
! !INPUT/OUTPUT PARAMETERS:
!   i1,i2 : tensor components of dielectric function tensor (in,integer)
!   eps   : zero frequency tolerance
!   n     : number of frequency points (in,integer)
!   w     : frequency grid (in,real(n))
!   im    : imaginary part of dielectric function tensor component (in,real(n))
!   re    : real part of dielectric function tensor component (out,real(n))
! !DESCRIPTION:
!   Performs a Kramers-Kronig transformation from the imaginary part to the
!   real part of the dielectric function.
!   Algorithm taken from routine {\tt linopt}.
!
! !REVISION HISTORY:
!   Created November 2007 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Integer, Intent (In) :: i1, i2
      Real (8), Intent (In) :: eps
      Integer, Intent (In) :: n
      Real (8), Intent (In) :: w (n), im (n)
      Real (8), Intent (Out) :: re (n)
  ! local variables
      Real (8), Parameter :: pi = 3.1415926535897932385d0
      Integer :: iw, jw
      Real (8) :: t1, t2
      Real (8), Allocatable :: fw (:), g (:), cf (:, :)
      Allocate (fw(n), g(n), cf(3, n))
      t1 = 0.d0
      g (:) = 0.d0
      If (i1 .Eq. i2) t1 = 1.d0
      Do iw = 1, n
         Do jw = 1, n
            t2 = w (jw) ** 2 - w (iw) ** 2
            !write(*,*) "iw, jw, t2", iw, jw, t2
        ! tolerance for range of integrand part w'/(w^2-w'^2)
            If (Abs(t2) .Gt. eps) Then
               fw (jw) = w (jw) * im (jw) / t2
               !write(*,*) "fw(jw)", fw(jw)
            Else
               fw (jw) = 0.d0
               !write(*,*) "fw(jw)", fw(jw)
            End If
         End Do
         Call fderiv (-1, n, w, fw, g, cf)
         !write(*,*) "g(n)", g(n)
         re (iw) = t1 + (2.d0/pi) * g (n)
      End Do
      Deallocate (fw, g, cf)
End Subroutine kramkron
