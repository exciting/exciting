!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_gensumrls
      Implicit None
Contains
!
!BOP
! !ROUTINE: gensumrls
! !INTERFACE:
!
!
      Subroutine gensumrls (w, eps, sumrls)
! !USES:
         Use modxs
         Use modmpi
! !INPUT/OUTPUT PARAMETERS:
!   w      : frequency grid (in,real(:))
!   eps    : dielectric function tensor component (in,complex(:))
!   sumrls : values of three different sumrules (out,real(3))
! !DESCRIPTION:
!   Expressions for the sumrules taken from C. Ambrosch-Draxl,
!   CPC 175 (2006) 1-14, p5, Eq. 26.
!
! !REVISION HISTORY:
!   Created March 2006 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! arguments
         Real (8), Intent (In) :: w (:)
         Complex (8), Intent (In) :: eps (:)
         Real (8), Intent (Out) :: sumrls (3)
    ! local variables
         Character (*), Parameter :: thisnam = 'gensumrls'
         Real (8), Allocatable :: f (:), cf (:, :), g (:)
         Integer :: n1 (1), n
!
         If (any(shape(w) .Ne. shape(eps))) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): input arr&
           &ays have diffenrent shape'
            Call terminate
         End If
!
         n1 = shape (w)
         n = n1 (1)
         Allocate (f(n), g(n), cf(3, n))
!
    ! zeroth frequency moment sumrule
         f (:) = aimag (eps(:))
         Call fderiv (-1, n, w, f, g, cf)
         sumrls (1) = g (n)
!
   ! ! first frequency moment sumrule
   !      f (:) = aimag (-1/eps(:)) * w (:)
   !      Call fderiv (-1, n, w, f, g, cf)
   !      sumrls (2) = g (n)
         sumrls(2) = 0
!
    ! one over frequency sumrule (pi half sumrule)
         f (1) = 0.d0
         If (n .Gt. 1) f (2:) = aimag (-1/eps(2:)) / w (2:)
         Call fderiv (-1, n, w, f, g, cf)
         sumrls (3) = g (n)
!
         Deallocate (f, g, cf)
!
      End Subroutine gensumrls
!EOC
!
End Module m_gensumrls
