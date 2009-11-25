!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_genloss
      Implicit None
Contains
!
!
      Subroutine genloss (eps, loss)
         Use modxs
         Implicit None
    ! arguments
         Complex (8), Intent (In) :: eps (:)
         Real (8), Intent (Out) :: loss (:)
    ! local variables
         Character (*), Parameter :: thisnam = 'genloss'
         If (any(shape(eps) .Ne. shape(loss))) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): input and&
           & output arrays have diffenrent shape'
            Call terminate
         End If
    ! loss function
         loss (:) = - aimag (1/eps(:))
      End Subroutine genloss
!
End Module m_genloss
