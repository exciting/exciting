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
      Subroutine genloss (eps, loss, tq0)
         Use modxs
         use modinput
         Implicit None
    ! arguments
         Complex (8), Intent (In) :: eps (:, :, :)
         Real (8), Intent (Out) :: loss (:, :, :)
         Logical :: tq0
    ! local variables
         integer :: iw
         complex(8) :: t3(3, 3)
         Character (*), Parameter :: thisnam = 'genloss'

    !
         If (any(shape(eps) .Ne. shape(loss))) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): input and&
           & output arrays have diffenrent shape'
            Call terminate
         End If
    ! loss function
! STK
         If(tq0) Then
            Do iw = 1, input%xs%energywindow%points
               call z3minv(eps(:, :, iw), t3(:, :))
               loss(:, :, iw) = -aimag(t3)
            enddo
         Else
            loss (1,1,:) = - aimag (1/eps(1,1,:))
         End If
      End Subroutine genloss
!
End Module m_genloss
