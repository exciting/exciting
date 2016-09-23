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
      Subroutine genloss (eps, loss, nc)
         Use modxs
         use modinput
         Implicit None
    ! arguments
         Complex (8), Intent (In) :: eps (:, :, :)
         Real (8), Intent (Out) :: loss (:, :, :)
         Integer, Intent (In) :: nc
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
         If(nc.Eq.1) Then
            loss (1,1,:) = - aimag (1/eps(1,1,:))
         Else
            Do iw = 1, input%xs%energywindow%points
               call z3minv(eps(:, :, iw), t3(:, :))
               loss(:, :, iw) = -aimag(t3)
            enddo
         End If
      End Subroutine genloss
!
End Module m_genloss
