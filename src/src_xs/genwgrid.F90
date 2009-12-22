!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_genwgrid
      Implicit None
!
Contains
!
!
      Subroutine genwgrid (n, intv, timag, brd, w_real, w_cmplx)
         Implicit None
    ! arguments
         Integer, Intent (In) :: n
         Real (8), Intent (In) :: intv (2)
    ! optional arguments
         Logical, Optional, Intent (In) :: timag
         Real (8), Optional, Intent (In) :: brd
         Real (8), Optional, Intent (Out) :: w_real (:)
         Complex (8), Optional, Intent (Out) :: w_cmplx (:)
    ! local variables
         Character (*), Parameter :: thisnam = 'genwgrid'
         Integer :: j, nerr
         Real (8) :: brdt, t1
         Complex (8) :: fac
!
    ! checking
         nerr = 0
         If (present(w_real) .And. present(w_cmplx)) Then
            Write (*, '("Error(", a, "both real and complex output grid&
           & specified")') thisnam
            nerr = nerr + 1
         End If
         If (nerr .Gt. 0) Stop
!
         If (present(w_real)) Then
       ! real grid
            t1 = (intv(2)-intv(1)) / dble (n)
            Do j = 1, n
               w_real (j) = t1 * dble (j-1) + intv (1)
            End Do
         Else
       ! complex grid
            brdt = 0.d0
            fac = (1.d0, 0.d0)
            If (present(brd)) brdt = brd
            If (present(timag)) Then
               If (timag) fac = (0.d0, 1.d0)
            End If
            t1 = (intv(2)-intv(1)) / dble (n)
            Do j = 1, n
               w_cmplx (j) = fac * (t1*dble(j-1)+intv(1)) + (0.d0, &
              & 1.d0) * brdt
            End Do
         End If
!
      End Subroutine genwgrid
!
End Module m_genwgrid
