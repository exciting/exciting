!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_zaxpyc
      Implicit None
Contains
!
!
      Subroutine zaxpyc (n, za, zx, incx, zy, incy)
         Complex (8), Intent (In) :: zx (:), za
         Complex (8), Intent (Inout) :: zy (:)
         Integer, Intent (In) :: incx, incy, n
    ! automatic arrays
         Complex (8) :: zt (size(zx))
         zt (:) = conjg (zx(:))
         Call zaxpy (n, za, zt, incx, zy, incy)
      End Subroutine zaxpyc
!
End Module m_zaxpyc
