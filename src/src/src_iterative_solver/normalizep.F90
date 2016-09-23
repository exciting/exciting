
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine normalizep (n, m, overlap, evecfv, ldv)
      Use modmain, Only: zone, zzero
!
      Implicit None
      Integer, Intent (In) :: n, m, ldv
      Complex (8), Intent (In) :: overlap (n*(n+1)/2)
      Complex (8), Intent (Out) :: evecfv (ldv, m)
      Complex (8) :: tmp (n), z
      Complex (8), External :: zdotc
      Integer :: i
      Do i = 1, m
         Call zhpmv ('U', n, zone, overlap(1), evecfv(1, i), 1, zzero, &
        & tmp, 1)
         z = Sqrt (dble(zdotc(n, evecfv(1, i), 1, tmp, 1)))
         z = 1.0 / z
         Call zscal (n, z, evecfv(1, i), 1)
      End Do
End Subroutine
