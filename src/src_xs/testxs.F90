!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine testxs
      Implicit None
      Call test_invert
End Subroutine testxs
!
!
Subroutine test_genscclieff
      Use modmain
      Implicit None
      Call init0
      Call init1
      Call genscclieff
End Subroutine test_genscclieff
!
!
Subroutine test_invert
      Use invert
      Implicit None
      Integer :: siz
      Real (8) :: ts0, ts1
      Real (8), Allocatable :: r (:, :)
      Complex (8), Allocatable :: a (:, :), b (:, :)
      Write (*,*) 'give matrix size:'
      Read (*,*) siz
      Allocate (r(siz, siz), a(siz, siz), b(siz, siz))
      a = 0.d0
      b = 0.d0
      Call random_number (r)
      a = a + r
      Call random_number (r)
      a = a + (0.d0, 1.d0) * r
      Call timesec (ts0)
      Call zinvert_lapack (a, b)
      Call timesec (ts1)
      Write (*,*) 'time for inversion of complex matrix using zgetrf/zg&
     &etri:', ts1 - ts0
      Deallocate (r, a, b)
End Subroutine test_invert
