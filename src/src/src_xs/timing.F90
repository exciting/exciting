!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine xstiming
      Use invert
      Use blaswrappers
      Use m_getunit
      Implicit None
  ! local variables
      Integer :: un, nlo, nup, ninc, j, k, npts
      Real (8) :: cpu0, cpu1, ts0, ts1, tc_mm, tc_inv, tc_diag, ts_mm, &
     & ts_inv, ts_diag
      Real (8) :: f
      Complex (8) :: alpha, beta
      Real (8), Allocatable :: v (:), mr (:, :), tc (:), x (:), c (:)
      Complex (8), Allocatable :: mz (:, :), mh (:, :), ma (:, :), mb &
     & (:, :)
      Real (8), External :: polynom
      nlo = 100
      nup = 2000
      ninc = 100
      alpha = (2.d0, 3.d0)
      beta = (4.d0, 5.d0)
      Call getunit (un)
      Open (un, File='TIMING.OUT', Form='formatted', Action='write', &
     & Status='replace')
      npts = (nup-nlo) / ninc + 1
      Allocate (x(npts), tc(npts), c(npts))
      k = 0
      Do j = nlo, nup, ninc
         k = k + 1
         Allocate (mr(j, j), v(j), mz(j, j), mh(j, j), ma(j, j), mb(j, &
        & j))
     ! set up random number matrix
         mz (:, :) = 0.d0
         Call random_number (mr)
         mz = mz + mr
         Call random_number (mr)
         mz = mz + (0.d0, 1.d0) * mr
         ma = mz
     ! test matrix multiplication
         Call timesec (ts0)
         Call cpu_time (cpu0)
         Call zgemm_wrap (mb, 'n', ma, 'n', mz, alpha, beta)
         Call timesec (ts1)
         Call cpu_time (cpu1)
         ts_mm = ts1 - ts0
         tc_mm = cpu1 - cpu0
     ! test inversion
         Call timesec (ts0)
         Call cpu_time (cpu0)
         Call zinvert_lapack (mz, mb)
         Call timesec (ts1)
         Call cpu_time (cpu1)
         ts_inv = ts1 - ts0
         tc_inv = cpu1 - cpu0
     ! test diagonalization
         Call timesec (ts0)
         Call cpu_time (cpu0)
         Call bsesoldiag (j, j, mz, v, mb)
         Call timesec (ts1)
         Call cpu_time (cpu1)
         ts_diag = ts1 - ts0
         tc_diag = cpu1 - cpu0
         Write (un, '(i8, 3g18.10)') j, tc_mm, tc_inv, tc_diag
         Deallocate (mr, v, mz, mh, ma, mb)
         x (k) = j
         tc (k) = tc_diag
      End Do
      Close (un)
      f = polynom (0, 5, x, tc, c, 1000.d0)
      Write (*,*) 'polynomial fit: coefficients a0, a1, a2, a3', c
      Write (*,*) 'extrapolation to 1000:', f
End Subroutine xstiming
