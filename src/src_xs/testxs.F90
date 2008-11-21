
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine testxs
  implicit none
  call test_invert
end subroutine testxs

subroutine test_genscclieff
  use modmain
  implicit none
  call init0
  call init1
  call genscclieff
end subroutine test_genscclieff

subroutine test_invert
  use invert
  implicit none
  integer :: siz
  real(8) :: ts0,ts1
  real(8), allocatable :: r(:,:)
  complex(8), allocatable :: a(:,:), b(:,:)
  write(*,*) 'give matrix size:'
  read(*,*) siz
  allocate(r(siz,siz),a(siz,siz),b(siz,siz))
  a=0.d0
  b=0.d0
  call random_number(r)
  a=a+r
  call random_number(r)
  a=a+(0.d0,1.d0)*r
  call timesec(ts0)
  call zinvert_lapack(a,b)
  call timesec(ts1)
  write(*,*) 'time for inversion of complex matrix using zgetrf/zgetri:',ts1-ts0
  deallocate(r,a,b)
end subroutine test_invert
