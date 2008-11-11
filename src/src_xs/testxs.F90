
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine testxs
  implicit none
  call test_iplocnr
end subroutine testxs

subroutine test_genscclieff
  use modmain
  implicit none
  call init0
  call init1
  call genscclieff
end subroutine test_genscclieff


subroutine test_iplocnr
  use modmain
  use modxs
  implicit none
  integer, external :: iplocnr
  task=440
  call init0
  call init1
  call init2
  write(*,*) iplocnr((/0,1,5/),(/2,2,16/))
  stop 'end of test_iplocnr'
end subroutine test_iplocnr

