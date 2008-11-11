
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine testxs
  implicit none
  call test_genscclieff
end subroutine testxs

subroutine test_genscclieff
  use modmain
  implicit none
  call init0
  call init1
  call genscclieff
end subroutine test_genscclieff
