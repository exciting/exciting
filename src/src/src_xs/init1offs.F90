!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine init1offs (vploff)
      Use modmain, Only:
      Use modinput
      Use modxs, Only: init1norealloc
      Implicit None
  ! arguments
      Real (8), Intent (In) :: vploff (3)
  ! local variables
      Real (8) :: vt (3)
      Logical :: lt
      lt = init1norealloc
      init1norealloc = .True.
      vt (:) = input%groundstate%vkloff(:)
      input%groundstate%vkloff (:) = vploff (:)
  ! call init1 without selected re-allocations and specified k-point offset
      Call init1
      init1norealloc = lt
      input%groundstate%vkloff (:) = vt (:)
End Subroutine init1offs
