!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine getevalsv0 (vpl, evalsvp)
      Use modmain
      Use modxs
  ! arguments
      Real (8), Intent (In) :: vpl (3)
      Real (8), Intent (Out) :: evalsvp (nstsv)
  ! local variables
      Real (8), Allocatable :: vklt (:, :)
      Character (256) :: filextt
  ! copy varialbes of k+(q=0) to default variables
      Allocate (vklt(3, nkptnr))
      vklt (:, :) = vkl (:, :)
      vkl (:, :) = vkl0 (:, :)
      filextt = filext
  ! call to getevalsv with changed (G+)k-point sets / matrix size
      Call genfilextread (task)
      Call getevalsv (vpl, evalsvp)
  ! restore original variables
      vkl (:, :) = vklt (:, :)
      filext = filextt
      Deallocate (vklt)
End Subroutine getevalsv0
