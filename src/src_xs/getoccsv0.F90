!
!
!
! Copyright (C) 2007-2008 S. Sagmeister J. K. Dewhurst, S. Sharma and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine getoccsv0 (vpl, occsvp)
      Use modmain
      Use modxs
  ! arguments
      Real (8), Intent (In) :: vpl (3)
      Real (8), Intent (Out) :: occsvp (nstsv)
  ! local variables
      Real (8), Allocatable :: vklt (:, :)
      Character (256) :: filextt
!
  ! copy varialbes of k+(q=0) to default variables
      Allocate (vklt(3, nkptnr))
      vklt (:, :) = vkl (:, :)
      vkl (:, :) = vkl0 (:, :)
      filextt = filext
!
  ! call to getevalsv with changed (G+)k-point sets / matrix size
      Call genfilextread (task)
      Call getoccsv (vpl, occsvp)
!
  ! restore original variables
      vkl (:, :) = vklt (:, :)
      filext = filextt
      Deallocate (vklt)
!
End Subroutine getoccsv0
