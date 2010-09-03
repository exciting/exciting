
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine getevecfv0 (vpl, vgpl, evecfvt)
      Use modmain
      Use modxs
  ! arguments
      Real (8), Intent (In) :: vpl (3)
      Real (8), Intent (In) :: vgpl (3, ngkmax)
      Complex (8), Intent (Out) :: evecfvt (nmatmax, nstfv, nspnfv)
  ! local variables
      Integer :: nmatmaxt, ngkmaxt
      Integer, Allocatable :: ngkt (:, :)
      Real (8), Allocatable :: vklt (:, :), vgklt (:, :, :, :)
      Character (256) :: filextt
  ! copy varialbes of k+(q=0) to default variables
      Allocate (ngkt(nspnfv, nkpt))
      Allocate (vklt(3, nkptnr))
      Allocate (vgklt(3, ngkmax, nspnfv, nkpt))
      nmatmaxt = nmatmax
      nmatmax = nmatmax0
      ngkmaxt = ngkmax
      ngkmax = ngkmax0
      ngkt (:, :) = ngk (:, :)
      ngk (:, :) = ngk0 (:, :)
      vklt (:, :) = vkl (:, :)
      vkl (:, :) = vkl0 (:, :)
      vgklt (:, :, :, :) = vgkl (:, :, :, :)
      ! re-allocate array
      deallocate(vgkl)
      allocate(vgkl(3, ngkmax0, nspnfv, nkpt))
      vgkl (:, :, :, :) = vgkl0 (:, :, :, :)
      filextt = filext
  ! call to getevecfv with changed (G+)k-point sets / matrix size
      Call genfilextread (task)
      Call getevecfv (vpl, vgpl, evecfvt)
  ! restore original variables
      nmatmax = nmatmaxt
      ngkmax = ngkmaxt
      ngk (:, :) = ngkt (:, :)
      vkl (:, :) = vklt (:, :)
      ! re-allocate array
      deallocate(vgkl)
      allocate(vgkl(3, ngkmax, nspnfv, nkpt))
      vgkl (:, :, :, :) = vgklt (:, :, :, :)
      filext = filextt
      Deallocate (ngkt, vklt, vgklt)
End Subroutine getevecfv0
