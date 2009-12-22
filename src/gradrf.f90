!
!
!
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine gradrf (rfmt, rfir, grfmt, grfir)
      Use modmain
      Use modinput
      Implicit None
      Real (8), Intent (In) :: rfmt (lmmaxvr, nrmtmax, natmtot)
      Real (8), Intent (In) :: rfir (ngrtot)
      Real (8), Intent (Out) :: grfmt (lmmaxvr, nrmtmax, natmtot, 3)
      Real (8), Intent (Out) :: grfir (ngrtot, 3)
! local variables
      Integer :: is, ia, ias, i, ig, ifg
! allocatable arrays
      Real (8), Allocatable :: grfmt1 (:, :, :)
      Complex (8), Allocatable :: zfft1 (:), zfft2 (:)
      Allocate (grfmt1(lmmaxvr, nrmtmax, 3))
      Allocate (zfft1(ngrtot), zfft2(ngrtot))
! muffin-tin gradient
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Call gradrfmt (input%groundstate%lmaxvr, nrmt(is), spr(:, &
           & is), lmmaxvr, nrmtmax, rfmt(:, :, ias), grfmt1)
            Do i = 1, 3
               grfmt (:, 1:nrmt(is), ias, i) = grfmt1 (:, 1:nrmt(is), &
              & i)
            End Do
         End Do
      End Do
! interstitial gradient
      zfft1 (:) = rfir (:)
      Call zfftifc (3, ngrid,-1, zfft1)
      Do i = 1, 3
         zfft2 (:) = 0.d0
         Do ig = 1, ngvec
            ifg = igfft (ig)
            zfft2 (ifg) = zi * vgc (i, ig) * zfft1 (ifg)
         End Do
         Call zfftifc (3, ngrid, 1, zfft2)
         grfir (:, i) = dble (zfft2(:))
      End Do
      Deallocate (grfmt1, zfft1, zfft2)
      Return
End Subroutine
