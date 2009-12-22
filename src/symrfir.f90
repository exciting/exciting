!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: symrfir
!
!
Subroutine symrfir (ngv, rfir)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngv  : number of G-vectors to be used for the Fourier space rotation
!          (in,integer)
!   rfir : real intersitial function (inout,real(ngrtot))
! !DESCRIPTION:
!   Symmetrises a real scalar interstitial function. The function is first
!   Fourier transformed to $G$-space, and then averaged over each symmetry by
!   rotating the Fourier coefficients and multiplying them by a phase factor
!   corresponding to the symmetry translation.
!
! !REVISION HISTORY:
!   Created July 2007 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: ngv
      Real (8), Intent (Inout) :: rfir (ngrtot)
! local variables
      Integer :: isym, lspl, ilspl, sym (3, 3)
      Integer :: iv (3), ig, jg, ifg, jfg
      Real (8) :: vtc (3), t1
      Complex (8) zt1
! allocatable arrays
      Complex (8), Allocatable :: zfft1 (:), zfft2 (:)
      Allocate (zfft1(ngrtot), zfft2(ngrtot))
! Fourier transform function to G-space
      zfft1 (:) = rfir (:)
      Call zfftifc (3, ngrid,-1, zfft1)
      zfft2 (:) = 0.d0
! loop over crystal symmetries
      Do isym = 1, nsymcrys
! translation in Cartesian coordinates
         Call r3mv (input%structure%crystal%basevect, vtlsymc(:, isym), &
        & vtc)
! index to lattice symmetry of spatial rotation
         lspl = lsplsymc (isym)
! inverse rotation required for rotation of G-vectors
         ilspl = isymlat (lspl)
         sym (:, :) = symlat (:, :, ilspl)
         Do ig = 1, ngv
            ifg = igfft (ig)
! multiply the transpose of the inverse symmetry matrix with the G-vector
            iv (1) = sym (1, 1) * ivg (1, ig) + sym (2, 1) * ivg (2, &
           & ig) + sym (3, 1) * ivg (3, ig)
            iv (2) = sym (1, 2) * ivg (1, ig) + sym (2, 2) * ivg (2, &
           & ig) + sym (3, 2) * ivg (3, ig)
            iv (3) = sym (1, 3) * ivg (1, ig) + sym (2, 3) * ivg (2, &
           & ig) + sym (3, 3) * ivg (3, ig)
            iv (:) = modulo (iv(:)-intgv(:, 1), ngrid(:)) + intgv (:, &
           & 1)
            jg = ivgig (iv(1), iv(2), iv(3))
            jfg = igfft (jg)
! complex phase factor for translation
            t1 = - dot_product (vgc(:, ig), vtc(:))
            zt1 = cmplx (Cos(t1), Sin(t1), 8)
            zfft2 (jfg) = zfft2 (jfg) + zt1 * zfft1 (ifg)
         End Do
      End Do
! Fourier transform to real-space and normalise
      Call zfftifc (3, ngrid, 1, zfft2)
      t1 = 1.d0 / dble (nsymcrys)
      rfir (:) = t1 * dble (zfft2(:))
      Deallocate (zfft1, zfft2)
      Return
End Subroutine
!EOC
