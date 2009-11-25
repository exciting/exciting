!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: symrvfir
!
!
Subroutine symrvfir (ngv, rvfir)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngv   : number of G-vectors to be used for the Fourier space rotation
!           (in,integer)
!   rvfir : real interstitial vector function (inout,real(ngrtot,ndmag))
! !DESCRIPTION:
!   Symmetrises a real interstitial vector function. See routines {\tt symrvf}
!   and {\tt symrfir} for details.
!
! !REVISION HISTORY:
!   Created July 2007 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: ngv
      Real (8), Intent (Inout) :: rvfir (ngrtot, ndmag)
! local variables
      Integer :: isym, lspl, ilspl, lspn
      Integer :: i, md, sym (3, 3), iv (3)
      Integer :: ig, ifg, jg, jfg
      Real (8) :: sc (3, 3), vtc (3), t1
      Complex (8) zv (3), zt1
! allocatable arrays
      Complex (8), Allocatable :: zfft1 (:, :), zfft2 (:, :)
      Allocate (zfft1(ngrtot, ndmag), zfft2(ngrtot, ndmag))
! Fourier transform vector function to G-space
      Do i = 1, ndmag
         zfft1 (:, i) = rvfir (:, i)
         Call zfftifc (3, ngrid,-1, zfft1(:, i))
      End Do
      zfft2 (:, :) = 0.d0
      Do isym = 1, nsymcrys
! translation vector in Cartesian coordinates
         Call r3mv (input%structure%crystal%basevect, vtlsymc(:, isym), &
        & vtc)
! index to spatial rotation lattice symmetry
         lspl = lsplsymc (isym)
! inverse rotation required for rotation of G-vectors
         ilspl = isymlat (lspl)
         sym (:, :) = symlat (:, :, ilspl)
! global spin proper rotation in Cartesian coordinates
         lspn = lspnsymc (isym)
         md = symlatd (lspn)
         sc (:, :) = dble (md) * symlatc (:, :, lspn)
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
! translation, spatial rotation and global spin rotation
            If (lspn .Eq. 1) Then
! global spin symmetry is the identity
               zfft2 (jfg, :) = zfft2 (jfg, :) + zt1 * zfft1 (ifg, :)
            Else
               If (ncmag) Then
! non-collinear case
                  zv (1) = sc (1, 1) * zfft1 (ifg, 1) + sc (1, 2) * &
                 & zfft1 (ifg, 2) + sc (1, 3) * zfft1 (ifg, 3)
                  zv (2) = sc (2, 1) * zfft1 (ifg, 1) + sc (2, 2) * &
                 & zfft1 (ifg, 2) + sc (2, 3) * zfft1 (ifg, 3)
                  zv (3) = sc (3, 1) * zfft1 (ifg, 1) + sc (3, 2) * &
                 & zfft1 (ifg, 2) + sc (3, 3) * zfft1 (ifg, 3)
                  zfft2 (jfg, :) = zfft2 (jfg, :) + zt1 * zv (:)
               Else
! collinear case
                  zfft2 (jfg, 1) = zfft2 (jfg, 1) + sc (3, 3) * zt1 * &
                 & zfft1 (ifg, 1)
               End If
            End If
         End Do
      End Do
! Fourier transform to real-space and normalise
      t1 = 1.d0 / dble (nsymcrys)
      Do i = 1, ndmag
         Call zfftifc (3, ngrid, 1, zfft2(:, i))
         rvfir (:, i) = t1 * dble (zfft2(:, i))
      End Do
      Deallocate (zfft1, zfft2)
      Return
End Subroutine
!EOC
