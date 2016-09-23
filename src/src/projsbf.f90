!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine projsbf
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: is, ia, ias, ir
      Integer :: idm, lmax, lm
      Real (8) :: t1
      Complex (8) zrho0
! automatic arrays
      Real (8) :: zn (nspecies)
! allocatable arrays
      Real (8), Allocatable :: rvfmt (:, :, :, :)
      Real (8), Allocatable :: rvfir (:, :)
      Real (8), Allocatable :: rfmt (:, :, :)
      Real (8), Allocatable :: rfir (:)
      Real (8), Allocatable :: grfmt (:, :, :, :)
      Real (8), Allocatable :: grfir (:, :)
      Real (8), Allocatable :: jlgr (:, :, :)
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zvclmt (:, :, :)
      Complex (8), Allocatable :: zvclir (:)
      Allocate (rvfmt(lmmaxvr, nrmtmax, natmtot, 3))
      Allocate (rvfir(ngrtot, 3))
      Allocate (rfmt(lmmaxvr, nrmtmax, natmtot))
      Allocate (rfir(ngrtot))
      Allocate (grfmt(lmmaxvr, nrmtmax, natmtot, 3))
      Allocate (grfir(ngrtot, 3))
      Allocate &
     & (jlgr(0:input%groundstate%lmaxvr+input%groundstate%npsden+1, &
     & ngvec, nspecies))
      Allocate (zrhomt(lmmaxvr, nrmtmax, natmtot))
      Allocate (zrhoir(ngrtot))
      Allocate (zvclmt(lmmaxvr, nrmtmax, natmtot))
      Allocate (zvclir(ngrtot))
      If ( .Not. associated(input%groundstate%spin)) Then
         Write (*,*)
         Write (*, '("Error(projsbf): spin-unpolarised field is zero")')
         Write (*,*)
         Stop
      End If
      If (ncmag) Then
! non-collinear
         rvfmt (:, :, :, :) = bxcmt (:, :, :, :)
         rvfir (:, :) = bxcir (:, :)
      Else
! collinear
         rvfmt (:, :, :, 1:2) = 0.d0
         rvfir (:, 1:2) = 0.d0
         rvfmt (:, :, :, 3) = bxcmt (:, :, :, 1)
         rvfir (:, 3) = bxcir (:, 1)
      End If
! compute the divergence of B-field
      rfmt (:, :, :) = 0.d0
      rfir (:) = 0.d0
      Do idm = 1, 3
         Call gradrf (rvfmt(:, :, :, idm), rvfir(:, idm), grfmt, grfir)
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ir = 1, nrmt (is)
                  rfmt (:, ir, ias) = rfmt (:, ir, ias) + grfmt (:, ir, &
                 & ias, idm)
               End Do
            End Do
         End Do
         rfir (:) = rfir (:) + grfir (:, idm)
      End Do
! divide by -4*pi
      t1 = - 1.d0 / fourpi
      rfmt (:, :, :) = t1 * rfmt (:, :, :)
      rfir (:) = t1 * rfir (:)
! convert real muffin-tin divergence to complex spherical harmonic expansion
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               Call rtozflm (input%groundstate%lmaxvr, rfmt(:, ir, &
              & ias), zrhomt(:, ir, ias))
            End Do
         End Do
      End Do
! store real interstitial divergence in a complex array
      zrhoir (:) = rfir (:)
! set the point charges to zero
      zn (:) = 0.d0
! compute the required spherical Bessel functions
      lmax = input%groundstate%lmaxvr + input%groundstate%npsden + 1
      Call genjlgpr (lmax, gc, jlgr)
! solve the complex Poisson's equation
      Call zpotcoul (nrmt, nrmtmax, spnrmax, spr, 1, gc, jlgr, ylmg, &
     & sfacg, zn, zrhomt, zrhoir, zvclmt, zvclir, zrho0)
! convert complex muffin-tin potential to real spherical harmonic expansion
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               Call ztorflm (input%groundstate%lmaxvr, zvclmt(:, ir, &
              & ias), rfmt(:, ir, ias))
            End Do
         End Do
      End Do
! store complex interstitial potential in real array
      rfir (:) = dble (zvclir(:))
! compute the gradient
      Call gradrf (rfmt, rfir, grfmt, grfir)
! subtract gradient from existing B-field
      If (ncmag) Then
! non-collinear
         bxcmt (:, :, :, :) = bxcmt (:, :, :, :) - grfmt (:, :, :, :)
         bxcir (:, :) = bxcir (:, :) - grfir (:, :)
      Else
! collinear
         bxcmt (:, :, :, 1) = bxcmt (:, :, :, 1) - grfmt (:, :, :, 3)
         bxcir (:, 1) = bxcir (:, 1) - grfir (:, 3)
      End If
! remove numerical noise from the muffin-tin B-field
      Do idm = 1, ndmag
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do lm = 1, lmmaxvr
                  Call fsmooth (10, nrmt(is), lmmaxvr, bxcmt(lm, 1, &
                 & ias, idm))
               End Do
            End Do
         End Do
      End Do
      Deallocate (rvfmt, rvfir, rfmt, rfir, grfmt, grfir, jlgr)
      Deallocate (zrhomt, zrhoir, zvclmt, zvclir)
      Return
End Subroutine
