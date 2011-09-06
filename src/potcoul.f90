!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: potcoul
! !INTERFACE:
!
!
Subroutine potcoul
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Calculates the Coulomb potential of the real charge density stored in the
!   global variables {\tt rhomt} and {\tt rhoir} by solving Poisson's equation.
!   These variables are converted to complex representations and passed to the
!   routine {\tt zpotcoul}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ir, lmax
      Complex (8) zrho0
! allocatable arrays
      Real (8), Allocatable :: jlgr (:, :, :)
      Complex (8), Allocatable :: zrhomt (:, :, :)
      Complex (8), Allocatable :: zrhoir (:)
      Complex (8), Allocatable :: zvclmt (:, :, :)
      Complex (8), Allocatable :: zvclir (:)
      Allocate &
     & (jlgr(0:input%groundstate%lmaxvr+input%groundstate%npsden+1, &
     & ngvec, nspecies))
      Allocate (zrhomt(lmmaxvr, nrmtmax, natmtot))
      Allocate (zrhoir(ngrtot))
      Allocate (zvclmt(lmmaxvr, nrmtmax, natmtot))
      Allocate (zvclir(ngrtot))
! convert real muffin-tin charge density to complex spherical harmonic expansion
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               Call rtozflm (input%groundstate%lmaxvr, rhomt(:, ir, &
              & ias), zrhomt(:, ir, ias))
            End Do
         End Do
      End Do
! store real interstitial charge density in complex array
      zrhoir (:) = rhoir (:)
! compute the required spherical Bessel functions
      lmax = input%groundstate%lmaxvr + input%groundstate%npsden + 1
      Call genjlgpr (lmax, gc, jlgr)
! solve the complex Poisson's equation
      Call zpotcoul (nrmt, nrmtmax, spnrmax, spr, 1, gc, jlgr, ylmg, &
     & sfacg, spzn, zrhomt, zrhoir, zvclmt, zvclir, zrho0)
! convert complex muffin-tin potential to real spherical harmonic expansion
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is)
               Call ztorflm (input%groundstate%lmaxvr, zvclmt(:, ir, &
              & ias), vclmt(:, ir, ias))
            End Do
         End Do
      End Do
! store complex interstitial potential in real array
      vclir (:) = dble (zvclir(:))
      Deallocate (jlgr, zrhomt, zrhoir, zvclmt, zvclir)
      Return
End Subroutine
!EOC
