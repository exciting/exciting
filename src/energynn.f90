!
!
!
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine energynn
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: is, ia, ias, lmax
      Real (8) :: vn, t1
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
! set the density to zero
      zrhomt (:, :, :) = 0.d0
      zrhoir (:) = 0.d0
! compute the required spherical Bessel functions
      lmax = input%groundstate%lmaxvr + input%groundstate%npsden + 1
      Call genjlgpr (lmax, gc, jlgr)
! solve the complex Poisson's equation
      Call zpotcoul (nrmt, nrmtmax, spnrmax, spr, 1, gc, jlgr, ylmg, &
     & sfacg, spzn, zrhomt, zrhoir, zvclmt, zvclir, zrho0)
! compute the nuclear-nuclear energy
      engynn = 0.d0
      Do is = 1, nspecies
! compute the bare nucleus potential at the origin
         Call potnucl (input%groundstate%ptnucl, 1, spr(:, is), &
        & spzn(is), vn)
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            t1 = dble (zvclmt(1, 1, ias)) * y00 - vn
            engynn = engynn + spzn (is) * t1
         End Do
      End Do
      engynn = 0.5d0 * engynn
      Deallocate (jlgr, zrhomt, zrhoir, zvclmt, zvclir)
      Return
End Subroutine
