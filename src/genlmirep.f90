!
!
!
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genlmirep (lmax, ld, elm, ulm)
      Use modmain
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: ld
      Real (8), Intent (Out) :: elm (ld, natmtot)
      Complex (8), Intent (Out) :: ulm (ld, ld, natmtot)
! local variables
      Integer :: is, ia, ias
      Integer :: lmmax, i, j, l, lm, n
      Integer :: isym, lspl, info, lwork
! allocatable arrays
      Real (8), Allocatable :: rwork (:)
      Complex (8), Allocatable :: ulat (:, :, :)
      Complex (8), Allocatable :: a (:, :), b (:, :)
      Complex (8), Allocatable :: h (:, :)
      Complex (8), Allocatable :: work (:)
      lmmax = (lmax+1) ** 2
      Allocate (rwork(3*lmmax))
      Allocate (ulat(lmmax, lmmax, nsymlat))
      Allocate (a(lmmax, lmmax), b(lmmax, lmmax))
      Allocate (h(lmmax, lmmax))
      lwork = 2 * lmmax
      Allocate (work(lwork))
! construct (l,m) rotation matrix for each lattice symmetry
      a (:, :) = 0.d0
      Do i = 1, lmmax
         a (i, i) = 1.d0
      End Do
      Do isym = 1, nsymlat
         Call rotzflm (symlatc(:, :, isym), lmax, lmmax, lmmax, a, &
        & ulat(:, :, isym))
      End Do
! set up quasi-random symmetric matrix H
      h (:, :) = 0.d0
      Do l = 0, lmax
         n = 2 * l + 1
         lm = idxlm (l,-l)
         Do i = lm, lm + n - 1
            Do j = i, lm + n - 1
               h (i, j) = dble (i*j)
               h (j, i) = h (i, j)
            End Do
         End Do
      End Do
! loop over species and atoms
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! symmetrise H with site symmetries
            b (:, :) = 0.d0
            Do isym = 1, nsymsite (ias)
! spatial rotation element in lattice point group
               lspl = lsplsyms (isym, ias)
! apply lattice symmetry as U*H*conjg(U')
               Call zgemm ('N', 'N', lmmax, lmmax, lmmax, zone, ulat(:, &
              & :, lspl), lmmax, h, lmmax, zzero, a, lmmax)
               Call zgemm ('N', 'C', lmmax, lmmax, lmmax, zone, a, &
              & lmmax, ulat(:, :, lspl), lmmax, zone, b, lmmax)
            End Do
! block diagonalise symmetrised H
            Do l = 0, lmax
               n = 2 * l + 1
               lm = idxlm (l,-l)
               Call zheev ('V', 'U', n, b(lm, lm), lmmax, elm(lm, ias), &
              & work, lwork, rwork, info)
            End Do
! the unitary matrix U is the transpose of the eigenvector array
            Do i = 1, lmmax
               Do j = 1, lmmax
                  ulm (i, j, ias) = b (j, i)
               End Do
            End Do
         End Do
      End Do
      Deallocate (rwork, ulat, a, b, h, work)
      Return
End Subroutine
