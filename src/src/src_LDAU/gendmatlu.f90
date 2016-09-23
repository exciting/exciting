!
!
!
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine gendmatlu
      Use modmain
      Implicit None
! local variables
      Integer :: ik, ispn, ist
      Integer :: is, ia, ias
      Real (8) :: t1
! allocatable arrays
      Complex (8), Allocatable :: apwalm (:, :, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: dmat (:, :, :, :, :)
! allocate local arrays
      Allocate (evecsv(nstsv, nstsv))
      Allocate (evecfv(nmatmax, nstfv, nspnfv))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
      Allocate (dmat(lmmaxlu, lmmaxlu, nspinor, nspinor, nstsv))
! zero the LDA+U density matrix
      dmatlu (:, :, :, :, :) = 0.d0
! begin loop over k-points
      Do ik = 1, nkpt
! get the eigenvectors and occupancies from file
         Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
         Call getevecsv (vkl(:, ik), evecsv)
         Call getoccsv (vkl(:, ik), occsv(:, ik))
! find the matching coefficients
         Do ispn = 1, nspnfv
            Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, &
           & ispn, ik), sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, &
           & ispn))
         End Do
! begin loop over atoms and species
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Call gendmat (.False., .False., 0, lmaxlu, is, ia, &
              & ngk(:, ik), apwalm, evecfv, evecsv, lmmaxlu, dmat)
               Do ist = 1, nstsv
                  t1 = wkpt (ik) * occsv (ist, ik)
                  dmatlu (:, :, :, :, ias) = dmatlu (:, :, :, :, ias) + &
                 & t1 * dmat (:, :, :, :, ist)
               End Do
            End Do
         End Do
      End Do
! symmetrise the density matrix
      Call symdmat (lmaxlu, lmmaxlu, dmatlu)
      Deallocate (evecfv, evecsv, apwalm, dmat)
      Return
End Subroutine
