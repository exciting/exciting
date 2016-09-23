!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: symrvf
! !INTERFACE:
!
!
Subroutine symrvf (lrstp, rvfmt, rvfir)
! !USES:
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   rvfmt : real muffin-tin vector field
!           (in,real(lmmaxvr,nrmtmax,natmtot,ndmag))
!   rvfir : real interstitial vector field
!           (in,real(ngrtot,ndmag))
! !DESCRIPTION:
!   Symmetrises a vector field defined over the entire unit cell using the full
!   set of crystal symmetries. If a particular symmetry involves rotating atom
!   1 into atom 2, then the spatial and spin rotations of that symmetry are
!   applied to the vector field in atom 2 (expressed in spherical harmonic
!   coefficients), which is then added to the field in atom 1. This is repeated
!   for all symmetry operations. The fully symmetrised field in atom 1 is then
!   rotated and copied to atom 2. Symmetrisation of the interstitial part of the
!   field is performed by {\tt symrvfir}. See also {\tt symrfmt} and
!   {\tt findsym}.
!
! !REVISION HISTORY:
!   Created May 2007 (JKD)
!   Fixed problem with improper rotations, February 2008 (L. Nordstrom,
!    F. Bultmark and F. Cricchio)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lrstp
      Real (8), Intent (Inout) :: rvfmt (lmmaxvr, nrmtmax, natmtot, &
     & ndmag)
      Real (8), Intent (Inout) :: rvfir (ngrtot, ndmag)
! local variables
      Integer :: is, ia, ja, ias, jas
      Integer :: isym, ir, lm, i, md
      Integer :: lspl, ilspl, lspn, ilspn
      Real (8) :: sc (3, 3), v (3), t1
! automatic arrays
      Logical :: done (natmmax)
! allocatable arrays
      Real (8), Allocatable :: rvfmt1 (:, :, :, :), rvfmt2 (:, :, :)
      Allocate (rvfmt1(lmmaxvr, nrmtmax, natmmax, ndmag))
      Allocate (rvfmt2(lmmaxvr, nrmtmax, ndmag))
      t1 = 1.d0 / dble (nsymcrys)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
      Do is = 1, nspecies
! make copy of vector field for all atoms of current species
         Do i = 1, ndmag
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Do ir = 1, nrmt (is), lrstp
                  rvfmt1 (:, ir, ia, i) = rvfmt (:, ir, ias, i)
               End Do
            End Do
         End Do
         done (:) = .False.
         Do ia = 1, natoms (is)
            If ( .Not. done(ia)) Then
               ias = idxas (ia, is)
               Do ir = 1, nrmt (is), lrstp
                  rvfmt (:, ir, ias, :) = 0.d0
               End Do
! begin loop over crystal symmetries
               Do isym = 1, nsymcrys
! equivalent atom
                  ja = ieqatom (ia, is, isym)
! parallel transport of vector field
                  lspl = lsplsymc (isym)
                  Do i = 1, ndmag
                     Call symrfmt (lrstp, is, symlatc(:, :, lspl), &
                    & rvfmt1(:, :, ja, i), rvfmt2(:, :, i))
                  End Do
! global spin proper rotation matrix in Cartesian coordinates
                  lspn = lspnsymc (isym)
                  md = symlatd (lspn)
                  sc (:, :) = dble (md) * symlatc (:, :, lspn)
! global spin rotation of vector field
                  If (ncmag) Then
! non-collinear case
                     Do ir = 1, nrmt (is), lrstp
                        Do lm = 1, lmmaxvr
                           v (1) = sc (1, 1) * rvfmt2 (lm, ir, 1) + sc &
                          & (1, 2) * rvfmt2 (lm, ir, 2) + sc (1, 3) * &
                          & rvfmt2 (lm, ir, 3)
                           v (2) = sc (2, 1) * rvfmt2 (lm, ir, 1) + sc &
                          & (2, 2) * rvfmt2 (lm, ir, 2) + sc (2, 3) * &
                          & rvfmt2 (lm, ir, 3)
                           v (3) = sc (3, 1) * rvfmt2 (lm, ir, 1) + sc &
                          & (3, 2) * rvfmt2 (lm, ir, 2) + sc (3, 3) * &
                          & rvfmt2 (lm, ir, 3)
                           rvfmt (lm, ir, ias, :) = rvfmt (lm, ir, ias, &
                          & :) + v (:)
                        End Do
                     End Do
                  Else
! collinear case
                     Do ir = 1, nrmt (is), lrstp
                        Do lm = 1, lmmaxvr
                           rvfmt (lm, ir, ias, 1) = rvfmt (lm, ir, ias, &
                          & 1) + sc (3, 3) * rvfmt2 (lm, ir, 1)
                        End Do
                     End Do
                  End If
! end loop over crystal symmetries
               End Do
! normalise
               Do ir = 1, nrmt (is), lrstp
                  rvfmt (:, ir, ias, :) = t1 * rvfmt (:, ir, ias, :)
               End Do
! mark atom as done
               done (ia) = .True.
! rotate into equivalent atoms
               Do isym = 1, nsymcrys
                  ja = ieqatom (ia, is, isym)
                  If ( .Not. done(ja)) Then
                     jas = idxas (ja, is)
! parallel transport of vector field (using operation inverse)
                     lspl = lsplsymc (isym)
                     ilspl = isymlat (lspl)
                     Do i = 1, ndmag
                        Call symrfmt (lrstp, is, symlatc(:, :, ilspl), &
                       & rvfmt(:, :, ias, i), rvfmt(:, :, jas, i))
                     End Do
! inverse of global proper rotation matrix in Cartesian coordinates
                     lspn = lspnsymc (isym)
                     ilspn = isymlat (lspn)
                     md = symlatd (ilspn)
                     sc (:, :) = dble (md) * symlatc (:, :, ilspn)
! global spin rotation of vector field
                     If (ncmag) Then
! non-collinear case
                        Do ir = 1, nrmt (is), lrstp
                           Do lm = 1, lmmaxvr
                              v (:) = rvfmt (lm, ir, jas, :)
                              rvfmt (lm, ir, jas, 1) = sc (1, 1) * v &
                             & (1) + sc (1, 2) * v (2) + sc (1, 3) * v &
                             & (3)
                              rvfmt (lm, ir, jas, 2) = sc (2, 1) * v &
                             & (1) + sc (2, 2) * v (2) + sc (2, 3) * v &
                             & (3)
                              rvfmt (lm, ir, jas, 3) = sc (3, 1) * v &
                             & (1) + sc (3, 2) * v (2) + sc (3, 3) * v &
                             & (3)
                           End Do
                        End Do
                     Else
! collinear case
                        Do ir = 1, nrmt (is), lrstp
                           Do lm = 1, lmmaxvr
                              rvfmt (lm, ir, jas, 1) = sc (3, 3) * &
                             & rvfmt (lm, ir, jas, 1)
                           End Do
                        End Do
                     End If
! mark atom as done
                     done (ja) = .True.
                  End If
               End Do
            End If
         End Do
      End Do
!---------------------------!
!     interstitial part     !
!---------------------------!
      Call symrvfir (ngvec, rvfir)
      Deallocate (rvfmt1, rvfmt2)
      Return
End Subroutine
!EOC
