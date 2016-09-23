!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: symrf
! !INTERFACE:
!
!
Subroutine symrf (lrstp, rfmt, rfir)
! !USES:
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   rfmt  : real muffin-tin function (inout,real(lmmaxvr,nrmtmax,natmtot))
!   rfir  : real intersitial function (inout,real(ngrtot))
! !DESCRIPTION:
!   Symmetrises a real scalar function defined over the entire unit cell using
!   the full set of crystal symmetries. In the muffin-tin of a particular atom
!   the spherical harmonic coefficients of every equivlent atom are rotated and
!   averaged. The interstitial part of the function is first Fourier transformed
!   to $G$-space, and then averaged over each symmetry by rotating the Fourier
!   coefficients and multiplying them by a phase factor corresponding to the
!   symmetry translation. See routines {\tt symrfmt} and {\tt symrfir}.
!
! !REVISION HISTORY:
!   Created May 2007 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lrstp
      Real (8), Intent (Inout) :: rfmt (lmmaxvr, nrmtmax, natmtot)
      Real (8), Intent (Inout) :: rfir (ngrtot)
! local variables
      Integer :: is, ia, ja, ias, jas, ir
      Integer :: isym, lspl, ilspl
      Real (8) :: t1
! automatic arrays
      Logical :: done (natmmax)
! allocatable arrays
      Real (8), Allocatable :: rfmt1 (:, :, :), rfmt2 (:, :)
      Allocate (rfmt1(lmmaxvr, nrmtmax, natmmax))
      Allocate (rfmt2(lmmaxvr, nrmtmax))
      t1 = 1.d0 / dble (nsymcrys)
!-------------------------!
!     muffin-tin part     !
!-------------------------!
      Do is = 1, nspecies
! make a copy of the input function
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nrmt (is), lrstp
               rfmt1 (:, ir, ia) = rfmt (:, ir, ias)
            End Do
         End Do
         done (:) = .False.
! loop over atoms
         Do ia = 1, natoms (is)
            If ( .Not. done(ia)) Then
               ias = idxas (ia, is)
               Do ir = 1, nrmt (is), lrstp
                  rfmt (:, ir, ias) = 0.d0
               End Do
! loop over crystal symmetries
               Do isym = 1, nsymcrys
! index to spatial rotation lattice symmetry
                  lspl = lsplsymc (isym)
! equivalent atom index (symmetry rotates atom ja into atom ia)
                  ja = ieqatom (ia, is, isym)
! apply the rotation to the muffin-tin function
                  Call symrfmt (lrstp, is, symlatc(:, :, lspl), &
                 & rfmt1(:, :, ja), rfmt2)
! accumulate in original function array
                  Do ir = 1, nrmt (is), lrstp
                     rfmt (:, ir, ias) = rfmt (:, ir, ias) + rfmt2 (:, &
                    & ir)
                  End Do
               End Do
! normalise
               Do ir = 1, nrmt (is), lrstp
                  rfmt (:, ir, ias) = t1 * rfmt (:, ir, ias)
               End Do
               done (ia) = .True.
! rotate into equivalent atoms
               Do isym = 1, nsymcrys
                  ja = ieqatom (ia, is, isym)
                  If ( .Not. done(ja)) Then
                     jas = idxas (ja, is)
                     lspl = lsplsymc (isym)
! inverse symmetry (which rotates atom ia into atom ja)
                     ilspl = isymlat (lspl)
! rotate symmetrised function into equivalent muffin-tin
                     Call symrfmt (lrstp, is, symlatc(:, :, ilspl), &
                    & rfmt(:, :, ias), rfmt(:, :, jas))
                     done (ja) = .True.
                  End If
               End Do
            End If
         End Do
      End Do
!---------------------------!
!     interstitial part     !
!---------------------------!
      Call symrfir (ngvec, rfir)
      Deallocate (rfmt1, rfmt2)
      Return
End Subroutine
!EOC
