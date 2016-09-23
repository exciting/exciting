
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: reformatdynamicalmatrices
! !INTERFACE:
subroutine reformatdynamicalmatrices
! !USES:
use modinput
Use mod_atoms
use mod_qpoint
use mod_constants
! !DESCRIPTION:
!   Collecting pieces of the dynamical matrices and assembling them in a nicer
!   way, such that $3\times 3$ matrices are displayed for each combination of
!   atoms and each $\bf{q}$ point.
!
! !REVISION HISTORY:
!   Created February 2010 (Sagmeister)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: n
! allocatable arrays
      Complex (8), Allocatable :: dyn (:, :)
      Complex (8), Allocatable :: dynq (:, :, :)
      Complex (8), Allocatable :: dynmatq(:,:,:,:,:,:,:)
! initialise universal variables
      Call init0
      Call init2
      n = 3 * natmtot
      allocate(dyn(3,n))
      Allocate (dynq(n, n, nqpt))
      Allocate (dynmatq(3,3,natmmax,nspecies,natmmax,nspecies,nqpt))
! read in the dynamical matrices
      Call readdyn (.false.,dynq)
! reorder dynamical matrices
      call reorderdynmat(dynq,dynmatq)
! write out unsymmetrized dynamical matrices
      call writedynamicalmatrices('DYNMAT.OUT',dynmatq)
! apply the acoustic sum rule
      Call sumrule (dynq)
! reorder dynamical matrices
      call reorderdynmat(dynq,dynmatq)
! write out unsymmetrized dynamical matrices including sumrule correction
      call writedynamicalmatrices('DYNMAT_SUMRULE.OUT',dynmatq)
! read in the dynamical matrices including symmetrization
      Call readdyn (.true.,dynq)
! reorder dynamical matrices
      call reorderdynmat(dynq,dynmatq)
! write out unsymmetrized dynamical matrices
      call writedynamicalmatrices('DYNMAT_SYM.OUT',dynmatq)
! apply the acoustic sum rule
      Call sumrule (dynq)
! reorder dynamical matrices
      call reorderdynmat(dynq,dynmatq)
! write out unsymmetrized dynamical matrices including sumrule correction
      call writedynamicalmatrices('DYNMAT_SYM_SUMRULE.OUT',dynmatq)
      Deallocate (dyn, dynq, dynmatq)
      Write (*,*)
      Write (*, '("Info(reformatdynamicalmatrices): reformatted dynamical matrices")')
      Write (*, '(" including symmetrization and/or application of accoustic sumrule")')
      Write (*, '(" written to DYNMAT.OUT, DYNMAT_SYM.OUT, DYNMAT_SUMRULE.OUT and DYNMAT_SYM_SUMRULE.OUT")')
      Write (*,*)
End Subroutine
!EOC
