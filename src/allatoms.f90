!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: allatoms
! !INTERFACE:
!
!
Subroutine allatoms
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Solves the Kohn-Sham-Dirac equations for each atom type in the solid and
!   finds the self-consistent radial wavefunctions, eigenvalues, charge
!   densities and potentials. The atomic densities can then be used to
!   initialise the crystal densities, and the atomic self-consistent potentials
!   can be appended to the muffin-tin potentials to solve for the core states.
!   Note that, irrespective of the value of {\tt xctype}, exchange-correlation
!   functional type 3 is used. See also {\tt atoms}, {\tt rhoinit},
!   {\tt gencore} and {\tt modxcifc}.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!   Modified for GGA, June 2007 (JKD)
!EOP
!BOC
      Implicit None
! always use LDA to setup atomic densities
      Integer, Parameter :: xctype_ = 3
      Integer, Parameter :: xcgrad_ = 0
      Integer :: is
! allocatable arrays
      Real (8), Allocatable :: rwf (:, :, :)
! allocate global species charge density and potential arrays
      If (allocated(sprho)) deallocate (sprho)
      Allocate (sprho(spnrmax, nspecies))
      If (allocated(spvr)) deallocate (spvr)
      Allocate (spvr(spnrmax, nspecies))
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rwf)
!$OMP DO
      Do is = 1, nspecies
         Allocate (rwf(spnrmax, 2, spnstmax))
         Call atom (input%groundstate%ptnucl, spzn(is), spnst(is), &
        & spn(:, is), spl(:, is), spk(:, is), spocc(:, is), xctype_, &
        & xcgrad_, input%groundstate%nprad, spnr(is), spr(:, is), &
        & speval(:, is), sprho(:, is), spvr(:, is), rwf)
         Deallocate (rwf)
      End Do
!$OMP END DO
!$OMP END PARALLEL
      Return
End Subroutine
!EOC
