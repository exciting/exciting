
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: allatoms
! !INTERFACE:
subroutine allatoms
! !USES:
use modmain
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
implicit none
! always use LDA to setup atomic densities
integer, parameter :: xctype_=3
integer, parameter :: xcgrad_=0
integer is
! allocatable arrays
real(8), allocatable :: rwf(:,:,:)
allocate(rwf(spnrmax,2,spnstmax))
! allocate global species charge density and potential arrays
if (allocated(sprho)) deallocate(sprho)
allocate(sprho(spnrmax,nspecies))
if (allocated(spvr)) deallocate(spvr)
allocate(spvr(spnrmax,nspecies))
do is=1,nspecies
  call atom(spzn(is),spnst(is),spn(1,is),spl(1,is),spk(1,is),spocc(1,is), &
   xctype_,xcgrad_,nprad,spnr(is),spr(1,is),speval(1,is),sprho(1,is), &
   spvr(1,is),rwf)
end do
deallocate(rwf)
return
end subroutine
!EOC

