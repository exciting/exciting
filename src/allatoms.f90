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
      Use FoX_wxml
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
      Integer :: is, i
      Logical :: dirac_eq
! allocatable arrays
      Real (8), Allocatable :: rwf (:, :, :)
      Type (xmlf_t), Save :: xf
      
      dirac_eq=(input%groundstate%CoreRelativity.eq."dirac")
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
        & xcgrad_, spnr(is), spr(:, is), &
        & speval(:, is), sprho(:, is), spvr(:, is), rwf,nrmt(is),dirac_eq)
         Deallocate (rwf)
      End Do
!$OMP END DO
!$OMP END PARALLEL

      Call xml_OpenFile ("atoms.xml", xf, replace=.True., pretty_print=.True.)
      Call xml_NewElement(xf,"atomlist")
      Call xml_NewElement (xf,"Hamiltonian")
      Call xml_AddAttribute (xf,"RelativityModel",trim(input%groundstate%CoreRelativity))
      Call xml_AddAttribute (xf,"xctype",xctype_)
      Call xml_EndElement (xf,"Hamiltonian")
      Do is = 1, nspecies
        Call xml_NewElement (xf,"atom")
        Call xml_AddAttribute (xf,"chemicalSymbol", trim(input%structure%speciesarray(is)%species%chemicalSymbol)) 
        Call xml_AddAttribute (xf,"species", trim(input%structure%speciesarray(is)%species%speciesfile))
        Call xml_NewElement (xf,"NumericalSetup")
        Call xml_AddAttribute (xf,"TotalNumberOfGridPoints",spnr(is))
        Call xml_AddAttribute (xf,"NumberOfMTGridPoints",nrmt(is))
        Call xml_AddAttribute (xf,"GridType",trim(input%groundstate%radialgridtype))
        Call xml_AddAttribute (xf,"rmin",spr(1, is))
        Call xml_AddAttribute (xf,"rmt",spr(nrmt(is), is))
        Call xml_AddAttribute (xf,"rmax",spr(spnr(is), is))
        Call xml_EndElement (xf,"NumericalSetup")
        Call xml_NewElement (xf,"spectrum")
        do i=1,spnst(is)
          Call xml_NewElement (xf,"state")
          Call xml_AddAttribute (xf,"n",spn(i, is))
          Call xml_AddAttribute (xf,"l",spl(i, is))
          Call xml_AddAttribute (xf,"kappa",spk(i, is))
          Call xml_AddAttribute (xf,"energy",speval(i, is))
          Call xml_EndElement (xf,"state")
        enddo
        Call xml_EndElement (xf,"spectrum")
        Call xml_EndElement (xf,"atom")
      End Do
      Call xml_EndElement (xf,"atomlist")
      Call xml_close (xf)

      Return
End Subroutine
!EOC
