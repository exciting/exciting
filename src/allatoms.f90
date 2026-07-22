!> Module for the calculation of the atomic densities and potentials
module allatoms
  use precision, only: dp, i32

  implicit none

  private 

  public :: calculate_allatoms
contains
!> Solves the Kohn-Sham-Dirac equations for each atom type in the solid and
!> finds the self-consistent radial wavefunctions, eigenvalues, charge
!> densities and potentials. The atomic densities can then be used to
!> initialise the crystal densities, and the atomic self-consistent potentials
!> can be appended to the muffin-tin potentials to solve for the core states.
!> Note that, irrespective of the value of `xctype`, exchange-correlation
!> functional type 3 is used. See also {\tt atoms}, {\tt rhoinit},
!> {\tt gencore} and {\tt modxcifc}.
subroutine calculate_allatoms(verbosity)
  use dfthalf, only: allocate_vhalf_global_arrays, dft_half_parameters, obtain_vhalf_potential
  use FoX_wxml, only: xmlf_t, xml_AddAttribute, xml_Close, xml_EndElement, &
                      xml_NewElement, xml_OpenFile
  use mod_atoms, only: natmtot, nspecies, speval, spk, spl, spn, spnrmax, &
    spnr, spnst, spnstmax, spocc, spr, sprho, spvr, spzn
  use mod_Gvector, only: ngrtot
  use mod_muffin_tin, only: lmmaxvr, nrmt, nrmtmax
  use modmpi, only : mpiglobal
  use modinput, only: input
  use to_char_conversion, only: to_char

  integer(i32), intent(in) :: verbosity
  ! always use LDA to setup atomic densities
  integer(i32), parameter :: xctype_ = 3
  integer(i32), parameter :: xcgrad_ = 0
  integer(i32) :: xctypearray(3)
  integer(i32) :: is, i
  logical :: dirac_eq, dft_half, point_nucleus
  character(100) :: buffer
  real(dp), allocatable :: rwf (:, :, :)
  Type (xmlf_t), Save :: xf
  type(dft_half_parameters) :: dft_half_params
  
  dirac_eq = (input%groundstate%CoreRelativity.eq."dirac")
  dft_half = associated( input%groundstate%dfthalf )
  point_nucleus = input%groundstate%ptnucl
  ! allocate global species charge density and potential arrays
  If (allocated(sprho)) deallocate (sprho)
  Allocate (sprho(spnrmax, nspecies))
  If (allocated(spvr)) deallocate (spvr)
  Allocate (spvr(spnrmax, nspecies))
  ! Allocate arrays for the DFT-1/2 part
  if ( dft_half ) call allocate_vhalf_global_arrays( ngrtot, &
    nrmtmax, lmmaxvr, natmtot )

  ! All libxc routines expect xctype(3) instead of a single integer
  xctypearray(1:3) = xctype_
  allocate (rwf(spnrmax, 2, spnstmax))
  Do is = 1, nspecies
    Call atom ( point_nucleus, spzn(is), spnst(is), &
      & spn(:, is), spl(:, is), spk(:, is), spocc(:, is), xctypearray, &
      & xcgrad_, spnr(is), spr(:, is), &
      & speval(:, is), sprho(:, is), spvr(:, is), rwf,nrmt(is),dirac_eq)
  end do

  if( dft_half ) then
    call dft_half_params%parse_input( input%groundstate%dfthalf, &
      input%structure%speciesarray, spnst )
    call dft_half_params%sanity_check( )
    call obtain_vhalf_potential( xctype_, xcgrad_, dirac_eq, point_nucleus, dft_half_params )
  end if

  if ((verbosity > 0) .and. mpiglobal%is_root ) then
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
        write(buffer,'(G22.12)') speval(i,is)
        Call xml_AddAttribute (xf,"energy",trim(adjustl(buffer)))
        Call xml_EndElement (xf,"state")
      enddo
      Call xml_EndElement (xf,"spectrum")
      Call xml_EndElement (xf,"atom")
    End Do
    Call xml_EndElement (xf,"atomlist")
    Call xml_close (xf)
  endif
End Subroutine
end module