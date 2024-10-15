! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! Created Jan 2021 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module implementing general initializations for RT-TDDFT
module rttddft_init
  use constants, only: zzero
  use m_gndstateq, only: gndstateq
  use MD, only: force, MD_input_keys
  use mod_APW_LO, only: apwordmax
  use mod_atoms, only: natmtot, nspecies, idxas, natoms
  use mod_bands, only: evalfv, nomax, numin, ikcbm, ikvbm, ikvcm
  use mod_core_states, only: ncg
  use mod_eigensystem, only: nmatmax
  use mod_eigenvalue_occupancy, only: occsv, nstfv, nstsv, efermi
  use mod_gkvector, only: ngk, ngkmax, vgkl, gkc, tpgkc, sfacgk
  use mod_kpoint, only: vkl, nkpt
  use mod_muffin_tin, only: lmmaxapw
  use mod_misc, only: filext
  use modgw, only: kset
  use modinput, only: input, getstructHybrid, emptynode
  use modmpi, only: rank, mpi_env_k, distribute_loop, terminate_if_false
  use modxs, only: isreadstate0
  use precision, only: dp, i32
  use rttddft_GlobalVariables, only: B_past, B_time, mathcalH, mathcalB
  use rttddft_HamiltonianOverlap, only: UpdateHam
  use rttddft_input, only: pmat_keys
  use rttddft_hybrids, only: hybrids_used, Set_Dimension_mixed_product_basis, set_barecoul_basis
  use rttddft_io, only: file_pmat_exists, read_pmat, write_pmat, &
                        file_pmat_mt_exists, read_pmat_mt, write_pmat_mt, write_file_info, &
                        write_file_info_fill_line_with_char, get_filename_pmat, get_filename_pmat_mt
  use rttddft_pmat, only: Obtain_Pmat_LAPWLOBasis
  use rttddft_VectorPotential, only: Vector_Potential

  implicit none
  
  private
  
  public :: initialize_rttddft

contains
!> This subroutine initializes many global variables in a RT-TDDFT calculation.
subroutine initialize_rttddft( input_pmat, predictorCorrector, vec_pot, molecular_dynamics, &
  evecfv_gnd, evecfv_time, evecfv_save, evecsv, &
  overlap, ham_time, ham_past, apwalm, pmat, pmatmt )
  !> Argument that encapsulates the input options of the element pmat
  type(pmat_keys), intent(in) :: input_pmat
  !> if `.True`, the predictor corrector loop is employed
  logical, intent(in) :: predictorCorrector
  !> type that encapsulates the vector potential
  type(Vector_Potential), intent(in) :: vec_pot
  !> variable that is an interface to the input keys defined in `input.xml` inside the `MD` block
  type(MD_input_keys), intent(in) :: molecular_dynamics
  !> Basis-expansion coefficients of the groundstate KS-WFs
  !> (nmatmax, nstfv, first_kpt : last_kpt)
  complex(dp), allocatable, intent(out) :: evecfv_gnd(:, :, :)
  !> Basis-expansion coefficients of the KS-WFs at time \(t\)
  !> (nmatmax, nstfv, first_kpt : last_kpt)
  complex(dp), allocatable, intent(out) :: evecfv_time(:, :, :)
  !> Basis-expansion coefficients of the KS-WFs at time \(t\) - auxiliary 
  !> variable used in the predictor-corrector loop
  !> (nmatmax, nstfv, first_kpt : last_kpt)
  complex(dp), allocatable, intent(out) :: evecfv_save(:, :, :)
  !> Basis-expansion coefficients of the KS-WFs: second-variational coefficients
  !> (nstfv, nstfv, first_kpt : last_kpt)
  complex(dp), allocatable, intent(out) :: evecsv(:, :, :)
  !> Overlap matrix (of basis functions)
  !> (nmatmax, nmatmax, first_kpt : last_kpt)
  complex(dp), allocatable, intent(out) :: overlap(:, :, :)
  !> Hamiltonian matrix at current time \(t\)
  !> (nmatmax, nmatmax, first_kpt : last_kpt)
  complex(dp), allocatable, intent(out) :: ham_time(:, :, :)
  !> Hamiltonian matrix at previous time \(t - \Delta t \)
  !> (nmatmax, nmatmax, first_kpt : last_kpt)
  complex(dp), allocatable, intent(out) :: ham_past(:, :, :)
  !> Matching coefficients of the (L)APWs
  !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
  complex(dp), allocatable, intent(out) :: apwalm(:, :, :, :, :)
  !> Momentum matrix elements (projected onto the (L)APW+LO basis elements)
  !> (nmatmax, nmatmax, 3, first_kpt : last_kpt)
  complex(dp), allocatable, intent(out)  :: pmat(:, :, :, :)
  !> Muffin-tin part of the Momentum matrix
  !> (nmatmax, nmatmax, 3, natmtot, first_kpt : last_kpt)
  complex(dp), allocatable, intent(out)  :: pmatmt(:, :, :, :, :)



  integer                     :: ik, first_kpt, last_kpt
  character(len=*), parameter :: new_line = achar(13)//achar(10)
  real(dp)                    :: voff(3)

  ! Backup groundstate variables
  call backup0
  call backup1

  !--------------------------------------------!
  !     map xs parameters associated to gs     !
  !--------------------------------------------!
  if (input%xs%rgkmax == 0.d0) input%xs%rgkmax = input%groundstate%rgkmax
  if (hybrids_used()) call adjustments_for_Hybrid_RTTDDFT()
  call mapxsparameters
  ! Initialize universal variables
  call init0
  call init1
  call init2

  if (hybrids_used()) call init_hybrids()

  call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)

  ! Interface with input variables
  voff(1:3) = input%xs%vkloff(1:3) 

  !> Print to RTTDDFT_INFO that we will start the single-shot GS calculation
  if (rank == 0) then
    call write_file_info_fill_line_with_char('=')
    call write_file_info('Non-self-consistent GS for TDDFT calculations - started'//new_line)
  end if

  ! Read from STATE.OUT exclusively
  isreadstate0 = .true.

  ! One-shot GS calculation
  ! Since an XS calculation with Hybrid functionals uses the GS parameters, a one shot GS calculation serves no purpose
  if (.not. hybrids_used()) call gndstateq(voff, '_RTTDDFT.OUT')

  allocate( evecfv_gnd(nmatmax, nstfv, first_kpt : last_kpt), source = zzero )
  allocate( evecfv_time(nmatmax, nstfv, first_kpt : last_kpt) )
  allocate( evecsv(nstsv, nstsv, first_kpt : last_kpt) )
  allocate( overlap(nmatmax, nmatmax, first_kpt : last_kpt), source = zzero )
  allocate( ham_time(nmatmax, nmatmax, first_kpt : last_kpt), source = zzero )
  allocate( ham_past(nmatmax, nmatmax, first_kpt : last_kpt), source = zzero )
  if ( predictorCorrector ) then
    allocate( evecfv_save(nmatmax, nstfv, first_kpt : last_kpt) )
  end if
  allocate( apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt) )
  allocate( pmat(nmatmax, nmatmax, 3, first_kpt : last_kpt) )
  if ( molecular_dynamics%valence_corrections .or. molecular_dynamics%basis_derivative ) &
  allocate( pmatmt(nmatmax, nmatmax, 3, natmtot, first_kpt : last_kpt) )



  call allocate_globals( first_kpt, last_kpt, molecular_dynamics%on, &
    allocate_mathcalH=molecular_dynamics%valence_corrections, &
    allocate_mathcalB=molecular_dynamics%valence_corrections .or. molecular_dynamics%basis_derivative,&
    allocate_B=molecular_dynamics%basis_derivative )
  
  if ( rank == 0 ) call write_to_info( molecular_dynamics%on, &
  predictorCorrector, evecfv_gnd, evecfv_time, evecfv_save, evecsv, &
  overlap, ham_time, ham_past, apwalm, pmat, pmatmt )

  call read_WF_potential_rttddft( first_kpt, evecfv_gnd, evecfv_time, evecsv )

  if ( hybrids_used() ) then
    if ( input%xs%realTimeTDDFT%calcNonlocalCurrentDensity ) then
      ! In the current implementation, the Coulomb potential used for the non local potential is calculated in plane wave basis
      ! For details, please refer to Eq. 61 in doi:10.1016/j.cpc.2012.09.018
      call terminate_if_false( input%groundstate%Hybrid%BasisBareCoulomb == "pw", "For RTTDDFT with hybrids only input%hybrid%barecoul%basis=pw is supported")
      call set_barecoul_basis()
    end if
  end if

  do ik = first_kpt, last_kpt
    ! Matching coefficients (apwalm)
    call match(ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, 1, ik), sfacgk(:, :, 1, ik), apwalm(:, :, :, :, ik))
  end do

  if( input_pmat%read_pmat_from_file ) then 
    call terminate_if_false( file_pmat_exists(), 'File:'//trim( get_filename_pmat() )//' not found')
    call read_pmat( first_kpt, pmat, mpi_env_k )
    if ( molecular_dynamics%on ) then 
      call terminate_if_false( file_pmat_mt_exists(), 'File:'//trim( get_filename_pmat_mt() )//' not found')
      call read_pmat_mt( first_kpt, pmatmt, mpi_env_k )
    end if
  else
    if( allocated(pmatmt) ) then 
      call Obtain_Pmat_LAPWLOBasis( first_kpt, input_pmat%force_pmat_hermitian, apwalm, pmat, pmatmt )
    else
      call Obtain_Pmat_LAPWLOBasis( first_kpt, input_pmat%force_pmat_hermitian, apwalm, pmat )
    end if
  end if
  if( input_pmat%write_pmat_to_file ) then
    call write_pmat( first_kpt, pmat, mpi_env_k )
    if ( molecular_dynamics%on ) call write_pmat_mt( first_kpt, pmatmt, mpi_env_k )
  end if

  ! Hamiltonian at time t=0
  call UpdateHam( first_kpt, vec_pot%a_tot, predcorr=.False., calculateOverlap=.True., forcePmatHermitian=input_pmat%force_pmat_hermitian, &
    overlap=overlap, ham_time=ham_time, ham_past=ham_past, apwalm=apwalm, pmat=pmat, pmatmt=pmatmt, &
    update_mathcalH=allocated(mathcalH), update_mathcalB=allocated(mathcalB), update_pmat=.False. )
  ham_past(:, :, :) = ham_time(:, :, :)

end subroutine

!> Allocate global arrays
subroutine allocate_globals(first_kpt, last_kpt, ionDynamics, allocate_mathcalH, &
                            allocate_mathcalB, allocate_B)
  !> index of the first `k-point` to be considered in the sum
  integer(i32), intent(in) :: first_kpt
  !> index of the last `k-point` considered
  integer(i32), intent(in) :: last_kpt
  !> if `.True`, we need to allocate arrays for Ehrenfest molecular dynamics
  logical, intent(in) :: ionDynamics
  !> if `.True`, we need to allocate the global array `mathcalH`
  logical, intent(in) :: allocate_mathcalH
  !> if `.True`, we need to allocate the global array `mathcalB`
  logical, intent(in) :: allocate_mathcalB
  !> if `.True`, we need to allocate the global arrays `B_time` and `B_past`
  logical, intent(in) :: allocate_B


  if (ionDynamics) then
    if (allocate_mathcalH) allocate (mathcalH(nmatmax, nmatmax, 3, natmtot, last_kpt))
    if (allocate_mathcalB) allocate (mathcalB(nmatmax, nmatmax, 3, natmtot, first_kpt:last_kpt))
    if (allocate_B) then
      allocate (B_time(nmatmax, nmatmax, first_kpt:last_kpt), source=zzero)
      allocate (B_past(nmatmax, nmatmax, first_kpt:last_kpt), source=zzero)
    end if
  end if

end subroutine

!> Output general information to `RTTDDFT_INFO.OUT`
subroutine write_to_info( ionDynamics, predictorCorrector, evecfv_gnd, &
  evecfv_time, evecfv_save, evecsv, overlap, ham_time, ham_past, &
  apwalm, pmat, pmatmt )
  !> Are we performing an MD calculation?
  logical, intent(in)         :: ionDynamics
  !> if `.True`, the predictor corrector loop is employed
  logical, intent(in) :: predictorCorrector
  !> Basis-expansion coefficients of the groundstate KS-WFs
  complex(dp), intent(in) :: evecfv_gnd(:, :, :)
  !> Basis-expansion coefficients of the KS-WFs at time \(t\)
  complex(dp), intent(in) :: evecfv_time(:, :, :)
  !> Basis-expansion coefficients of the KS-WFs at time \(t\) - auxiliary 
  !> variable used in the predictor-corrector loop
  complex(dp), allocatable, intent(in) :: evecfv_save(:, :, :)
  !> Basis-expansion coefficients of the KS-WFs: second-variational coefficients
  complex(dp), intent(in) :: evecsv(:, :, :)
    !> Overlap matrix (of basis functions)
  complex(dp), intent(in) :: overlap(:, :, :)
  !> Hamiltonian matrix at current time \(t\)
  complex(dp), intent(in) :: ham_time(:, :, :)
  !> Hamiltonian matrix at previous time \(t - \Delta t \)
  complex(dp), intent(in) :: ham_past(:, :, :)
  !> Matching coefficients of the (L)APWs
  complex(dp), intent(in) :: apwalm(:, :, :, :, :)
  !> Momentum matrix elements (projected onto the (L)APW+LO basis elements)
  !> (nmatmax, nmatmax, 3, first_kpt : last_kpt)
  complex(dp), intent(in) :: pmat(:, :, :, :)
  !> Muffin-tin part of the Momentum matrix
  !> (nmatmax, nmatmax, 3, natmtot, first_kpt : last_kpt)
  complex(dp), allocatable, intent(in) :: pmatmt(:, :, :, :, :)


  character(len=100)          :: string
  character(len=*), parameter :: formatMemory = '(A40,F12.1)'
  integer(i32), parameter     :: MB = 1048576

  call write_file_info('Non-self-consistent GS for TDDFT calculations - finished')
  call write_file_info_fill_line_with_char('=')
  call write_file_info('Allocated memory (MiB per MPI process)')
  write (string, formatMemory) 'Coefficients to match LAPW functions:', dble(sizeof(apwalm)/MB)
  call write_file_info(string)
  write (string, formatMemory) 'Wavefunctions:', &
      dble((sizeof(evecfv_gnd) + sizeof(evecfv_time) + sizeof(evecsv))/MB)
  call write_file_info(string)
  write (string, formatMemory) 'Hamiltonian and Overlap matrices:', &
      dble((sizeof(overlap) + sizeof(ham_time) + sizeof(ham_past))/MB)
  call write_file_info(string)
  if (predictorCorrector) then
    write (string, formatMemory) 'Extra storage (predictor-corrector):', &
      dble((sizeof(ham_time) + sizeof(evecfv_save))/MB)
    call write_file_info(string)
  end if
  write (string, formatMemory) 'Momentum matrix:', dble((sizeof(pmat))/MB)
  call write_file_info(string)
  if (ionDynamics) then
    call write_file_info(string)
    write (string, formatMemory) 'Molecular Dynamics - Muffin-tin aux. matrices:', &
      dble((sizeof(pmatmt) + sizeof(mathcalH) + sizeof(mathcalB) + sizeof(B_time) + sizeof(B_past))/MB)
    call write_file_info(string)
  end if
  call write_file_info_fill_line_with_char('=')
  ! General info to be printed to RTTDDFT_INFO
  call write_file_info('Important output files: AVEC.OUT, PVEC.OUT, JIND.OUT.')
  call write_file_info('JIND.OUT contains the x, y, and z components of the current density.')
  call write_file_info('PVEC.OUT contains the x, y, and z components of the polarization vector.')
  call write_file_info('AVEC.OUT contains in each line 6 elements:')
  call write_file_info(': the x components of the induced and the total vector potential.')
  call write_file_info(': the y components of the induced and the total vector potential.')
  call write_file_info(': the z components of the induced and the total vector potential.')
end subroutine

!> checks for consistency between gs hybrid calculation and rttddft and initializes pointer
subroutine adjustments_for_Hybrid_RTTDDFT()
  use rttddft_hybrids, only: hybrids_used
  logical :: is_compatible
  !> Tetrahedron method is used for Hybrid calculations
  call terminate_if_false(input%groundstate%stypenumber==-1, "stypenumber in the groundstate element must be set to libbzint to use To run RTTDDFT on top of hybrid functional calculations in the xs element")
  call terminate_if_false(hybrids_used(), "The non local current density can be computed only when hybrid functionals are used")
  is_compatible = .false.
  call is_gs_input_compatible_with_xs(input, is_compatible)
  if( is_compatible .eqv. .false.) call terminate_if_false( is_compatible, 'ERROR(rttddft_init): Parameters form GS not compatible with XS!')
  ! This pointer needs to get associated for when the calculation is started on top of a hybrid gs calculation
  if (.not. associated(input%groundstate%Hybrid)) &
      input%groundstate%Hybrid => getstructHybrid(emptynode)
end subroutine

!> For RT-TDDFT with hybrid funcionals, the same parameters for xs and the gs should be taken. This subroutine checks for that
subroutine is_gs_input_compatible_with_xs( inp, is_compatible)
  use modinput, only: input_type
  !> Information from the input file
  type(input_type), intent(in) :: inp
  !> True, if the input is compatible with rttddft and hybrid functionals (otherwise false)
  logical, intent(out) :: is_compatible
  !> Incompatibility message rttddft and hybrid functionals
  integer :: i

  is_compatible = .true.

  if (inp%xs%nosym .eqv. inp%groundstate%nosym) is_compatible = .false.
  do i = 1, 3
    if (inp%xs%ngridk(i) == inp%groundstate%ngridk(i)) is_compatible = .false.
    if (inp%xs%vkloff(i) == inp%groundstate%vkloff(i)) is_compatible = .false.
  end do
  if (inp%xs%reducek .eqv. inp%groundstate%reducek) is_compatible = .false.
  if (inp%xs%rgkmax == inp%groundstate%rgkmax) is_compatible = .false.
  if (inp%xs%swidth == inp%groundstate%swidth) is_compatible = .false.
  if (inp%xs%nempty == inp%groundstate%nempty) is_compatible = .false.

end subroutine

!> read WF and potential from potential gs run. For hybrid functionals, the parameters are read from the PBE run
subroutine read_WF_potential_rttddft( first_kpt, evecfv_gnd, evecfv_time, evecsv )
  !> First k-point treated by this (MPI)rank
  integer, intent(in)         :: first_kpt
  !> Basis-expansion coefficients of the groundstate KS-WFs
  !> (nmatmax, nstfv, first_kpt : last_kpt)
  complex(dp), intent(out) :: evecfv_gnd(:, :, first_kpt :)
  !> Basis-expansion coefficients of the KS-WFs at time \(t\)
  !> (nmatmax, nstfv, first_kpt : last_kpt)
  complex(dp), intent(out) :: evecfv_time(:, :, first_kpt :) 
  !> Basis-expansion coefficients of the KS-WFs: second-variational coefficients
  !> (nstfv, nstfv, first_kpt : last_kpt)
  complex(dp), intent(out) :: evecsv(:, :, first_kpt :)

  integer                     :: ik, last_kpt
  logical                     :: file_exists
  character(len=50)           :: string

  last_kpt = ubound( evecfv_gnd, 3 )
  if (hybrids_used()) then
    inquire (File='STATE_PBE.OUT', Exist=file_exists)
    call terminate_if_false(file_exists, 'ERROR(rttddft_init): Start from GS calculation is not possible, STATE_PBE.OUT is missing!')
    isreadstate0 = .false. ! We read not only from STATE.OUT
    string = filext
    filext = '_PBE.OUT'
  else
    string = filext
    filext = '_RTTDDFT.OUT'
  end if
  call readstate        ! read the density and potentials from file
  call gencore          ! generate the core wavefunctions and densities
  call genmeffig
  call linengy          ! find the new linearization energies
  call genapwfr         ! generate the APW radial functions
  call genlofr          ! generate the local-orbital radial functions
  call olprad           ! compute the overlap radial integrals
  if ( hybrids_used() ) then
    filext = string
    call energykncr()       ! core kinetic energy
    call init_product_basis()
    call readstate()
    call readfermi()
    call read_vxnl()        !The non local potential is read out from file
    call genmeffig

    !----------------------------------------
    ! Read KS eigenvalues from file EVALFV.OUT
    !----------------------------------------
    if (allocated(evalfv)) deallocate (evalfv)
    allocate (evalfv(nstfv, kset%nkpt))
    evalfv(:, :) = 0.d0
    do ik = 1, kset%nkpt
      call getevalfv(kset%vkl(:, ik), evalfv(:, ik))
    end do

    ! VB / CB state index
    call find_vbm_cbm(1, nstfv, kset%nkpt, evalfv, efermi, nomax, numin, ikvbm, ikcbm, ikvcm)

    ! The matrix sizes of mixed product basis quantities depend on if the core electrons are treated as valence
    ! This code block sets the dimension accordingly and is later read out in UpdateNonlocalCurrentDensity
    if ((input%gw%coreflag == 'all') .or. (input%gw%coreflag == 'xal')) then
      call Set_Dimension_mixed_product_basis(nomax + ncg)
    else
      call Set_Dimension_mixed_product_basis(nomax)
    end if

    ! Set BZ integration weights
    call kintw()

    deallocate (evalfv)

  end if

  ! Get the eigenvectors and occupations from file
  do ik = first_kpt, last_kpt
    ! Eigenvectors (first and second-variational components)
    call getevecfv(vkl(:, ik), vgkl(:, :, :, ik), evecfv_gnd(:, :, ik))
    evecfv_time(1:nmatmax, 1:nstfv, ik) = evecfv_gnd(1:nmatmax, 1:nstfv, ik)
    call getevecsv(vkl(:, ik), evecsv(:, :, ik))
    call getoccsv(vkl(:, ik), occsv(:, ik))
  end do

  filext = string

end subroutine

end module rttddft_init
