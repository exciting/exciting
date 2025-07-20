!> Module designed for the task QPEigenvalues
module task_QPEigenvalues
  use asserts, only: assert
  use constants, only: zzero
  use exciting_mpi, only: mpiinfo
  use gw_info, only: write_to_gwinfo, write_to_gwinfo_boxmessage
  use mod_bands, only: evalfv, bandstructure_analysis
  use mod_eigenvalue_occupancy, only: efermi
  use mod_kpointset, only: k_set
  use mod_kqpts, only: kpoints_sets
  use mod_selfenergy, only: evalks, evalqp, selfex, selfec, freq_selfc, sigc, znorm, read_selfec_from_files, &
    read_selfex_from_files, generate_frequency_grid_for_correlation_self_energy, plot_selfc
  use mod_vxc, only: read_vxcnn, vxcnn
  use modgw, only: nvelgw
  use modinput, only: input, gw_type
  use modmpi, only: terminate_if_false, terminate, mpiglobal
  use precision, only: dp, i32
  use quasiparticle_energies, only: write_qp_energies_text_format
  use to_char_conversion, only: to_char

  implicit none 

  private

  character(len=*), parameter :: task_name = "QPEigenvalues"
  character(len=*), parameter :: file_name_qp_energies = 'EVALQP'
  character(len=*), parameter :: extension_binary_format ='.OUT'

  enum, bind(C)
    enumerator :: method_to_obtain_Fermi_level
    enumerator :: libbzint, DFT_VBM_CBm
  end enum

  !> Interface to the parameters defined in the input file
  type task_QPEigenvalues_parameters
    private
    type(kpoints_sets) :: k_points
    logical :: print_sigma_c
    integer(kind(method_to_obtain_Fermi_level)) :: method_Fermi_level
  contains
    procedure :: parse_input, sanity_checks => task_QPEigenvalues_sanity_checks
  end type

  public :: execute_task_QPEigenvalues

contains
!> Convert a string to the enum [[method_to_obtain_Fermi_level]]
pure function string_to_method_to_obtain_Fermi_level( string ) result(method)
  !> String to be converted
  character(len=*), intent(in) :: string
  integer(kind(method_to_obtain_Fermi_level)) :: method

  method = DFT_VBM_CBm
  select case( trim( string ) )
    case("from_DFT_VBM_CBM_indexes")
      method = DFT_VBM_CBm
    case("from_QP_Eigenvalues")
      method = libbzint
  end select
end function

!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_kpt )
  class(task_QPEigenvalues_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in):: gw_inp
  !> maximum number of k-points
  integer(i32), intent(in) :: n_kpt

  call this%sanity_checks( gw_inp )
  call this%k_points%parse_input( gw_inp%taskGroup%QPEigenvalues%kpointsarray, n_kpt )
  this%print_sigma_c = gw_inp%printSelfC
  this%method_Fermi_level = string_to_method_to_obtain_Fermi_level( gw_inp%taskGroup%QPEigenvalues%FermiLevel )
end subroutine


!> Perform sanity checks on the input parameters in the `gw` element
subroutine task_QPEigenvalues_sanity_checks( this, gw_inp )
  class(task_QPEigenvalues_parameters), intent(in) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in) :: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%QPEigenvalues), &
    'Element QPEigenvalues must be present when executing '//'"'//task_name//'"' )

end subroutine

!> Execute a task QPEigenvalues calculation
subroutine execute_task_QPEigenvalues( first_band, last_band, k_points_irreducible, file_format )
  !> Index of the first KS band for which vxc is evaluated
  integer(i32), intent(in) :: first_band
  !> Index of the last KS band for which vxc is evaluated
  integer(i32), intent(in) :: last_band
  !> Structure with the full set of irreducible k-points
  type(k_set), intent(in) :: k_points_irreducible
  !> Format of input/ouput files
  character(len=*), intent(in) :: file_format

  integer(i32) :: i_start, i_end, i_VBM, i_CBm, unit
  real(dp) :: E_Fermi_GW, E_gap, dos_fermi
  type(task_QPEigenvalues_parameters) :: input_parameters
  type(k_set) :: k_points_used
  logical :: myrank_writes_to_outputs

  myrank_writes_to_outputs = ( mpiglobal%is_root )
  call input_parameters%parse_input( input%gw, k_points_irreducible%nkpt )
  call input_parameters%k_points%obtain_list_of_indexes()
  if( myrank_writes_to_outputs ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )

  associate( list_kpt => input_parameters%k_points%list_of_indexes )
    i_start = 1
    i_end = size( list_kpt )
    ! Consistency check for the method to calculate the Fermi level
    if ( input_parameters%method_Fermi_level==libbzint .and. size( list_kpt ) < k_points_irreducible%nkpt ) &
      call terminate( "Fermi level using libbzint only possible when calculation includes all k-points" )
    ! Attention: we are using `k_points_used` only to pack the following data from `k_points_irreducible`
    ! nkpt, vkl, and wkpt according to `list_kpt`
    k_points_used%nkpt = i_end - i_start + 1
    k_points_used%vkl = k_points_irreducible%vkl(:, list_kpt(i_start:i_end))
    k_points_used%wkpt = k_points_irreducible%wkpt(list_kpt(i_start:i_end))

    call read_vxcnn( file_format )
    
    call read_selfec_from_files( list_kpt(i_start:i_end), file_format )
    call read_selfex_from_files( list_kpt(i_start:i_end), file_format )
    if ( trim(input%gw%selfenergy%method) == "ac" ) then
      ! Analytical continuation of the correlation self-energy from the complex to the real frequency axis
      call generate_frequency_grid_for_correlation_self_energy( input%gw )
      call calcselfc_ac()
    end if

    if( input_parameters%print_sigma_c .and. myrank_writes_to_outputs) call plot_selfc( freq_selfc%freqs, list_kpt(i_start:i_end), selfec, first_band )
  
    if ( allocated(sigc) ) deallocate(sigc)
    allocate( sigc(first_band:last_band, i_start:i_end), source=zzero )
    if ( allocated(znorm) ) deallocate(znorm)
    allocate( znorm(first_band:last_band, i_start:i_end), source=0._dp )
    if( allocated(evalqp) ) deallocate( evalqp )
    allocate( evalqp(first_band:last_band, i_start:i_end), source=0._dp )

    ! Trick: wrap into evalfv only the eigenvalues belonging to the selected k-points
    if( allocated(evalks) ) deallocate( evalks )
    allocate( evalks(first_band:last_band, i_start:i_end), source=evalfv(first_band:last_band, list_kpt(i_start:i_end)) )
    deallocate( evalfv )
    allocate( evalfv, source=evalks )

    call solve_QP_equation()
    if( myrank_writes_to_outputs ) then 
      call write_qp_energies_text_format( list_kpt, k_points_used%vkl, k_points_used%wkpt, &
        first_band, evalks, evalqp, real( vxcnn%diag_elements, dp ), selfex, sigc, znorm )
      if( input_parameters%method_Fermi_level==libbzint ) then
        call fermi_exciting( .false., nvelgw, last_band-first_band+1, k_points_irreducible%nkpt, &
          evalqp(first_band:last_band,:), k_points_irreducible%ntet, k_points_irreducible%tnodes, &
          k_points_irreducible%wtet, k_points_irreducible%tvol, E_Fermi_GW, E_gap, dos_fermi )
      else
        i_VBM = index_VBM( evalks, efermi ) + first_band - 1
        i_CBm = index_CBm( evalks, efermi ) + first_band - 1
        E_Fermi_GW = 0.5_dp * ( maxval( evalqp(i_VBM, :) ) + minval( evalqp(i_CBm, :) ) )
      end if
      call bandstructure_analysis( 'G0W0 band structure, assuming Efermi = ' // to_char(E_Fermi_GW), &
        first_band, evalqp, E_Fermi_GW, .false., list_kpt(i_start:i_end), k_points_used%vkl )
      call putevalqp( file_name_qp_energies//extension_binary_format, &
        k_points_used, first_band, last_band, evalks, efermi, evalqp, E_Fermi_GW )
    end if

  end associate

end subroutine

pure integer(i32) function index_VBM( eigs, E_Fermi )
  real(dp), intent(in) :: eigs(:, :)
  real(dp), intent(in) :: E_Fermi

  index_VBM = maxval( maxloc( eigs, dim = 1, mask = eigs < E_Fermi ) )
end function

pure integer(i32) function index_CBm( eigs, E_Fermi ) 
  real(dp), intent(in) :: eigs(:, :)
  real(dp), intent(in) :: E_Fermi

  index_CBm = minval( minloc( eigs, dim=1, mask = eigs > E_Fermi ) )
end function

end module