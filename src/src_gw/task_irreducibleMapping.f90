!> This module contains a task to use symmetry operations to map the inverse dielectric
!> matrix computed in the irreducible wedge into the FBZ. 
!> Further information regarding the use of symmetry operations can be found in mod_gw_symmetry.f90
module task_irreducibleMapping

  use constants,    only: zzero
  use precision,    only: i32, dp, max_length => str_32, str_128
  use exciting_mpi, only: mpiinfo
  use modgw,        only: freq, kset, kqset
  use mod_kqpts,    only: kpoints_sets
  use modinput,     only: input, gw_type
  use modmpi,       only: terminate_if_false, mpiglobal, distribute_loop
  use mod_misc_gw,  only: Gamma, gammapoint
  use to_char_conversion, only: to_char

  implicit none 

  private

  public :: execute_task_irreducibleMapping

  !> Task name
  character(len=*), parameter :: task_name = "irreducibleMapping"

  !> Interface to the parameters defined in the input file
  type task_irreducibleMapping_parameters
    private
    type(kpoints_sets) :: q_points
    integer(i32) :: n_omega
    character(len=max_length) :: output_format
    real(dp) :: eigenvalue_cutoff_Coulomb_matrix
  contains
    procedure :: parse_input, sanity_checks
  end type task_irreducibleMapping_parameters

contains

  !> Parses the input into the task_irreducibleMapping_parameters type
  subroutine parse_input( this, gw_inp, n_qpt )
    class(task_irreducibleMapping_parameters), intent(inout) :: this
    !> type with the variables given in the input file
    type(gw_type), intent(in) :: gw_inp
    !> maximum number of q-points (in the irreducible wedge)
    integer(i32), intent(in) :: n_qpt

    call this%sanity_checks( gw_inp )
    call this%q_points%parse_input( gw_inp%taskGroup%irreducibleMapping%qpointsarray, n_qpt )
    this%n_omega = gw_inp%freqgrid%nomeg
    this%output_format = trim( adjustl( gw_inp%taskGroup%outputFormat ) )
    this%eigenvalue_cutoff_Coulomb_matrix = gw_inp%barecoul%barcevtol
  end subroutine


  !> Perform sanity checks on the input parameters in the `gw` element
  subroutine sanity_checks( this, gw_inp )
    class(task_irreducibleMapping_parameters), intent(in) :: this
    !> type with the variables given in the input file
    type(gw_type), intent(in) :: gw_inp

    call terminate_if_false( associated(gw_inp%taskGroup%irreducibleMapping), &
      'Element irreducibleMapping must be present when executing '//'"'//task_name//'"' )
    call terminate_if_false( associated(gw_inp%freqgrid), &
      'Element freqgrid must be present when executing '//'"'//task_name//'"' )

  end subroutine sanity_checks

  !> Task driver for the symmetry regeneration of the inverse dielectric matrix
  subroutine execute_task_irreducibleMapping(n_qpoints_max, file_format)

    use mod_symmetry,    only: symlat, lsplsymc, nsymcrys, find_equivalent_wavevectors
    use mod_gw_symmetry, only: rotate_from_qa_to_qb_matrices
    use gw_info,         only: write_to_gwinfo, write_to_gwinfo_boxmessage
    use mod_dielectric_function, only: epsilon, read_inverse_epsilon_from_file, write_inverse_epsilon_to_file

    !> The number of points in the irreducible wedge
    integer(i32), intent(in) :: n_qpoints_max
    !> The file format
    character(len=*), intent(in) :: file_format

    integer(i32) :: iq, i, i_start, i_end, omega_i, omega_f
    integer(i32) :: iqfbz_representative
    type(task_irreducibleMapping_parameters) :: input_parameters
    type(mpiinfo) :: mpi_environment_qpoints

    integer(i32), allocatable :: iqeq_list(:), isymeq_list(:)
    integer(i32) :: n_iqeq, iqeq, sys_error
    complex(dp), allocatable  :: epsilon_iq(:,:,:)

    if( mpiglobal%rank == 0 ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )
    call input_parameters%parse_input( input%gw, n_qpoints_max )
    call input_parameters%q_points%obtain_list_of_indexes()
    call mpi_environment_qpoints%init( mpiglobal%comm )
    call distribute_loop( mpi_environment_qpoints, size(input_parameters%q_points%list_of_indexes), i_start, i_end)

    omega_i = 1
    omega_f = input_parameters%n_omega
 
    do i = i_start, i_end

      ! Get the iq point in the irreducible list
      iq = input_parameters%q_points%list_of_indexes(i)

      ! Get the representative of the equivalence class
      iqfbz_representative = kset%ikp2ik(iq)

      if( mpiglobal%rank == 0 ) call write_to_gwinfo( '('//task_name//'): q-point cycle, iq (IBZ), iq (FBZ) = ' &
                                                      // to_char(iq) // ' ' // to_char(iqfbz_representative) )

      ! Create a symbolic link for the representative point
      call execute_command_line('ln -sf INVERSE-EPSILON_IQ'//to_char(iq)// &
        '.OUT INVERSE-EPSILON_Q'//to_char(iqfbz_representative)//'.OUT > /dev/null 2>&1', exitstat=sys_error)

      call terminate_if_false(sys_error == 0, "Error(execute_task_irreducibleMapping): error creating symbolic links.")

      ! Check if the representative is Gamma
      Gamma = gammapoint(kqset%vqc(:,iqfbz_representative))
      
      ! Gamma point does not generate any other point, so it can be ignored 
      if (Gamma) cycle

      ! Get the inverse dielectric matrix
      call read_inverse_epsilon_from_file(iq, Gamma, file_format, .true.)
      call move_alloc(epsilon, epsilon_iq)

      ! Get the points in the equivalance class
      call find_equivalent_wavevectors( 3, kqset%vql(:, iqfbz_representative), kqset%vql(:,:), kqset%nkpt, &
                                        symlat(:, :, lsplsymc(1:nsymcrys)), nsymcrys, iqeq_list, isymeq_list, unique_only=.true.)
      
      ! Iterate over the equivalent points and regenerate using symmetry the inverse dielectric matrix
      do n_iqeq = 1, size(iqeq_list)
        
        ! Get equivalent point
        iqeq = iqeq_list(n_iqeq)

        ! Ignore identity
        if (iqeq == iqfbz_representative) cycle

        if( mpiglobal%rank == 0 ) call write_to_gwinfo( '     - Mapping inverse dielectric matrix : ' &
                                                        // to_char(iqfbz_representative) // ' => ' // to_char(iqeq) )

        ! Rotate the inverse dielectric matrix
        call rotate_from_qa_to_qb_matrices(iqfbz_representative, iqeq, input_parameters%eigenvalue_cutoff_Coulomb_matrix, &
                                           file_format, epsilon_iq, epsilon)

        ! Save the inverse of epsilon computed from the representative of the class using crystal symmetry
        call write_inverse_epsilon_to_file(iqeq, Gamma, file_format)

        ! Free epsilon
        deallocate(epsilon)

      end do

    end do 

  end subroutine execute_task_irreducibleMapping

end module task_irreducibleMapping
