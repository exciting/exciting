!> This module contains the quatities, and the information
!> regarding self consistent GW
module mod_selfconsistent_gw

    use precision,  only: i32, dp, str_1024
    use modmpi,     only: terminate_if_false, mpiglobal, barrier, terminate
    use mod_mpi_gw, only: mpi_sum_array
    use gw_io,      only: read_from_file, write_to_file, open_file
    use mod_qsgw,   only: prepare_next_iteration_qsgw, prepare_current_iteration_qsgw, check_convergence_qsgw

    implicit none
    
    private

    public :: no_selfconsistent_gw, qsgw, &
              initialize_selfconsitent_gw, prepare_next_iteration, prepare_current_iteration, &
              is_gw_selfconsistent_flavour, gw_first_iteration, check_convergence

    !> Enumeration for different GW schemes
    enum, bind(c)
        enumerator :: no_selfconsistent_gw = 0
        enumerator :: qsgw = 1
        !enumerator :: evgw0 = 2
    end enum

    !> Selected GW scheme (default = no_selfconsistent_gw)
    integer, public, protected :: selfconsistent_type = no_selfconsistent_gw
    !> Iteration counter
    integer, public, protected, target :: iteration = 0
    !> Condition for convergence
    real(dp), public, protected, target :: selfconsistent_gw_eps = 1.0e-4_dp

    !> The file for the iteration
    character(len=*), public, parameter :: file_name_iteration = 'GW_ITERATION.OUT'

contains

    !> Checks the existance of a folder by creating a temporary
    !> file. This is done in this way because inquiring folders 
    !> is not consistent along compilers.
    logical function check_directory_existence(dir_path) result(exists)
        character(len=*), intent(in) :: dir_path

        character(len=str_1024) :: test_file
        integer :: unit_num, ios

        test_file = trim(dir_path) // '/.dirtest.tmp'

        open(newunit=unit_num, file=trim(test_file), status='new', action='write', iostat=ios)

        if (ios == 0) then
            exists = .true.
            close(unit_num, status='delete', iostat=ios)
        else
            exists = .false.
        end if

    end function check_directory_existence

    !> Initializes the module information
    !> and controls the self-consitent GW 
    !> files. 
    !> It call specialized functions depending on
    !> the SCF flavour.
    subroutine initialize_selfconsitent_gw()
        
        use modinput, only: input, gw_type

        !> type with the variables given in the input file
        type(gw_type), pointer :: gw_inp

        gw_inp => input%gw

        ! If GW block is not associated ignore
        if (.not. associated(gw_inp)) return

        ! If not defined ignore
        if (.not. associated(gw_inp%selfconsistency)) return

        ! If associated but using the old GW implementation
        ! raise a fatal error
        call terminate_if_false(gw_inp%taskname == 'taskGroup' .or. gw_inp%taskname == 'skip', &
                                "Error(initialize_selfconsitent_gw): Self-consistent " // &
                                "GW can only be used via 'taskGroup'. 'skip' is also accepted.")

        ! Select the kind of self-consistency
        select case(gw_inp%selfconsistency%type)
        case("oneshot")
            selfconsistent_type = no_selfconsistent_gw
        case("QSGW")
            selfconsistent_type = qsgw
        end select

        selfconsistent_gw_eps = gw_inp%selfconsistency%eps

    end subroutine initialize_selfconsitent_gw

    !> Creates a folder with files for the next
    !> GW iteration.
    subroutine prepare_next_iteration(file_format)

        use modinput, only: input
        !> The format of the large files to IO
        character(len=*), intent(in) :: file_format

        integer(i32) :: unit, ierror, warning_unit
        logical :: directory_already_exists

        ! In case one arrives to this function from G0W0 run
        ! simply ignore it
        if (selfconsistent_type == no_selfconsistent_gw) return

        ! Generate a folder to save the data for the next iteration
        ! and generate a file inidication the iteration that will run
        if (mpiglobal%rank == 0) then

            directory_already_exists = check_directory_existence("next_iteration")

            if (directory_already_exists) call terminate("Error(QSGW): next_iteration folder already exists. We stop to prevent overwritting.")
            
            call execute_command_line("mkdir -p next_iteration", cmdstat=ierror)
            call terminate_if_false(ierror == 0, "Error when creating the next_iteration folder.")
            
            call open_file("next_iteration/" // file_name_iteration, 'write', 'text', unit)
            write(unit, *) iteration + 1
            close(unit)
            
        end if

        call barrier(mpiglobal)

        ! Specific saves by code
        select case(selfconsistent_type)
        case(qsgw)
            call prepare_next_iteration_qsgw(iteration, file_format)
        case default
            call terminate("Error(save_iteration_information): shouldn't be here")
        end select

    end subroutine prepare_next_iteration

    !> Modifies globals to account for the previous self-consistent iteration
    !> In this way the self consistent calculations are transparent to the users
    subroutine prepare_current_iteration()

        integer(i32) :: unit, ierror
        logical :: not_first_iteration

        ! For non selfconsisten gw cases, or
        ! any other run ignore
        if (selfconsistent_type == no_selfconsistent_gw) return

        ! First check if the iteration file exist
        inquire(file=file_name_iteration, exist=not_first_iteration)

        if (not_first_iteration) then
            call open_file(file_name_iteration, 'read', 'text', unit)
            read(unit, *) iteration
            close(unit)
        end if

        select case(selfconsistent_type)
        case(qsgw)
            call prepare_current_iteration_qsgw(iteration)
        case default
            call terminate("Error(prepare_current_iteration): shouldn't be here")
        end select 

    end subroutine prepare_current_iteration

    !> Checks the convergence of the self-consistent GW scheme
    subroutine check_convergence(iteration, convergence_tolerance)

        integer(i32), intent(in) :: iteration
        real(dp),     intent(in) :: convergence_tolerance

        select case(selfconsistent_type)
        case(qsgw)
            call check_convergence_qsgw(iteration, convergence_tolerance)
        case(no_selfconsistent_gw)
            continue
        case default
            call terminate("Error(check_convergence): shouldn't be here")
        end select 

    end subroutine check_convergence

    !> 

    !> Checks if the input flavour is equal to the flavour
    !> of the GW self consistent.
    pure logical function is_gw_selfconsistent_flavour(flavour)
        integer, intent(in) :: flavour
        is_gw_selfconsistent_flavour = ( flavour == selfconsistent_type)
    end function is_gw_selfconsistent_flavour

    !> Checks if it is first iteration
    pure logical function gw_first_iteration() 
        gw_first_iteration = ( iteration == 0)
    end function gw_first_iteration

end module mod_selfconsistent_gw
