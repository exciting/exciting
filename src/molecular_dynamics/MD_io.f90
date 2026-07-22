module MD_io
#include "asserts.fpp"
  use exciting_mpi, only: mpiinfo
  use file_utils, only: add_default_extension, copy_text_file, read_last_and_penultimate_lines_from_file
  use math_utils, only: all_zero
  use MD, only: force, trajectory
  use mod_misc, only: filext
  use modmpi, only: terminate_if_false
  use os_utils, only: join_paths
  use precision, only: dp, i32, str_64, str_128, str_256
  use rttddft_file_formats, only: file_handler, hdf5, binary
  use rttddft_io_hdf5, only: read_array_hdf5, write_array_hdf5
  use to_char_conversion, only: to_char
  
  implicit none
  private
  public :: read_trajectory, write_trajectory

  type :: basic_io
  private
    integer(i32) :: file_unit
    character(len=str_64)  :: name
  contains
  private
    procedure :: set_name => basic_io_set_name
    procedure :: open_file => basic_io_open_file
    procedure :: close_file => basic_io_close_file
    procedure :: write_to_file => basic_io_write_to_file
  end type

  type, public :: MD_out
    private
    logical                      :: all_force_contributions
    type(basic_io), allocatable  :: positions_velocities_forces(:)
    type(basic_io), allocatable  :: F_EXT(:)
    type(basic_io), allocatable  :: F_HF(:)
    type(basic_io), allocatable  :: F_core(:)
    type(basic_io), allocatable  :: F_val(:)
  contains 
    private
    procedure, public :: write_to_files => write_to_MD_outs
    procedure, public :: open_files => open_MD_outs
    procedure, public :: close_files => close_MD_outs
    procedure, public :: copy_files => copy_MD_outs
    procedure, public :: read_time_and_forces_from_files => read_last_time_and_forces
  end type

  integer(i32), parameter :: n_cartesian = 3

  character(len=*), parameter :: basic_name_pos_vel_forces = 'ATOM_'
  character(len=*), parameter :: basic_name_F_core = 'FCR_'
  character(len=*), parameter :: basic_name_F_EXT = 'FEXT_'
  character(len=*), parameter :: basic_name_F_HF = 'FHF_'
  character(len=*), parameter :: basic_name_F_val = 'FVAL_'

  character(len=*), parameter :: filename_trajectory = 'TRAJECTORY'

  character(len=*), parameter :: format_default = 'F20.10'
  character(len=*), parameter :: format_time = 'F12.4'
  character(len=*), parameter :: format_position = format_default
  character(len=*), parameter :: format_velocity = format_default
  character(len=*), parameter :: format_force = format_default

  contains
  !> Set the file name as `name`
  subroutine basic_io_set_name( this, name )
    class(basic_io), intent(inout)  :: this
    !> Name to be given to the file
    character(len=*), intent(in)    :: name
    this%name = trim( name )
  end subroutine  

  !> Open file
  subroutine basic_io_open_file( this, name, is_new )
    class(basic_io), intent(inout) :: this
    !> Name to be given to the file
    character(len=*), intent(in)    :: name
    !> If `.true.`, a new file should be created, erasing old ones (if they exist).   
    !> If `.false.`, append to an existing file
    logical, intent(in) :: is_new

    character(len=*), parameter :: status_new = "replace"
    character(len=*), parameter :: status_old = "old"
    character(len=*), parameter :: position_new = "rewind"
    character(len=*), parameter :: position_old = "append"

    call this%set_name( name )
    if( is_new ) then
      open( newunit=this%file_unit, file=trim( name ), status=status_new, position=position_new )
    else
      open( newunit=this%file_unit, file=trim( name ), status=status_old, position=position_old )
    end if
  end subroutine

  !> Write string to file
  subroutine basic_io_write_to_file( this, string )
    class(basic_io), intent(in) :: this
    character(len=*), intent(in) :: string
    write( this%file_unit, '(A)' ) trim( string )
  end subroutine

  !> Close file
  subroutine basic_io_close_file( this )
    class(basic_io), intent(in) :: this
    close( this%file_unit )
  end subroutine

  !> Open standard MD output files
  subroutine open_MD_outs( this, new, n, all_force_contributions )
    class(MD_out), intent(inout) :: this
    !> Same meaning as in [[basic_io_open_file]]
    logical, intent(in) :: new 
    !> Number of files in each member class
    integer(i32), intent(in)     :: n
    !> Output the individual contributions to the total force
    logical, intent(in)          :: all_force_contributions

    integer(i32) :: i

    allocate( this%positions_velocities_forces(n) )

    do i = 1, n
      call this%positions_velocities_forces(i)%open_file( &
        add_integer_and_extension_to_filename(basic_name_pos_vel_forces, i ), new )
    end do

    this%all_force_contributions = all_force_contributions
    
    if( all_force_contributions ) then
      allocate( this%F_core(n), this%F_EXT(n), this%F_HF(n), this%F_val(n) )
      do i = 1, n
        call this%F_core(i)%open_file( add_integer_and_extension_to_filename( basic_name_F_core, i ), new )
        call this%F_EXT(i)%open_file( add_integer_and_extension_to_filename( basic_name_F_EXT, i ), new )
        call this%F_HF(i)%open_file( add_integer_and_extension_to_filename( basic_name_F_HF, i ), new )
        call this%F_val(i)%open_file( add_integer_and_extension_to_filename( basic_name_F_val, i ), new )
      end do
    end if
  end subroutine

  !> Write to standard MD output files
  subroutine write_to_MD_outs( this, time, nuclei_motion, forces )
    class(MD_out), intent(in) :: this
    !> Time \( t \)
    real(dp), intent(in) :: time
    !> This argument packs nuclei positions and velocities at time \( t \)
    class(trajectory), intent(inout) :: nuclei_motion
    !> Forces on each atom
    type(force), intent(in) :: forces
    
    integer(i32) :: i, n
    character(len=str_128) :: format_string
    character(len=str_256) :: string

    n = size( this%positions_velocities_forces, 1 )
    call nuclei_motion%assert_consistency( )
    CALL_ASSERT( size(nuclei_motion%velocities, 2) == n, 'velocities must have n elements' )
    CALL_ASSERT( size(forces%total, 1) == 3, 'atom_forces must have 3 coordinates' )
    CALL_ASSERT( size(forces%total, 2) == n, 'atom_forces must have n elements' )
    
    do i = 1, n
      format_string = '('//trim(format_time)// &
        ',3'//trim(format_position)// &
        ',3'//trim(format_velocity)// &
        ',3'//trim(format_force)//')'
      write( string, format_string ) time, nuclei_motion%positions(:, i), &
        nuclei_motion%velocities(:, i), forces%total(:, i)
      call this%positions_velocities_forces(i)%write_to_file( string )
    end do

    if( this%all_force_contributions ) then
      do i = 1, n
        format_string = '('//trim(format_time)//',3'//trim(format_force)//')'
        write( string, format_string ) time, forces%core(:, i)
        call this%F_core(i)%write_to_file( string )
        write( string, format_string ) time, forces%EXT(:, i)
        call this%F_EXT(i)%write_to_file( string )
        write( string, format_string ) time, forces%HF(:, i)
        call this%F_HF(i)%write_to_file( string )
        write( string, format_string ) time, forces%val(:, i)
        call this%F_val(i)%write_to_file( string )
      end do
    end if
  end subroutine

  !> Close standard MD output files
  subroutine close_MD_outs( this )
    class(MD_out), intent(in) :: this
    
    integer(i32) :: i, n
    
    n = size( this%positions_velocities_forces, 1 )
    do i = 1, n
      call this%positions_velocities_forces(i)%close_file()
    end do
    if( this%all_force_contributions ) then
      do i = 1, n
        call this%F_core(i)%close_file()
        call this%F_EXT(i)%close_file()
        call this%F_HF(i)%close_file()
        call this%F_val(i)%close_file()
      end do
    end if
  end subroutine

  !> Read the last line of each file that stores positions, velocities, and forces
  subroutine read_last_time_and_forces( this, time, forces )
    class(MD_out), intent(in) :: this
    !> Time \(t\) to be read from the last line
    real(dp), intent(out) :: time
    !> Force on each atom, to be read from last and penultimate lines
    class(force), intent(inout) :: forces

    integer(i32) :: i, n
    character(len=str_256) :: last_line, penultimate_line
    real(dp) :: positions_ignore(n_cartesian), velocities_ignore(n_cartesian), t_aux
    real(dp), allocatable :: t(:), dt(:)

    n = find_number_of_files( basic_name_pos_vel_forces )
    if( allocated( forces%total ) ) deallocate( forces%total )
    if( allocated( forces%total_save ) ) deallocate( forces%total_save )
    allocate( forces%total(n_cartesian, n), forces%total_save(n_cartesian, n), t(n), dt(n) )
    do i = 1, n
      call read_last_and_penultimate_lines_from_file( add_integer_and_extension_to_filename( basic_name_pos_vel_forces, i ), &
        last_line, penultimate_line )
      read( last_line, * ) t(i), positions_ignore, velocities_ignore, forces%total(:, i)
      read( penultimate_line, * ) t_aux, positions_ignore, velocities_ignore, forces%total_save(:, i)
      dt(i) = t(i) - t_aux
    end do
    call terminate_if_false( all_zero( t-t(1) ), "Last time is not the same across the multiple output files")
    call terminate_if_false( all_zero( dt-dt(1) ), "Time step is not the same across the multiple output files")
    time = t(1)
  end subroutine

  !> Copy standard MD output files that have the additional extension `source_additional_extension`
  !> to files with the same base name but without the extension.
  subroutine copy_MD_outs( this, source_additional_extension, all_force_contributions )
    class(MD_out), intent(in) :: this
    !> Additional extension of source files (used to differentiate them from destination files)
    character(len=*), intent(in) :: source_additional_extension
    !> Output the individual contributions to the total force
    logical, intent(in)          :: all_force_contributions

    integer(i32) :: i, n

    n = find_number_of_files( basic_name_pos_vel_forces, source_additional_extension )
    do i = 1, n
      call wrapper_copy_txt_file_generic( basic_name_pos_vel_forces, i, source_additional_extension )
    end do

    if( all_force_contributions ) then
      do i = 1, n
        call wrapper_copy_txt_file_generic( basic_name_F_core, i, source_additional_extension )
        call wrapper_copy_txt_file_generic( basic_name_F_EXT, i, source_additional_extension )
        call wrapper_copy_txt_file_generic( basic_name_F_HF, i, source_additional_extension )
        call wrapper_copy_txt_file_generic( basic_name_F_val, i, source_additional_extension )
      end do
    end if
  contains
    subroutine wrapper_copy_txt_file_generic( file_name, m, src_extra_extension )
      character(len=*), intent(in) :: file_name
      integer(i32), intent(in) :: m
      character(len=*), intent(in) :: src_extra_extension

      call copy_text_file( source_name=add_integer_and_extension_to_filename( file_name, m ) // trim( src_extra_extension ), &
        destination_name=add_integer_and_extension_to_filename( file_name, m ) )
    end subroutine
  end subroutine

  !> Write positions and velocities to an output file
  subroutine write_trajectory( nuclei_motion, handler, mpi_env )
    !> This argument packs nuclei positions and velocities
    class(trajectory), intent(in) :: nuclei_motion
    !> File handler
    type(file_handler), intent(in) :: handler
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    integer(i32) :: unit

    call nuclei_motion%assert_consistency()

    select case( handler%file_format )
      case( binary )
        if( mpi_env%is_root ) then
          open( newunit=unit, file=add_default_extension( filename_trajectory ), action="write", form="unformatted", access="stream" )
          write( unit ) nuclei_motion%positions, nuclei_motion%velocities
          close( unit )
        end if
      case( hdf5 )
        call handler%assert_consistency( )
        call write_array_hdf5( handler%file_name, handler%path, filename_trajectory // "-" // add_default_extension("nuclei_positions"), &
          nuclei_motion%positions, mpi_env, serial_access=.true. )
        call write_array_hdf5( handler%file_name, handler%path, filename_trajectory // "-" // add_default_extension("nuclei_velocities"), &
          nuclei_motion%velocities, mpi_env, serial_access=.true. )
      case default
        CALL_ASSERT( .false., "Unrecognized format" )
    end select
  end subroutine

  !> Read positions and velocities to an output file
  subroutine read_trajectory( nuclei_motion, handler, mpi_env )
    !> This argument packs nuclei positions and velocities
    class(trajectory), intent(inout) :: nuclei_motion
    !> File handler
    type(file_handler), intent(in) :: handler
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    integer(i32) :: unit

    call nuclei_motion%assert_consistency()

    select case( handler%file_format )
      case( binary )
        open( newunit=unit, file=add_default_extension( filename_trajectory ), action="read", form="unformatted", access="stream" )
        read( unit ) nuclei_motion%positions, nuclei_motion%velocities
        close( unit )
      case( hdf5 )
        call handler%assert_consistency( )
        call read_array_hdf5( handler%file_name, handler%path, filename_trajectory // "-" // add_default_extension("nuclei_positions"), &
          nuclei_motion%positions, mpi_env )
        call read_array_hdf5( handler%file_name, handler%path, filename_trajectory // "-" // add_default_extension("nuclei_velocities"), &
          nuclei_motion%velocities, mpi_env )
      case default
    end select
  end subroutine

  !> Find the number of files in the current directory with name matching the pattern: 
  !> `add_integer_and_extension_to_filename( base_name, n )`.
  !> If `additional_extension` is present, then search for 
  !> `add_integer_and_extension_to_filename( base_name, n ) // trim( additional_extension )`
  integer(i32) function find_number_of_files( base_name, additional_extension ) result(n)
    !> Base name of the files to find
    character(len=*), intent(in) :: base_name
    !> Additional extension of the files to find
    character(len=*), optional, intent(in) :: additional_extension

    logical :: file_exists
    character(len=*), parameter :: no_ending = ""
    character(len=:), allocatable :: ending
    
    if( present( additional_extension ) ) then
      ending = trim( additional_extension )
    else
      ending = no_ending
    end if

    file_exists = .true.
    n = 0
    do while ( file_exists )
      n = n + 1
      inquire( file=add_integer_and_extension_to_filename( base_name, n ) // ending, exist=file_exists)
    end do
    n = n - 1
  end function

  !> Return a string made up from `file_name` and the integer `i`
  pure function combine_string_and_int( file_name, i ) result( name )
    !> base file name
    character(len=*), intent(in)  :: file_name
    !> integer to be combined
    integer(i32), intent(in) :: i
    !> file name combined with the integer `i`
    character(len=:), allocatable :: name

    integer(i32), parameter :: length = 4
    character(len=length) :: string

    write(string,'(I' // to_char( length ) // '.' // to_char( length )// ')') i
    name = trim( file_name )//trim( string )
  end function

  !> Return a string made up from `file_name`, the integer `i` and the default extension
  pure function add_integer_and_extension_to_filename( file_name, i ) result( name )
    !> base file name
    character(len=*), intent(in)  :: file_name
    !> integer to be combined
    integer(i32), intent(in) :: i
    !> file name with the integer `i` and the default extension
    character(len=:), allocatable :: name
    name = add_default_extension( combine_string_and_int( file_name, i ) )
  end function

end module