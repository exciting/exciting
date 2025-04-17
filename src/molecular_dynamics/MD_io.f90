module MD_io
  use asserts, only: assert
  use file_utils, only: add_default_extension, copy_text_file, read_last_and_penultimate_lines_from_file
  use math_utils, only: all_zero
  use MD, only: force, trajectory
  use mod_misc, only: filext
  use modmpi, only: terminate_if_false
  use precision, only: dp, i32, str_64, str_128, str_256
  use to_char_conversion, only: to_char
  
  implicit none
  private

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
  subroutine basic_io_open_file( this, name )
    class(basic_io), intent(inout) :: this
    !> Name to be given to the file
    character(len=*), intent(in)    :: name

    call this%set_name( name )
    open( newunit=this%file_unit, file=trim( this%name ), status='replace' )
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

  subroutine open_MD_outs( this, n, all_force_contributions )
    class(MD_out), intent(inout) :: this
    !> Number of files in each member class
    integer(i32), intent(in)     :: n
    !> Output the individual contributions to the total force
    logical, intent(in)          :: all_force_contributions

    integer(i32) :: i

    allocate( this%positions_velocities_forces(n) )

    do i = 1, n
      call this%positions_velocities_forces(i)%open_file( &
        add_integer_and_extension_to_filename(basic_name_pos_vel_forces, i ) )
    end do

    this%all_force_contributions = all_force_contributions
    
    if( all_force_contributions ) then
      allocate( this%F_core(n), this%F_EXT(n), this%F_HF(n), this%F_val(n) )
      do i = 1, n
        call this%F_core(i)%open_file( add_integer_and_extension_to_filename( basic_name_F_core, i ) )
        call this%F_EXT(i)%open_file( add_integer_and_extension_to_filename( basic_name_F_EXT, i ) )
        call this%F_HF(i)%open_file( add_integer_and_extension_to_filename( basic_name_F_HF, i ) )
        call this%F_val(i)%open_file( add_integer_and_extension_to_filename( basic_name_F_val, i ) )
      end do
    end if
  end subroutine

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
    call assert( size(nuclei_motion%velocities, 2) == n, 'velocities must have n elements' )
    call assert( size(forces%total, 1) == 3, 'atom_forces must have 3 coordinates' )
    call assert( size(forces%total, 2) == n, 'atom_forces must have n elements' )
    
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

  subroutine close_MD_outs( this )
    class(MD_out), intent(in) :: this
    ! Local variables
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