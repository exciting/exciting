module MD_io
  use asserts, only: assert
  use m_getunit, only: getunit
  use MD, only: force
  use mod_misc, only: filext
  use precision, only: dp, i32
  
  implicit none
  private

  type basic_io
    integer(i32)       :: file_unit
    character(len=50)  :: name

    contains

    procedure, private :: find_available_unit => basic_io_find_available_unit
    procedure, private :: set_name => basic_io_set_name
    procedure, private :: open_file => basic_io_open_file
    procedure, private :: close_file => basic_io_close_file
    procedure, private :: write_to_file => basic_io_write_to_file
  end type

  type, public :: MD_out
    logical                      :: all_force_contributions
    type(basic_io), allocatable  :: positions_velocities_forces(:)
    type(basic_io), allocatable  :: F_EXT(:)
    type(basic_io), allocatable  :: F_HF(:)
    type(basic_io), allocatable  :: F_core(:)
    type(basic_io), allocatable  :: F_val(:)

    contains 

    procedure, public :: write_to_files => write_to_MD_outs
    procedure, public :: open_files => open_MD_outs
    procedure, public :: close_files => close_MD_outs
  end type

  character(len=*), parameter :: basic_name_pos_vel_forces = 'ATOM_'
  character(len=*), parameter :: basic_name_F_core = 'FCR_'
  character(len=*), parameter :: basic_name_F_EXT = 'FEXT_'
  character(len=*), parameter :: basic_name_F_HF = 'FHF_'
  character(len=*), parameter :: basic_name_F_val = 'FVAL_'

  character(len=*), parameter :: format_default = 'F20.10'
  character(len=*), parameter :: format_time = 'F12.4'
  character(len=*), parameter :: format_position = format_default
  character(len=*), parameter :: format_velocity = format_default
  character(len=*), parameter :: format_force = format_default

  contains
  !> Find unused integer to be used as unit as store it in `file_unit`
  subroutine basic_io_find_available_unit( this )
    class(basic_io), intent(inout) :: this
    call getunit( this%file_unit )
  end subroutine

  !> Set name (will be used as the file name)
  subroutine basic_io_set_name( this, name )
    class(basic_io), intent(inout)  :: this
    !> Name to be given to the file
    character(len=*), intent(in)    :: name
    this%name = trim( name )
  end subroutine  

  !> Open file
  subroutine basic_io_open_file( this )
    class(basic_io), intent(in) :: this
    open( this%file_unit, file=trim( this%name ), status='replace', &
      position='append' )
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
    ! Local variables
    integer(i32) :: i
    character(len=50) :: string

    allocate( this%positions_velocities_forces(n) )

    do i = 1, n
      call this%positions_velocities_forces(i)%find_available_unit()
      write(string,'(I4.4)') i
      call this%positions_velocities_forces(i)%set_name( &
        basic_name_pos_vel_forces//trim(string)//trim(filext) )
      call this%positions_velocities_forces(i)%open_file()
    end do

    this%all_force_contributions = all_force_contributions
    
    if( all_force_contributions ) then
      allocate( this%F_core(n), this%F_EXT(n), this%F_HF(n), this%F_val(n) )
      do i = 1, n
        write(string,'(I4.4)') i
        call this%F_core(i)%find_available_unit()
        call this%F_core(i)%set_name( basic_name_F_core//trim(string)//trim(filext) )
        call this%F_core(i)%open_file()
        call this%F_EXT(i)%find_available_unit()
        call this%F_EXT(i)%set_name( basic_name_F_EXT//trim(string)//trim(filext) )
        call this%F_EXT(i)%open_file()
        call this%F_HF(i)%find_available_unit()
        call this%F_HF(i)%set_name( basic_name_F_HF//trim(string)//trim(filext) )
        call this%F_HF(i)%open_file()
        call this%F_val(i)%find_available_unit()
        call this%F_val(i)%set_name( basic_name_F_val//trim(string)//trim(filext) )
        call this%F_val(i)%open_file()
      end do
    end if
  end subroutine

  subroutine write_to_MD_outs( this, time, atom_positions, atom_velocities, &
      atom_forces, forces_contr )
    class(MD_out), intent(in) :: this
    real(dp), intent(in) :: time
    real(dp), intent(in) :: atom_positions(:, :)
    real(dp), intent(in) :: atom_velocities(:, :)
    real(dp), intent(in) :: atom_forces(:, :)
    type(force), intent(in), optional :: forces_contr
    ! Local variables
    integer(i32) :: i, n
    character(len=100) :: format_string
    character(len=200) :: string

    n = size( this%positions_velocities_forces, 1 )
    call assert( size(atom_positions, 1) == 3, 'atom_positions must have 3 coordinates' )
    call assert( size(atom_velocities, 1) == 3, 'atom_velocities must have 3 coordinates' )
    call assert( size(atom_forces, 1) == 3, 'atom_forces must have 3 coordinates' )
    call assert( size(atom_positions, 2) == n, 'atom_positions must have n elements' )
    call assert( size(atom_velocities, 2) == n, 'atom_velocities must have n elements' )
    call assert( size(atom_forces, 2) == n, 'atom_forces must have n elements' )
    call assert( present(forces_contr) .eqv. this%all_force_contributions, &
      'present(forces_contr) must be equivalent to all_force_contributions' )

    do i = 1, n
      format_string = '('//trim(format_time)// &
        ',3'//trim(format_position)// &
        ',3'//trim(format_velocity)// &
        ',3'//trim(format_force)//')'
      write( string, format_string ) time, atom_positions(:, i), &
        atom_velocities(:, i), atom_forces(:, i)
      call this%positions_velocities_forces(i)%write_to_file( string )
    end do

    if( this%all_force_contributions ) then
      do i = 1, n
        format_string = '('//trim(format_time)//',3'//trim(format_force)//')'
        write( string, format_string ) time, forces_contr%core(:, i)
        call this%F_core(i)%write_to_file( string )
        write( string, format_string ) time, forces_contr%EXT(:, i)
        call this%F_EXT(i)%write_to_file( string )
        write( string, format_string ) time, forces_contr%HF(:, i)
        call this%F_HF(i)%write_to_file( string )
        write( string, format_string ) time, forces_contr%val(:, i)
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

end module