!> Module for writing files in the `yaml` format. 
!> All writing has an instantanious effect on the file on the disc.
module xyaml
  use asserts, only: assert
  use precision, only: sp, dp, str_64, str_256
  use to_char_conversion, only: to_char
  use modmpi, only: terminate_if_false

  implicit none

  private
  public :: yaml_type

  character(2), parameter :: default_ident = "  "

  !> Type for managing writing data to disk in `yaml` format.
  type yaml_type
    !> Output unit of the file
    integer, private :: file_unit = -999
    !> Name of the file
    character(str_64), private :: file_name = ' '
    !> Identifier for the file being already initialized.
    !> Will set to `.true.` after intializing a new file is initialized and
    !> never mutated after wards.
    logical, private :: is_initialized = .false.
    !> Identifier for the file being open or closed
    logical, private :: is_open = .false.
    !> Element depth of the file
    integer, private :: depth = 0

  contains

    procedure :: open_file, close_file
    procedure :: open_element, close_element

    procedure :: write_field_scalar, write_field_vector, write_field_matrix
    generic :: write_field => write_field_scalar, write_field_vector, write_field_matrix
  end type yaml_type

  contains

  !> Open or create a `yaml` file for writing to it. 
  !>
  !> If the file already exists but is closed, it is opened again. Name and unit are not mutated! 
  !> If the file does not exist yet, it is created and opened. 
  !> If the file is already open, the routine has no effect at all.
  subroutine open_file(this, file_unit, file_name)
    !> `yaml` file container
    class(yaml_type), intent(inout) :: this
    !> Output unit of the file
    integer, intent(in), optional :: file_unit
    !> File name; it should end with `.yaml` or `.yml`
    character(*), intent(in), optional :: file_name


    integer :: ios
    character(str_256) :: err_msg

    if (this%is_open) then
      return
    
    else if (this%is_initialized) then
      open(unit=this%file_unit, file=this%file_name, iostat=ios, access='APPEND')

    else
      call assert(present(file_unit), 'file_unit is not given for uninitialized yaml file.')
      call assert(file_unit > 0, 'file_unit < 0')
      this%file_unit = file_unit

      call assert(present(file_name), 'file_name is not given for uninitialized yaml file.')
      this%file_name = trim(adjustl(file_name))

      open(unit=this%file_unit,  file=this%file_name, iostat=ios)
      write(err_msg, '("File ", A, " cannot be opened in unit ", I3, ".")') file_name, file_unit
      call terminate_if_false(ios == 0, err_msg)
      this%is_open = .true.
      this%is_initialized = .true.
    end if 
  end subroutine open_file

  !> Close an open `yaml` file. 
  !>
  !> If the file is already closed, the routine has no effect.
  subroutine close_file(this)
    !> `yaml` file container
    class(yaml_type), intent(inout) :: this

    if (this%is_open) then
      close(unit=this%file_unit)
      this%is_open = .false.
    end if
  end subroutine close_file

  !> Open a new element in a yaml file. 
  !>
  !> All fields or elements added after this command before closing it will be part of this element.
  subroutine open_element(this, element_name)
    !> `yaml` file container
    class(yaml_type), intent(inout) :: this
    !> Name of the element.
    character(*), intent(in) :: element_name

    call idention(this%file_unit, this%depth)
    write(this%file_unit, '(A, ": ")') element_name
    this%depth = this%depth + 1
  end subroutine open_element

  !> Close the element most recently opened. 
  !>
  !> After calling this routine, the element can not be opened again and nothing can be added to this element.
  subroutine close_element(this)
    !> `yaml` file container
    class(yaml_type), intent(inout) :: this

    this%depth = this%depth - 1
  end subroutine close_element

  !> Write a scalar to a yaml file field.
  !> 
  !> All native fortran types are supported.
  subroutine write_field_scalar(this, element_name, input)
    !> `yaml` file container
    class(yaml_type), intent(in) :: this
    !> Name of the element
    character(*), intent(in) :: element_name
    !> Scalar variable to write to file.
    class(*), intent(in) :: input

    character(:), allocatable :: input_string

    select type(input)
    type is(logical)
      input_string = to_char(input)
    type is(integer)
      input_string = to_char(input)
    type is(real(sp))
      input_string = to_char(input)
    type is(real(dp))
      input_string = to_char(input)
    type is(complex(sp))
      input_string = to_char(input)
    type is(complex(dp))
      input_string = to_char(input)
    type is(character(*))
      input_string = input
    class default
      call assert(.false., 'Input type is not supported.')
    end select

    call idention(this%file_unit, this%depth)
    write(this%file_unit, '(A,": ",A)') element_name, input_string
  end subroutine write_field_scalar

  !> Write a vector to a yaml file field.
  !> 
  !> All native fortran types are supported expect **character** vectors.
  subroutine write_field_vector(this, element_name, input)
    !> `yaml` file container
    class(yaml_type), intent(in) :: this
    !> Name of the element
    character(*), intent(in) :: element_name
    !> Vector to write to file
    class(*), intent(in) :: input(:)

    character(:), allocatable :: input_string

    select type(input)
    type is(logical)
      input_string = to_char(input)
    type is(integer)
      input_string = to_char(input)
    type is(real(sp))
      input_string = to_char(input)
    type is(real(dp))
      input_string = to_char(input)
    type is(complex(sp))
      input_string = to_char(input)
    type is(complex(dp))
      input_string = to_char(input)
    class default
      call assert(.false., 'Input type is not supported.')
    end select

    call idention(this%file_unit, this%depth)
    write(this%file_unit, '(A,": ",A)') element_name, input_string
  end subroutine write_field_vector

  !> Write a matrix to a yaml file field.
  !> 
  !> Only **real(dp)** and **complex(dp)** are supported (Issue #128).
  subroutine write_field_matrix(this, element_name, input)
    !> `yaml` file container
    class(yaml_type), intent(in) :: this
    !> Name of the element
    character(*), intent(in) :: element_name
    !> Matrix to write to file
    class(*), intent(in) :: input(:, :)

    character(:), allocatable :: input_string

    select type(input)
    type is(real(dp))
      input_string = to_char(input)
    type is(complex(dp))
      input_string = to_char(input)
    class default
      call assert(.false., 'Input type is not supported.')
    end select

    call idention(this%file_unit, this%depth)
    write(this%file_unit, '(A,": ",A)') element_name, input_string
  end subroutine write_field_matrix

  !> Ident the writing correctly for a given depth.
  subroutine idention(file_unit, depth)
    !> Output unit of the file
    integer, intent(in) :: file_unit
    !> Depth of the idention
    integer, intent(in) ::depth
    integer :: i
    do i=1, depth
      write(file_unit, '(A)', advance="no") default_ident
    end do
  end subroutine idention

end module xyaml