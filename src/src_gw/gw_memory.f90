!> Module for reporting memory usage in GW calculations
module gw_memory
#include "asserts.fpp"
  use precision, only: dp, i32, long_int
  use xyaml, only: yaml_type

  implicit none

  private

  type(yaml_type) :: yaml_file
  character(len=*), parameter :: file_name = "GW_MEMORY.yaml"
  character(len=*), parameter :: string_sgi = "sgi"
  character(len=*), parameter :: string_mpwipw = "mpwipw" 
  character(len=*), parameter :: string_vmat = "vmat" 
  character(len=*), parameter :: string_epsilon = "epsilon" 
  character(len=*), parameter :: string_fnm = "fnm" 
  character(len=*), parameter :: string_minmmat = "minmmat" 
  character(len=*), parameter :: string_minm = "minm" 
  character(len=*), parameter :: string_temp_calcminm2 = "temp_calcminm2"

  enum, bind(C)
    enumerator :: valid_field
    enumerator :: field_sgi
    enumerator :: field_mpwipw
    enumerator :: field_vmat
    enumerator :: field_epsilon
    enumerator :: field_fnm
    enumerator :: field_minmmat
    enumerator :: field_minm
    enumerator :: field_temp_calcminm2
  end enum

  type, public :: field_and_values
    private
    integer(kind(valid_field)) :: field
    integer(i32), allocatable :: dims(:)
    integer(i32) :: unit_multiplier
  contains
    procedure :: init => initialize_field_and_values
  end type

  ! Public procedures
  public :: close_file_memory_usage, &
            open_file_memory_usage, &
            write_memory_usage
  ! Enum
  public :: valid_field, field_sgi, field_mpwipw, field_vmat, field_epsilon, &
    field_fnm, field_minmmat, field_minm, field_temp_calcminm2

contains

pure subroutine initialize_field_and_values( this, field, dims, unit_multiplier )
  class(field_and_values), intent(inout) :: this
  integer(kind(valid_field)), intent(in) :: field
  integer(i32), contiguous, intent(in) :: dims(:)
  integer(i32), intent(in) :: unit_multiplier

  this%field = field
  this%dims = dims
  this%unit_multiplier = unit_multiplier
end subroutine

function field_to_string( field ) result(string)
  integer(kind(valid_field)), intent(in) :: field
  character(len=:), allocatable :: string

  select case( field )
    case(field_sgi)
      string = string_sgi
    case(field_mpwipw)
      string = string_mpwipw
    case(field_vmat)
      string = string_vmat
    case(field_epsilon)
      string = string_epsilon
    case(field_fnm)
      string = string_fnm
    case(field_minmmat)
      string = string_minmmat
    case(field_minm)
      string = string_minm
    case(field_temp_calcminm2)
      string = string_temp_calcminm2
    case default
      CALL_ASSERT(.false., "undefined field")
  end select
end function

subroutine open_file_memory_usage( )
  call yaml_file%open_file( file_name=file_name )
end subroutine

subroutine write_memory_usage( element_name, fields )
  character(len=*), intent(in) :: element_name
  type(field_and_values), intent(in) :: fields(:)

  character(len=*), parameter :: string_complex_dp = "complex(dp)"
  integer(i32), parameter :: complex_dp_bytes = sizeof( cmplx(0._dp, 0._dp, kind=dp) )

  integer(i32) :: i, n
  call yaml_file%open_element( element_name )
  ! Each field is a new element with subfields
  n = size( fields )
  do i = 1, n
    call yaml_file%open_element( field_to_string( fields(i)%field ) )
    call yaml_file%write_field( "type", string_complex_dp )
    call yaml_file%write_field( "dims", fields(i)%dims )
    call yaml_file%write_field( "memory", (product(int(fields(i)%dims, kind=long_int))*complex_dp_bytes/fields(i)%unit_multiplier) )
    call yaml_file%close_element( )
  end do
  call yaml_file%close_element( )
end subroutine

subroutine close_file_memory_usage( )
  call yaml_file%close_file()
end subroutine

end module