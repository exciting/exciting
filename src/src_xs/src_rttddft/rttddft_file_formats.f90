module rttddft_file_formats
#include "asserts.fpp"
  implicit none

  private

  public :: restart_format, binary, hdf5

  !> Enumerator for the supported restart formats
  enum, bind(C)
    !> Variables of this enum type must be declared as `integer(kind(restart_format))`
    enumerator :: restart_format
    !> List of possible values for the enumerator
    enumerator :: binary, hdf5
  end enum

  !> Type that encapsulates metadata and configuration for HDF5 restart files
  type, public :: file_handler
    !> File format, according to the enum [[restart_format]]
    integer(kind(restart_format)) :: file_format
    !> Name of the HDF5 file (e.g., "rt.h5")
    character(len=:), allocatable :: file_name
    !> Path within the HDF5 file (typically the root group "/")
    character(len=:), allocatable :: path
  contains
    procedure :: assert_consistency => file_handler_assert_consistency
  end type

contains

  !> Check if the strings `file_name` and `path` and have been allocated
  subroutine file_handler_assert_consistency( this )
    class( file_handler ), intent(in) :: this

    CALL_ASSERT( allocated( this%file_name ), "string file_name must be allocated" )
    CALL_ASSERT( allocated( this%path ), "string path must be allocated" )
  end subroutine
end module