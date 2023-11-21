!> Global definitions for the xhdf5 module
module xhdf5_globals

#ifdef _HDF5_  
  use hdf5, only: HID_T, HSIZE_T, HSSIZE_T
#endif

  use precision, only: sp, dp 


  implicit none


  private


  !> Root group of any HDF5 file
  character(*), parameter, public :: h5group_root = './'

  !> HDF5 integer kinds
#ifdef _HDF5_
  integer, parameter, public :: hdf5_id = HID_T
  integer, parameter, public :: hdf5_size = HSIZE_T
  integer, parameter, public :: hdf5_ssize = HSSIZE_T
#else 
  integer, parameter, public :: hdf5_id = sp
  integer, parameter, public :: hdf5_size = sp
  integer, parameter, public :: hdf5_ssize = sp
#endif

  !> Value for undefined HDF5 id.
  integer(hdf5_id), parameter, public  :: h5id_undefined = -1
  !> Value for undefined HDF5 size.
  integer(hdf5_size), parameter, public :: h5size_undefined = -1

  !> Blank, to fill non needed character length.
  character(1), parameter, public :: hdf5_blank = ' '
  !> String to represent a `.true.` value.
  character(*), parameter, public :: true_string = 'true'
  !> String to represent a `.false.` value.
  character(*), parameter, public :: false_string = 'false'

end module xhdf5_globals