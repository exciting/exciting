module hdf5_dataset_utils

#ifdef _HDF5_  
  use hdf5
#endif

  use xhdf5_error_handling, only: handle_hdf5_error, xhdf5_assert
  use xhdf5_globals, only: hdf5_id, hdf5_size, h5id_undefined
  use mpi_utils, only: mpi_comm_type

  implicit none 

  private
  public :: hdf5_get_dataset_shape


  contains
  

  !> Create a new `group` at `h5path`, associated to a link in an HDF5 file or group associated with `h5id`.
  subroutine hdf5_get_dataset_shape(mpi_comm, h5id_file, h5path, dataset, dataset_shape)
    !> MPI communicator handle
    type(mpi_comm_type), intent(in) :: mpi_comm
    !> Identifier of the file or group, used by HDF5.
    integer(hdf5_id), intent(in) :: h5id_file
    !> Path to an object in an HDF5 file. Can be given as `path/to/object`.
    character(*), intent(in) :: h5path
    !> Name of the dataset to obtain shape of.
    character(*), intent(in) :: dataset
    !> Shape of the dataset.
    integer(hdf5_size), allocatable, intent(out) :: dataset_shape(:)

    integer :: rank
    integer(hdf5_size), allocatable :: dataset_shape_max(:)

#ifdef _HDF5_
    integer :: h5err
    integer(hdf5_id) :: h5id_group, h5id_dset, h5id_dspace, h5id_dtype

    call xhdf5_assert(mpi_comm, h5id_file /= h5id_undefined, &
            'Error(xhdf5%write): HDF5 file is not initialized.')

    call h5gopen_f(h5id_file, h5path, h5id_group, h5err)
    call handle_hdf5_error(mpi_comm, 'h5gopen_f' ,h5err)

    call h5dopen_f(h5id_group, dataset, h5id_dset, h5err)
    call handle_hdf5_error(mpi_comm, 'h5dopen_f', h5err)

    call h5dget_space_f(h5id_dset, h5id_dspace, h5err)
    call handle_hdf5_error(mpi_comm, 'h5dget_space_f', h5err)

    call h5sget_simple_extent_ndims_f(h5id_dspace, rank, h5err)
    call handle_hdf5_error(mpi_comm, 'h5dget_space_f', h5err)

    allocate(dataset_shape(rank))
    allocate(dataset_shape_max(rank))
    call H5Sget_simple_extent_dims_f(h5id_dspace, dataset_shape, dataset_shape, h5err)
    if (h5err == rank) h5err = 0
    call handle_hdf5_error(mpi_comm, 'H5Sget_simple_extent_dims', h5err)

    call h5dclose_f(h5id_dset, h5err)
    call handle_hdf5_error(mpi_comm, 'h5dclose_f', h5err)

    call h5gclose_f(h5id_group, h5err)
    call handle_hdf5_error(mpi_comm, 'h5gclose_f', h5err)
#endif    
  end subroutine hdf5_get_dataset_shape

end module hdf5_dataset_utils  
