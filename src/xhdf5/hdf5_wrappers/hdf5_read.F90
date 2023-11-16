!> Wrapper for reading any data set from an HDF5 file
module hdf5_read
  use iso_c_binding, only: c_ptr

#ifdef _HDF5_  
  use hdf5
#endif

  use xhdf5_error_handling, only: handle_hdf5_error, xhdf5_assert
  use xhdf5_globals, only: hdf5_id, hdf5_size, hdf5_ssize, h5id_undefined
  use mpi_utils, only: mpi_comm_type
  use datainfo_utils, only: datainfo_type


  implicit none
  

  private
  public :: hdf5_read_dataset


  contains 
  

  !> Read a data set or chunk of an data set of any type and shape from file. 
  subroutine hdf5_read_dataset(mpi_comm, h5id_file, h5path, dataset, dataset_type, dataset_ptr, datainfo, serial_access)
    !> MPI communicator handle
    type(mpi_comm_type), intent(in) :: mpi_comm
    !> Identifier of the file or group, used by HDF5.
    integer(hdf5_id), intent(in) :: h5id_file
    !> Path to the group to create the data set.
    character(*), intent(in) :: h5path
    !> Name of the data set to read.
    character(*), intent(in) :: dataset
    !> Type id of the data set.
    integer(hdf5_id), intent(in) :: dataset_type
    !> Pointer to the first element of the data set to be written.
    type(c_ptr), intent(inout) :: dataset_ptr
    !> Type that holds information about the global and local data set shape and offset.
    type(datainfo_type) :: datainfo
    !> Set to `.true.` if only serial access is possible.
    logical, intent(in) :: serial_access

#ifdef _HDF5_
    integer(hdf5_id) :: h5id_group, h5id_fspace, h5id_dset, h5id_dspace, h5id_plist
    integer :: h5err

    call xhdf5_assert(mpi_comm, h5id_file /= h5id_undefined, &
            'Error(xhdf5%write): HDF5 file is not initialized.')

    if(serial_access) then
      call xhdf5_assert(mpi_comm, all(datainfo%offset == 0), 'serial_access is .true. but all(datainfo%offset == 0) is .false.')
    end if 

    ! Open group to read the data set from
    call h5gopen_f(h5id_file, h5path, h5id_group, h5err)
    call handle_hdf5_error(mpi_comm, 'h5gopen_f', h5err)

    call h5dopen_f(h5id_group, dataset, h5id_dset, h5err)
    call handle_hdf5_error(mpi_comm, 'h5dopen_f', h5err)

    ! Create a data set for the local chunk of the array and select the hyperslap to read from
    call h5screate_simple_f(datainfo%rank, datainfo%local_shape, h5id_dspace, h5err)
    call handle_hdf5_error(mpi_comm, 'h5screate_simple_f', h5err)

    call h5dget_space_f(h5id_dset, h5id_fspace, h5err)
    call handle_hdf5_error(mpi_comm, 'h5dget_space_f', h5err)

    call h5sselect_hyperslab_f(h5id_fspace, H5S_SELECT_SET_F, datainfo%offset, datainfo%local_shape, h5err)
    call handle_hdf5_error(mpi_comm, 'h5sselect_hyperslab_f', h5err)

    ! Create a property list to enable parallel writing, if MPI is used.
    call h5pcreate_f(H5P_DATASET_XFER_F, h5id_plist, h5err)
    call handle_hdf5_error(mpi_comm, 'h5pcreate_f', h5err)

#ifdef MPI 
    if(.not. serial_access) then 
      call h5pset_dxpl_mpio_f(h5id_plist, H5FD_MPIO_COLLECTIVE_F, h5err)
      call handle_hdf5_error(mpi_comm, 'h5pset_dxpl_mpio_f', h5err)
    end if
#endif

    ! read data from file
    call h5dread_f(h5id_dset, dataset_type, dataset_ptr, h5err, &
                    file_space_id = h5id_fspace, mem_space_id = h5id_dspace, xfer_prp = h5id_plist)
    call handle_hdf5_error(mpi_comm, 'h5dread_f', h5err)

    ! Close open objects
    call h5pclose_f(h5id_plist, h5err)
    call handle_hdf5_error(mpi_comm, 'h5pclose_f', h5err)

    call h5sclose_f(h5id_fspace, h5err)
    call handle_hdf5_error(mpi_comm, 'h5sclose_f', h5err)
    
    call h5sclose_f(h5id_dspace, h5err)
    call handle_hdf5_error(mpi_comm, 'h5sclose_f', h5err)

    call h5dclose_f(h5id_dset, h5err)
    call handle_hdf5_error(mpi_comm, 'h5dclose_f', h5err)
  
    call h5gclose_f(h5id_group, h5err)
    call handle_hdf5_error(mpi_comm, 'h5gclose_f', h5err)
#endif _HDF5_    
  end subroutine hdf5_read_dataset

end module hdf5_read  