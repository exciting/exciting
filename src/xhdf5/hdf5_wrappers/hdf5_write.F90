!> Wrapper for writing any dataset to an HDF5 file
module hdf5_write
  use iso_c_binding, only: c_ptr
  use hdf5_utils
  use asserts, only: assert

  implicit none

  private
  public :: hdf5_write_dataset


  contains 
  

  !> Write a dataset or chunk of an dataset of any type and shape to file. 
  subroutine hdf5_write_dataset(h5id_file, h5path, dataset, data_type_id, dataset_chunk_ptr, dataset_rank, dataset_chunk_shape, dataset_shape, offset, serial_access)
    !> Identifier of the file or group, used by HDF5.
    integer(hdf5_id), intent(in) :: h5id_file
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Type id of the dataset.
    integer(hdf5_id), intent(in) :: data_type_id
    !> Pointer to the first element of the dataset to be written.
    type(c_ptr), intent(inout) :: dataset_chunk_ptr
    !> Rank of the dataset.
    integer, intent(in) :: dataset_rank
    !> Shape of the local dataset chunk.
    integer(hdf5_size), intent(in) :: dataset_chunk_shape(dataset_rank)
    !> Shape of the total dataset.
    integer(hdf5_size), intent(in) :: dataset_shape(dataset_rank)
    !> Offset in the total dataset of the local dataset chunk.
    integer(hdf5_ssize) , intent(in) :: offset(dataset_rank)
    !> Set to `.true.` if only serial access is possible.
    logical, intent(in) :: serial_access

#ifdef XHDF5
    integer(hdf5_id) :: h5id_group, h5id_fspace, h5id_dset, h5id_dspace, h5id_plist
    integer :: h5err 

    if(serial_access) then
      call assert(all(offset == 0), 'serial_access is .true. but all(offset == 0) is .false.')
    end if 

    ! Open group to write the dataset in
    call h5gopen_f(h5id_file, h5path, h5id_group, h5err)
    call handle_hdf5_error('h5gopen_f' ,h5err)

    ! Create a dataset for the whole array
    call h5screate_simple_f(dataset_rank, dataset_shape, h5id_fspace, h5err)
    call handle_hdf5_error('h5screate_simple_f', h5err)

    call h5dcreate_f(h5id_group, dataset, data_type_id, h5id_fspace, h5id_dset, h5err)
    call handle_hdf5_error('h5dcreate', h5err)

    call h5sclose_f(h5id_fspace, h5err)
    call handle_hdf5_error('h5sclose_f', h5err)

    ! Create a dataset for the local chunk of the array and select the hyperslap to write into
    call h5screate_simple_f(dataset_rank, dataset_chunk_shape, h5id_dspace, h5err)
    call handle_hdf5_error('h5screate_simple_f', h5err)

    call h5dget_space_f(h5id_dset, h5id_fspace, h5err)
    call handle_hdf5_error('h5dget_space_f', h5err)

    call h5sselect_hyperslab_f(h5id_fspace, H5S_SELECT_SET_F, offset, dataset_chunk_shape, h5err)
    call handle_hdf5_error('h5sselect_hyperslab_f', h5err)

    ! Create a property list to enable parallel writing, if MPI is used.
    call h5pcreate_f(H5P_DATASET_XFER_F, h5id_plist, h5err)
    call handle_hdf5_error('h5pcreate_f', h5err)

#ifdef MPI
    if(.not. serial_access) then 
      call h5pset_dxpl_mpio_f(h5id_plist, H5FD_MPIO_COLLECTIVE_F, h5err)
      call handle_hdf5_error('h5pset_dxpl_mpio_f', h5err)
    end if
#endif

    ! Write data to file
    call h5dwrite_f(h5id_dset, data_type_id, dataset_chunk_ptr, h5err, &
                    file_space_id = h5id_fspace, mem_space_id = h5id_dspace, xfer_prp = h5id_plist)
    call handle_hdf5_error('h5dwrite_f', h5err)

    ! Close open objects
    call h5pclose_f(h5id_plist, h5err)
    call handle_hdf5_error('h5pclose_f', h5err)

    call h5sclose_f(h5id_fspace, h5err)
    call handle_hdf5_error('h5sclose_f', h5err)
    
    call h5sclose_f(h5id_dspace, h5err)
    call handle_hdf5_error('h5sclose_f', h5err)

    call h5dclose_f(h5id_dset, h5err)
    call handle_hdf5_error('h5dclose_f', h5err)
  
    call h5gclose_f(h5id_group, h5err)
    call handle_hdf5_error('h5gclose_f', h5err)
#endif    
  end subroutine hdf5_write_dataset

end module hdf5_write