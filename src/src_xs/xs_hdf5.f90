module xs_hdf5

  !> HDF5 file for all data related to the xs part.
  character(*), parameter :: h5file_xs = 'xs.h5'
  !> Dataset name for wave function plot (in real space).
  character(*), parameter :: h5ds_wfplot = 'wfplot' 
  !> Dataset name for \(\mathbf{p}\)-matrix elements.
  character(*), parameter :: h5ds_pmat = 'pmat'

  !> HDF5 file for all data related to the screening part of the xs part.
  character(*), parameter :: h5file_screening = 'xs_screening.h5'
  
  
  character(*), parameter :: h5file_results = 'xs_results.h5'

end module xs_hdf5