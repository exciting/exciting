module fastBSE_file_strings

  !> Name of the filet where the \(\mathbf{G}\)-vectors are read from.
  character(*), parameter :: g_grid_file = 'GQPOINTS_QMT001.OUT'

  !> Name of the group that holds the results
  character(*), parameter :: result_group = 'fastBSE_results'
  !> Name of the dataset containing \(\omega\).
  character(*), parameter :: omega_dataset = 'omega'
  !> Name of the dataset containing \(\Im(\epsilon)\).
  character(*), parameter :: absspec_dataset = 'absorption_spectrum'
  !> Name of the dataset containing the Gauss quadrature energies.
  character(*), parameter :: exc_evals_dataset = 'exciton_energies'
  !> Name of the dataset containing the Gauss quadrature weights.
  character(*), parameter :: oscstr_dataset = 'oscillator_strengths'
  !> Name of the dataset containing the exciton evecs_tridiag.
  character(*), parameter :: exc_evecs_dataset = 'exciton_eigenvectors'
  !> Name of the dataset containing the independent particles band gap.
  character(*), parameter :: ip_bandgap_dataset = 'ip_band_gap'
  !> Name of the dataset containing the number of valid excitons per cartesian dimension.
  character(*), parameter :: n_exc_dataset = 'n_excitons'
  !> Name of the dataset containing the number lanczos iterations per cartesian dimension.
  character(*), parameter :: n_its_dataset = 'n_lanczos_iterations'

  !> Text file name for the absorption spectrum file.
  character(*), parameter :: fname_calculate_absorption_spectrum = 'fastBSE_absorption_spectrum.out'
  !> Text file name for the exciton eneriges.
  character(*), parameter :: fname_exciton_energies = 'fastBSE_exciton_energies.out'
  !> Text file name for the oscillator strength from gauss quadrature text file.
  character(*), parameter :: fname_oscillator_strengths = 'fastBSE_oscillator_strengths.out'

  !> Group name for the isdf objects for fastBSE
  character(*), parameter :: isdf_group = "ISDF"
  !> Group name for the isdf objects to compose the fastBSE exchange potential
  character(*), parameter :: vexc_ou_group = "vexc_ou"
  !> Group name for the isdf objects corresponding to the occupied wave functions
  !> to compose the fastBSE screened potential
  character(*), parameter :: wscr_oo_group = "wscr_oo"
  !> Group name for the isdf objects corresponding to the unoccupied wave functions
  !> to compose the fastBSE screened potential
  character(*), parameter :: wscr_uu_group = "wscr_uu"
  !> Dataset name for the ISDF coordinates 
  character(*), parameter :: rspace_coordinates_dataset = "r_cartesian"
  !> Dataset name for the ISDF interpolation point indices in the full grid
  character(*), parameter :: isdf_indices_dataset = "r_isdf_indices"
  !> Dataset name for the ISDF coefficients
  character(*), parameter :: zeta_dataset = "zeta"
  !> Dataset name of the periodic part of the occupied states wave functions evaluated on the 
  !> ISDF coordinates
  character(*), parameter :: u_o_isdf_dataset = "u_o_isdf"
  !> Dataset name of the periodic part of the unoccupied states wave functions evaluated on the 
  !> ISDF coordinates
  character(*), parameter :: u_u_isdf_dataset = "u_u_isdf"

    !> Group name in hdf5 file to save transitions for fastBSE.
  character(*), parameter, public :: groundstate_properties_group = "groundstate_properties"
  !> Dataset name for the band index lookup table
  character(*), parameter, public :: band_index_dataset = "band_index_lookup_table"
  !> Dataset name for the list of bands
  character(*), parameter, public :: band_list_dataset = 'band_list'
  !> Dataset name for the eigen energies.
  character(*), parameter, public :: eigen_energies_dataset = "eigen_energies"
  !> Dataset name for the list of \(\mathbf k\)-points.
  character(*), parameter, public :: k_list_dataset = 'k_list'
  !> Dataset name for lattice vectors.
  character(*), parameter, public :: lattice_vectors_dataset = 'lattice_vectors'
  !> Dataset name for the mask of the transitions
  character(*), parameter, public :: transition_mask_dataset = "mask"
  !> Dataset name for the matrix elements
  character(*), parameter, public :: matrix_elements_dataset = "matrix_elements"
  !> Dataset name for the number of \(\mathbf k\)-points per dimension
  character(*), parameter, public :: ngridk_dataset = 'ngridk'
  !> Dataset name for the \(\mathbf r\)-sampling
  character(*), parameter, public :: ngridr_dataset = 'ngridr'
  !> Dataset name for the transition energies.
  character(*), parameter, public :: transition_energies_dataset = "transition_energies"
  !> Dataset name for the periodic part of the wavefunctions
  character(*), parameter, public :: u_dataset = 'u'
  !> Dataset name for the unoccupied and occupied band limits per \(\mathbf k\)-point
  character(*), parameter, public :: uo_limits_dataset = "uo_limits"

end module fastBSE_file_strings