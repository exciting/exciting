module fastBSE
  use precision, only: dp
  use constants, only: pi, zzero
  use asserts, only: assert 
  use modinput, only: input_type
  use modmpi, only: mpiinfo, terminate_if_false, distribute_loop
  use seed_generation, only: set_seed
  use math_utils, only: random_order
  use distributions, only: lorentzian
  use grid_utils, only: n_grid_diff, partial_grid, mesh_1d, linspace
  use xlapack, only: norm, diagonalize_symtridiag, matrix_multiply
  use os_utils, only: join_paths, make_directory
  use m_getunit, only: getunit

  use xfftw, only: fft_type, FFTW_FORWARD, abort_if_not_fftw3
  use unit_cell_utils, only: reciprocal_lattice, volume_parallelepiped
  use dynamic_indices, only: dynindex_type
  use xgrid, only: regular_grid_type, setup_unitcell_grid, setup_fft_grid
  
  use xisdf, only: cvt, qrcp, isdf

  use xhdf5, only: xhdf5_type, abort_if_not_hdf5
  use formatted_file_parsers, only: read_eigen_energies, read_grid_coordinates, read_QP_energies

  use bethe_salpeter_hamiltonian, only: bsh_type, vexc_isdf_kernel_type, wscr_isdf_kernel_type
  use bse_diagonal, only: setup_transition_energies
  use iterative_solver, only: lanczos_iteration
  use bse_post_processing, only: absorption_spectrum, symmetrize_quantity

  use bse_utils, only: bse_type_to_bool

  use bse_transitions, only: transition_type

  

  implicit none 

  private
  public :: fastBSE_main, fastBSE_sanity_checks

  
  ! Formated file names

  !> Brand name
  character(*), parameter :: brand_name = 'fastBSE'

  !> Name of the filet where the \(\mathbf{G}\)-vectors are read from.
  character(*), parameter :: g_grid_file = 'GQPOINTS_QMT001.OUT'
  !> File with Kohn-Sham eigen energies.
  character(*), parameter :: ip_energy_file = 'EIGVAL_QMT001.OUT'
  !> File with GW eigen energies.
  character(*), parameter :: qp_energy_file = 'EVALQP.OUT'


  ! Paths in the HDF5 file to the input datasets

  !> Name of the dataset containing the wavefunctions.
  character(*), parameter :: wfplot_dataset = 'wfplot'
  !> Name of the dataset containing the fourier transform of the screened interactions.
  character(*), parameter :: wscreened_dataset = 'screened_coulomb_non_reduced_q/interaction'
  !> Name of the dataset containing the \mathbf{q}-vectors.
  character(*), parameter :: qvectors_dataset  = 'screened_coulomb_non_reduced_q/q_vectors_cartesian'
  !> Name of the dataset containing the \mathbf{G+q}-vectors.
  character(*), parameter :: gqvectors_dataset = 'screened_coulomb_non_reduced_q/G+q_vectors_cartesian'
  !> Name of the dataset containing the number of \(\mathbf{G}\)-vectors.
  character(*), parameter :: number_of_G_per_q_dataset = 'screened_coulomb_non_reduced_q/number_of_G_per_q'
  !> Name of the dataset containing the momentum matrix.
  character(*), parameter :: pmatrix_dataset = 'pmat_xs'

  ! Names of the result datasets

  !> Name of the dataset containing \(\omega\).
  character(*), parameter :: omega_dataset = 'omega'
  !> Name of the dataset containing \(\Im(\epsilon)\).
  character(*), parameter :: eps_im_dataset = 'eps_im'
  !> Name of the dataset containing the exciton energies.
  character(*), parameter :: exciton_energies_dataset = 'exciton_eval'
  !> Name of the dataset containing the exciton eigenvectors.
  character(*), parameter :: exciton_eigenvectors_dataset = 'exciton_evec'
  !> Name of the dataset containing the Gauss quadrature energies.
  character(*), parameter :: gaussquad_energies_dataset = 'gaussquad_energies'
  !> Name of the dataset containing the Gauss quadrature weights.
  character(*), parameter :: gaussquad_weights_dataset = 'gaussquad_weights'

  contains 

  !> Main function for solving the fastBSE BSE as propesed by Henneke et al. in 
  !> *Communications in Applied Mathematics and Computational Science, (2019), 89-113, 15(1)*
  !> The function is called in [[xsmain]] and organizes the workflow as following:
  !>
  !> 1. Initialize all needed input data from the input file and priliminary calculations.
  !>    Also see [[initialize_input]].
  !>
  !> 2. Calclate ISDF for the exchange kernel and calculate the fastBSE exchange kernel.
  !>    Also see [[calculate_vexc_isdf_kernel]] and [[vexc_isdf_kernel/intialize]]
  !>
  !> 3. Calclate ISDF for the screened kernel and calculate the fastBSE screened kernel.
  !>    Also see [[calculate_wscr_isdf_kernel]] and [[wscr_isdf_kernel/intialize]]
  !>
  !> 4. Diagonalize the fastBSE BSH with the Lanczos algorithm, calculate the absorption spectrum 
  !>    and the approximation to the exciton eigen system.
  !>    Also see [[diagonalize]].
  !>
  !> 5. Write the results to an HDF5 file.
  !>    Also see [[write_results]]. 
  subroutine fastBSE_main(mpi_env, input, h5file, h5group, info_unit)
    !> MPI environment.
     type(mpiinfo), intent(inout) :: mpi_env
    !> Input file container.
    type(input_type) :: input
    !> Path to the HDF5 file. It contains input data as wave functions and
    !> is used to store the results of the run.
    character(*), intent(in) :: h5file
    !> Group in the HDF5 file with the input data and to store the results to.
    character(*), intent(in) :: h5group
    !> Unit of the info file.
    integer, intent(in) :: info_unit

    ! BSH objects

    !> Container to store the fastBSE parts of the BSH.
    type(bsh_type) :: bsh
    !> Energy differences of the bse_transitions from occupied to unoccupied states on all \(\mathbf{k}\)-points 
    !> (`[[k_grid]]`):
    !> \[
    !>  \Delta E_{vck} = E_{ck} - E_{vk}  
    !> \]
    !> Where \(E_{vk}\) and \(E_{ck}\) are the band energies of the \(v\)'th occupied and \(c\)'th unoccupied state on
    !> the \(k\)'th \(\mathbf{k}\)-point. The band energies are read in from a preliminary calculation.
    real(dp), allocatable :: transition_energies(:) ! n_occupied_bands x n_unoccupied_bands x n_kpoints
    !> Scissor shift energy, to add to the unoccupied band energies.
    type(vexc_isdf_kernel_type) :: vexc
    !> Screened interaction kernel in fastBSE form.
    type(wscr_isdf_kernel_type) :: wscr
    !> Dipole matrix elements
    complex(dp), allocatable :: dipole_matrix_elements(:, :)

    type(transition_type) :: transitions

    type(regular_grid_type) :: r_grid, g_grid

    logical :: calculate_vexc, calculate_wscr 

    integer :: r_sampling(3)
    real(dp) :: lattice(3, 3), r_offset(3), epslat, omega

    

    call bse_type_to_bool(input%xs%BSE%bsetype, calculate_vexc, calculate_wscr)

    lattice    = input%structure%crystal%basevect
    omega      = volume_parallelepiped(lattice)
    epslat     = input%structure%epslat
    r_offset   = spread(0._dp, 1, 3)
    r_sampling = input%xs%fastBSE%rsampling
    r_grid     = setup_unitcell_grid(r_sampling, r_offset, lattice, epslat)
    g_grid     = setup_fft_grid(r_grid)
    
    call setup_transitions(mpi_env, input, info_unit, h5file, h5group, transition_energies, dipole_matrix_elements, transitions)
    call bsh%set(transition_energies)


    if (calculate_vexc) then
      
      call setup_exchange_kernel(mpi_env, info_unit, h5file, h5group, g_grid, omega, transitions, vexc)
      call bsh%set(vexc)

    end if

    if (calculate_wscr) then

      call setup_screened_kernel(mpi_env, input, info_unit, h5file, h5group, g_grid, omega, transitions, wscr)
      call bsh%set(wscr)
    
    end if

    call diagonalize(mpi_env, input, info_unit, h5file, h5group, bsh, dipole_matrix_elements)

    deallocate(transition_energies)
    call bsh%finalize()

  end subroutine fastBSE_main


!================================================================================================================================


  !> Setup the transitions to be considered for fastBSE.
  subroutine setup_transitions(mpi_env, input, info_unit, h5file, h5group, transition_energies, matrix_elements, transitions)
    use fastBSE_transitions, only: h5group_transitions, h5ds_energies, h5ds_dmat, h5ds_band_idx, h5ds_mask
    !> MPI environment.      
    type(mpiinfo), intent(inout) :: mpi_env
    !> Input file container.
    type(input_type), intent(in) :: input 
    !> Unit of info file.
    integer, intent(in) :: info_unit
    !> Name of the HDF5 file to read transitions from.
    character(*), intent(in) :: h5file
    !> Name of the HDF5 group to read transitions from.
    character(*), intent(in) :: h5group
    !> Transition energies of the electrons, either from DFT or GW.
    real(dp), intent(out), allocatable :: transition_energies(:)
    !> Renormalization of the matrix elements. 
    complex(dp), intent(out), allocatable :: matrix_elements(:, :)
    !> Transition container.
    type(transition_type), intent(out) :: transitions

    integer, parameter :: iqmt = 1
    character(*), parameter :: thisname = 'setup_fastBSE_diagonal_exciting'

    integer, allocatable :: dims(:), band_idx(:, :), transition_mask_int(:)
    
    type(xhdf5_type) :: h5
    character(:), allocatable :: group

    call h5%initialize(h5file, mpi_env)
    group = join_paths(h5group, h5group_transitions)

    ! Read eigen energies
    call h5%dataset_shape(group, h5ds_energies, dims)
    allocate(transition_energies(dims(1)))
    call h5%read(group, h5ds_energies, transition_energies, [1])

    ! Read matrix elements
    allocate(matrix_elements(dims(1), 3))
    call h5%read(group, h5ds_dmat, matrix_elements, [1, 1])
    deallocate(dims)

    ! Read transitions
    call h5%dataset_shape(group, h5ds_band_idx, dims)
    allocate(band_idx(dims(1), dims(2)))
    call h5%read(group, h5ds_band_idx, band_idx, [1, 1])
    
    call assert(dims(1) == 6, &
            'Error: fastBSE: First dim of band_idx is not 6.')
    call assert(dims(2) == product(input%xs%ngridk), &
            'Error: fastBSE: Second index of band_idx is not the number of k-points.')
    deallocate(dims)
    
    ! Read allowed transitions
    call h5%dataset_shape(group, h5ds_mask, dims)
    allocate(transition_mask_int(dims(1)))
    call h5%read(group, h5ds_mask, transition_mask_int, [1])
    call h5%finalize()

    call transitions%initialize(band_idx, transition_mask_int)
  end subroutine setup_transitions


  !> Setup decomposed \(V_{x}\) as needed for fastBSE.
  subroutine setup_exchange_kernel(mpi_env, info_unit, h5file, h5group, g_grid, omega, transitions, vexc)

    ! Globals
    use fastBSE_isdf, only: h5group_isdf_vexc, h5ds_zeta, h5ds_wfplot_isdf_o, h5ds_wfplot_isdf_u
    use fastBSE_transitions, only: h5group_transitions

    !> MPI environment.      
    type(mpiinfo), intent(inout) :: mpi_env
    !> Unit of info file.
    integer, intent(in) :: info_unit
    !> Name of the HDF5 file to read ISDF from.
    character(*), intent(in) :: h5file
    !> Name of the HDF5 group to read ISDF from.
    character(*), intent(in) :: h5group
    !> \(\mathbf{G}-grid\) corresponding to the real space grid the wave functions are defined on.
    type(regular_grid_type), intent(in) :: g_grid
    !> Volume of the unit cell
    real(dp), intent(in) :: omega
    !> Transitions container.
    type(transition_type), intent(in) :: transitions
    !> Container for compressed exchange interaction kernel.
    type(vexc_isdf_kernel_type), intent(out) :: vexc
    

    character(:), allocatable :: group

    type(xhdf5_type) :: h5
    type(fft_type) :: fft

    integer :: n_r, n_o, n_u, n_isdf, r_sampling(3), i_isdf
    real(dp) :: t_start, t_end 

    integer, allocatable :: shape_u_o(:), shape_u_u(:), shape_zeta(:)
    real(dp), allocatable :: G_vecs(:, :), v_coulomb_hat(:)
    complex(dp), allocatable :: u_o_isdf(:, :), u_u_isdf(:, :), zeta(:, :)

    call timesec(t_start)

    ! Read dataset shapes
    call h5%initialize(h5file, mpi_env)
    
    group = join_paths(h5group, h5group_isdf_vexc)
    call h5%dataset_shape(group, h5ds_zeta, shape_zeta, .true.)
    call h5%dataset_shape(group, h5ds_wfplot_isdf_o, shape_u_o, .true.)
    call h5%dataset_shape(group, h5ds_wfplot_isdf_u, shape_u_u, .true.)

    ! Initialize arrays
    n_r = shape_zeta(1)
    n_isdf = shape_zeta(2)
    n_o = shape_u_o(2)
    n_u = shape_u_u(2)

    allocate(zeta(g_grid%number_of_points(), n_isdf))
    allocate(u_o_isdf(n_isdf, n_o))
    allocate(u_u_isdf(n_isdf, n_u))
    call h5%read(group, h5ds_zeta, zeta, [1, 1])
    call h5%read(group, h5ds_wfplot_isdf_o, u_o_isdf, [1, 1])
    call h5%read(group, h5ds_wfplot_isdf_u, u_u_isdf, [1, 1])
    call h5%finalize()    

    ! Fourier tranform zeta to zeta
    call fft%initialize(g_grid%sampling, FFTW_FORWARD, zeta(:, 1))
    do i_isdf=1, n_isdf
      call fft%execute(zeta(:, i_isdf))
    end do 
    call fft%delete()

    ! Setup coulomb interaction in Fourier space
    call read_grid_coordinates(g_grid_file, G_vecs) ! Read in G vectors, ordered by length and |G + k| <= gqmax
    G_vecs = G_vecs(:, 2:) ! ignore long range part
    v_coulomb_hat = omega / transitions%n_k * 4._dp * pi / sum(G_vecs ** 2, 1)

    ! Choose all G with |G + k| <= gqmax
    zeta = zeta(g_grid%composite_index(G_vecs), :)

    ! Setup exchange kernel
    call vexc%initialize(transitions, u_o_isdf, u_u_isdf, zeta, v_coulomb_hat)

    call timesec(t_end)

    write(info_unit, '(A)')        'Setup exchange interaction kernel.'
    write(info_unit, '(A, F15.6)') 'Time (s): ', t_end-t_start
    write(info_unit, *)
  end subroutine


  !> Setup decomposed screened interaction kernel for fastBSE.
  subroutine setup_screened_kernel(mpi_env, input, info_unit, h5file, h5group, g_grid, omega, transitions, wscr)
    ! h5 string names for isdf 
    use fastBSE_isdf, only: h5group_isdf_wscr_o, h5group_isdf_wscr_u, h5ds_zeta, &
                              h5ds_wfplot_isdf_o, h5ds_wfplot_isdf_u

    !> MPI environment.      
    type(mpiinfo), intent(inout) :: mpi_env
    !> Input file container.
    type(input_type), intent(in) :: input 
    !> Unit of info file.
    integer, intent(in) :: info_unit
    !> Name of the HDF5 file to read ISDF from.
    character(*), intent(in) :: h5file
    !> Name of the HDF5 group to read ISDF from.
    character(*), intent(in) :: h5group
    !> \(\mathbf{G}-grid\) corresponding to the real space grid the wave functions are defined on.
    type(regular_grid_type), intent(in) :: g_grid
    !> Volume of the unit cell
    real(dp), intent(in) :: omega
    !> Transitions container.
    type(transition_type), intent(in) :: transitions
    !> Container for compressed screened interaction kernel.
    type(wscr_isdf_kernel_type), intent(out) :: wscr

    character(:), allocatable :: group

    type(xhdf5_type) :: h5
    type(fft_type) :: fft
    type(dynindex_type) :: gq_indices

    integer :: n_r, n_k, n_isdf_o, n_isdf_u, r_sampling(3), n_gq_max, k_sampling(3), i_isdf, n_o, n_u
    integer :: ik_first, ik_last, n_k_local
    real(dp) :: t_start, t_end

    integer, allocatable :: shape_u_o(:), shape_u_u(:), shape_zeta_o(:), shape_zeta_u(:)
    real(dp), allocatable :: G_vecs(:, :), q_vecs(:, :), gq_coordinates(:, :, :)
    complex(dp), allocatable :: u_o_isdf(:, :), u_u_isdf(:, :), zeta_o(:, :), zeta_u(:, :)
    complex(dp), allocatable :: w_screened_hat(:, :, :)
    integer, allocatable :: n_g_per_q(:), index_matrix(:, :)
    
   call timesec(t_start)

   call h5%initialize(h5file, mpi_env)

   ! Read dataset shapes
   group = join_paths(h5group, h5group_isdf_wscr_o)
   call h5%dataset_shape(group, h5ds_zeta, shape_zeta_o, .true.)
   call h5%dataset_shape(group, h5ds_wfplot_isdf_o, shape_u_o, .true.)

   group = join_paths(h5group, h5group_isdf_wscr_u)
   call h5%dataset_shape(group, h5ds_zeta, shape_zeta_u, .true.)
   call h5%dataset_shape(group, h5ds_wfplot_isdf_u, shape_u_u, .true.)
   

   ! Initialize array sizes  
   n_isdf_o   = shape_zeta_o(2)
   n_isdf_u   = shape_zeta_u(2)
   n_r        = shape_zeta_o(1)
   
   n_u    = shape_u_u(2)
   n_o    = shape_u_o(2)
   k_sampling = input%xs%ngridk
   n_k    = product(k_sampling)

   ! Setup G+q grid
   allocate(n_g_per_q(n_k))
   call h5%read(h5group, number_of_G_per_q_dataset, n_g_per_q, [1])
   n_gq_max = maxval(n_g_per_q)

   allocate(gq_coordinates(3, n_gq_max, n_k))
   allocate(q_vecs(3, n_k))
   call h5%read(h5group, gqvectors_dataset, gq_coordinates, [1, 1, 1])
   call h5%read(h5group, qvectors_dataset, q_vecs, [1, 1])

   
   G_vecs = reshape(gq_coordinates - spread(q_vecs, 2, n_gq_max), [3, n_gq_max * n_k])
   index_matrix = reshape(g_grid%composite_index(G_vecs), [n_gq_max, n_k])
   call gq_indices%init(index_matrix, n_g_per_q)

   ! Read screened interaction 
   allocate(w_screened_hat(n_gq_max, n_gq_max, n_k))
   call h5%read(h5group, wscreened_dataset, w_screened_hat, [1, 1, 1])

   ! Read ISDF datasets
   group = join_paths(h5group, h5group_isdf_wscr_o)
   allocate(zeta_o(n_r, n_isdf_o))
   allocate(u_o_isdf(n_isdf_o, n_o))
   call h5%read(group, h5ds_zeta, zeta_o, [1, 1])
   call h5%read(group, h5ds_wfplot_isdf_o, u_o_isdf, [1, 1, 1])

   group = join_paths(h5group, h5group_isdf_wscr_u)
   allocate(zeta_u(n_r, n_isdf_u))
   allocate(u_u_isdf(n_isdf_u, n_u))
   call h5%read(group, h5ds_zeta, zeta_u, [1, 1])
   call h5%read(group, h5ds_wfplot_isdf_u, u_u_isdf, [1, 1, 1])
   call h5%finalize()

   ! Fourier tranform zeta to zeta
   call fft%initialize(g_grid%sampling, FFTW_FORWARD)
   do i_isdf=1, n_isdf_o
     call fft%execute(zeta_o(:, i_isdf))
   end do 

   do i_isdf=1, n_isdf_u
     call fft%execute(zeta_u(:, i_isdf))
   end do 
   call fft%delete()

   ! Setup screened kernel
   call wscr%intialize(mpi_env, transitions, u_o_isdf, u_u_isdf, zeta_o, zeta_u, omega/n_k, &
           w_screened_hat, gq_indices, g_grid, k_sampling)

   call timesec(t_end)
   
   write(info_unit, '(A)')        'Setup screened interaction kernel.'
   write(info_unit, '(A, F15.6)') 'Time (s): ', t_end-t_start
   write(info_unit, *)
  end subroutine 

  !> Diagonlaize the fastBSE BSH and calculate the results:
  !>
  !> 1. Diagonalize the decompoed BSH with a Lanczos iteration.
  !>
  !> 2. Calculate the absorption spectrum.
  !>
  !> 3. Calculate the approximation to the exciton eigen energies and vectors.
  subroutine diagonalize(mpi_env, input, info_unit, h5file, h5group, bsh, dipole_matrix_elements)
    use unit_conversion, only: hartree_to_ev
    !> MPI environment.      
    type(mpiinfo), intent(inout) :: mpi_env
    !> Input file container.
    type(input_type), intent(in) :: input 
    !> Unit of info file.
    integer, intent(in) :: info_unit
    !> Name of the HDF5 file to write results to.
    character(*), intent(in) :: h5file
    !> Name of the HDF5 group to write results to.
    character(*), intent(in) :: h5group
    !> fastBSE Bethe-Salpeter Hamiltonian.
    type(bsh_type), intent(inout) :: bsh
    !> Optical absorption vector.
    complex(dp), intent(in) :: dipole_matrix_elements(:, :)

    integer :: lanczosmaxits
    real(dp) :: scaling, broadening
    complex(dp), allocatable :: exciton_evec(:, :, :)
    real(dp), allocatable ::  omega(:), eps_im(:, :), exciton_eval(:, :), gaussquad_energies(:, :), gaussquad_weights(:, :)

    type(xhdf5_type) :: h5
    real(dp) :: t_start, t_end

    character(:), allocatable :: bse_type, result_group, text_file
    integer :: i_dim, iterations, n_k, n_o, n_u, n_transitions, n_omega_points, ierr
    real(dp) :: omega_vol, lattice(3, 3), omega_intervall(2), energy_conversion
    real(dp), allocatable ::  alpha(:), beta(:), eigenvalues(:), eigenvectors(:, :), tridiag_vec(:, :)
    complex(dp), allocatable :: Q_k(:, :)

    call timesec(t_start)

    ! Intitialize parameters from input file
    lattice         = input%structure%crystal%basevect
    omega_vol       = volume_parallelepiped(lattice)
    n_k             = product(input%xs%ngridk)
    n_o             = input%xs%BSE%nstlbse(2) - input%xs%BSE%nstlbse(1) + 1
    n_u             = input%xs%BSE%nstlbse(4) - input%xs%BSE%nstlbse(3) + 1
    n_transitions   = size(dipole_matrix_elements, 1)
    broadening      = input%xs%broad
    scaling         = 8 * pi**2 / omega_vol / n_k
    omega_intervall = input%xs%energywindow%intv
    n_omega_points  = input%xs%energywindow%points
    bse_type        = trim(adjustl(input%xs%BSE%bsetype))
    lanczosmaxits       = min(input%xs%fastBSE%lanczosmaxits, n_transitions)
    
    ! Setup omega
    omega = linspace(omega_intervall, n_omega_points)

    ! Initialize arrays
    allocate(eigenvectors(2 * lanczosmaxits - 1, 2 * lanczosmaxits - 1))
    allocate(tridiag_vec(lanczosmaxits, lanczosmaxits))
    allocate(gaussquad_energies(2 * lanczosmaxits - 1, 3), source = 0._dp)
    allocate(gaussquad_weights(2 * lanczosmaxits - 1, 3), source = 0._dp)
    allocate(eps_im(n_omega_points, 3), source = 0._dp)
    allocate(exciton_eval(lanczosmaxits, 3), source = 0._dp)
    allocate(exciton_evec(n_transitions, lanczosmaxits, 3), source = zzero)
    
    energy_conversion = 1._dp
    if(input%xs%tevout) energy_conversion = hartree_to_ev
    

    do i_dim=1, 3

      ! Run Lanczos iteration
      call lanczos_iteration(lanczosmaxits, bsh_times_vector, dipole_matrix_elements(:, i_dim), alpha, beta, Q_k)
      iterations = size(alpha)

      ! Calculate eigenvalues and oscillator strengths for the spectrum via gauss quadrature
      call diagonalize_symtridiag( &
              [alpha(mesh_1d(1, iterations)), alpha(mesh_1d(iterations - 1, 1))], &
              [beta(mesh_1d(1, iterations)), beta(mesh_1d(iterations - 2, 1))], &
              eigenvalues, eigenvectors &
            )
      gaussquad_energies(:, i_dim) = eigenvalues
      gaussquad_weights(:, i_dim) = eigenvectors(1, :) ** 2
      

      ! Calculate spectrum
      call absorption_spectrum( &
              pack(eigenvalues, eigenvalues >= 0._dp), &
              pack(eigenvectors(1, :) ** 2, eigenvalues >= 0._dp), &
              omega, &
              broadening, &
              lorentzian, &
              eps_im(:, i_dim) &
            )
      eps_im(:, i_dim) = scaling * norm(dipole_matrix_elements(:, i_dim))**2 * eps_im(:, i_dim)


      ! Calculate approximation to the eigen system by calculating the ritz values and vectors
      call diagonalize_symtridiag(alpha, beta(: iterations - 1), eigenvalues, tridiag_vec)
      exciton_eval(:, i_dim) = eigenvalues
      call matrix_multiply(Q_k, tridiag_vec, exciton_evec(:, :, i_dim))

      deallocate(alpha, beta, Q_k)
    end do

    ! Symmetrize quantities
    if(.not. input%xs%BSE%nosymspec) then
      call symmetrize_quantity(eps_im)
      call symmetrize_quantity(gaussquad_energies)
      call symmetrize_quantity(gaussquad_weights)
    end if


    

    result_group = brand_name // '_' // bse_type 
    call h5%initialize(h5file, mpi_env)
    call h5%initialize_group(h5group, result_group)
    result_group = trim(adjustl(join_paths(h5group, result_group)))

    call h5%write(result_group, omega_dataset, omega, [1], shape(omega))
    call h5%write(result_group, eps_im_dataset, eps_im, [1, 1], shape(eps_im))
    call h5%write(result_group, exciton_energies_dataset, exciton_eval, [1, 1], shape(exciton_eval))
    call h5%write(result_group, exciton_eigenvectors_dataset, exciton_evec, [1, 1, 1], shape(exciton_evec))
    call h5%write(result_group, gaussquad_energies_dataset, gaussquad_energies, [1, 1], shape(gaussquad_energies))
    call h5%write(result_group, gaussquad_weights_dataset, gaussquad_weights, [1, 1], shape(gaussquad_weights))
    call h5%finalize()

    text_file = brand_name // '_' // bse_type // '_' // eps_im_dataset // '.out'
    call write_eps_im_textfile(text_file , omega, eps_im, energy_conversion)

    deallocate(omega, eps_im, exciton_eval, exciton_evec, gaussquad_energies, gaussquad_weights, eigenvalues, eigenvectors, tridiag_vec)

    call timesec(t_end)

    write(info_unit, '(A)')        'Diagonalization done.'
    write(info_unit, '(A, I8)')    'Lanczos iterations done: ', iterations 
    write(info_unit, '(A, F15.6)') 'Time (s): ', t_end - t_start
    write(info_unit, *)    

    contains 

    subroutine bsh_times_vector(vector_in, vector_out)
      complex(dp), intent(in) :: vector_in(:)
      complex(dp), intent(out) :: vector_out(:)
      call bsh%times_vector(vector_in, vector_out)
    end subroutine

    !> Write the imaginary part of the diagonal of \(\epsilon^M\) to a text file.
    !> The first column contains \(\omega\), the following \(\epsilon^M_{11}), \(\epsilon^M_{22}) and \(\epsilon^M_{33}).
    subroutine write_eps_im_textfile(fname, omega, eps_im, energy_conversion)
      !> Name of the file
      character(*), intent(in) :: fname
      !> Energy grid for \(\omega\)
      real(dp), intent(in) :: omega(:)
      !> Imaginary part of the diagonal of \(\epsilon^M\) with dimensions stored columnwise.
      real(dp), intent(in) :: eps_im(:, :)
      !> Energy scaling. As usual 1.0 means the energy is in Hartree.
      real(dp), intent(in) :: energy_conversion

      integer :: n_omega, i_omega, unit 


      n_omega = size(omega)

      call getunit(unit)

      open(unit, file=fname, form='formatted', action='write', status='replace')

      write(unit, '(A)') '# fastBSE imaginary macroscopic dielectric function'
      write(unit, '(A)') '# '
      write(unit, '(A, f12.6, A)')'# Energy scaling: ', energy_conversion, ' Hartree'
      write(unit, '(A, f12.6, A, f12.6, A)')'# Broadening: ', broadening, ' x ', energy_conversion, ' Hartree'
      write(unit, '(A, A22, 1x, A23, 1x, A23, 1x, A23)') '#',  'omega', 'oc11', 'oc22', 'oc33'
      write(unit, '(SP, E23.16, 1x, E23.16, 1x, E23.16, 1x, E23.16)') &
              (omega(i_omega) * energy_conversion, eps_im(i_omega, 1), eps_im(i_omega, 2), eps_im(i_omega, 3), i_omega=1, n_omega)
      
      close(unit)
    end subroutine

  end subroutine

  !> Terminate exciting if the input file is not compatible with `fastBSE`.
  subroutine fastBSE_sanity_checks(mpi_env, input)
    use modmpi, only: terminate_mpi_env
    !> MPI environment to terminate.
    type(mpiinfo), intent(inout) :: mpi_env
    !> Input file container to check.
    type(input_type), intent(in) :: input 

    call abort_if_not_fftw3(mpi_env, "Error(fastBSE): exciting needs to be linked to FFTW3 for running fastBSE.")
    call abort_if_not_hdf5(mpi_env, "Error(fastBSE): exciting needs to be compiled with HDF5 to run fastBSE module.")

    if(input%xs%BSE%coupling) then 
      call terminate_mpi_env(mpi_env, 'Error(fastBSE): fastBSE only supports TDA.')
    end if 

    if(input%xs%BSE%xas) then
      call terminate_mpi_env(mpi_env, 'Error(fastBSE): fastBSE only supports valence excitations.')
    end if

    if(associated(input%xs%excitonPlot)) then
      call terminate_mpi_env(mpi_env, 'Error(fastBSE): fastBSE does not support the element excitonPlot.')
    end if

    if(associated(input%xs%writekpathweights)) then
      call terminate_mpi_env(mpi_env, 'Error(fastBSE): fastBSE does not support the element writekpathweights.')
    end if

    if(associated(input%xs%excitonPlot)) then
      call terminate_mpi_env(mpi_env, 'Error(fastBSE): fastBSE does not support the element writeexcitons.')
    end if

  end subroutine

end module fastBSE   
