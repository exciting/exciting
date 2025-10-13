module fastBSE
  use precision, only: dp, i32
  use constants, only: pi, zzero
  use asserts, only: assert 
  use modinput, only: input_type
  use modmpi, only: mpiinfo, terminate_if_false, distribute_loop, warn_if_false
  use seed_generation, only: set_seed
  use math_utils, only: random_order
  use distributions, only: lorentzian
  use grid_utils, only: n_grid_diff, partial_grid, mesh_1d, linspace
  use xlapack, only: norm, diagonalize_symtridiag, matrix_multiply
  use os_utils, only: join_paths
  use unit_conversion, only: hartree_to_ev
  use xstring
  use xfftw, only: fft_type, FFTW_FORWARD, abort_if_not_fftw3
  use unit_cell_utils, only: reciprocal_lattice, volume_parallelepiped
  use dynamic_indices, only: dynindex_type
  use xgrid, only: regular_grid_type, setup_unitcell_grid, setup_fft_grid
  use xhdf5, only: xhdf5_type, abort_if_not_hdf5
  use formatted_file_parsers, only: read_eigen_energies, read_grid_coordinates, read_QP_energies
  use bethe_salpeter_hamiltonian, only: bsh_type, vexc_isdf_kernel_type, wscr_isdf_kernel_type
  use bse_diagonal, only: setup_transition_energies
  use iterative_solver, only: lanczos
  use bse_post_processing, only: calculate_absorption_spectrum, setup_symmetric_matrix
  use bse_utils, only: bse_type_to_bool
  use fastBSE_groundstate_properties, only: read_transitions_hdf5
  use fastBSE_isdf, only: read_isdf_hdf5
  use write_screening, only: read_screened_coulomb_hdf5
  use bse_transitions, only: transition_type
  use fastBSE_file_strings


  implicit none 


  private
  public :: fastBSE_main, fastBSE_sanity_checks, fastBSE_human_readable_output, read_fastBSE_excitons_hdf5


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
  subroutine fastBSE_main(mpi_env, input, h5file, h5path, info_unit)
    !> MPI environment.
     type(mpiinfo), intent(inout) :: mpi_env
    !> Input file container.
    type(input_type) :: input
    !> Path to the HDF5 file. It contains input data as wave functions and
    !> is used to store the results of the run.
    character(*), intent(in) :: h5file
    !> Group in the HDF5 file with the input data and to store the results to.
    character(*), intent(in) :: h5path
    !> Unit of the info file.
    integer(i32), intent(in) :: info_unit

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
    !> Exchange interaction kernel in fastBSE form.
    type(vexc_isdf_kernel_type) :: vexc
    !> Screened interaction kernel in fastBSE form.
    type(wscr_isdf_kernel_type) :: wscr
    !> Dipole matrix elements
    complex(dp), allocatable :: dipole_matrix_elements(:, :)

    type(transition_type) :: transitions

    type(regular_grid_type) :: r_grid, g_grid

    logical :: calculate_vexc, calculate_wscr 

    integer(i32) :: r_sampling(3)
    real(dp) :: lattice(3, 3), r_offset(3), epslat, omega

    

    call bse_type_to_bool(input%xs%BSE%bsetype, calculate_vexc, calculate_wscr)

    lattice    = input%structure%crystal%basevect
    omega      = volume_parallelepiped(lattice)
    epslat     = input%structure%epslat
    r_offset   = spread(0._dp, 1, 3)
    r_sampling = input%xs%fastBSE%ngridr
    r_grid     = setup_unitcell_grid(r_sampling, r_offset, lattice, epslat)
    g_grid     = setup_fft_grid(r_grid)
    
    call setup_transitions(mpi_env, input, info_unit, h5file, h5path, transition_energies, dipole_matrix_elements, transitions)
    call bsh%set(transition_energies)

    if (calculate_vexc) then
      call setup_exchange_kernel(mpi_env, info_unit, h5file, h5path, g_grid, omega, transitions, vexc)
      call bsh%set(vexc)
    end if

    if (calculate_wscr) then
      call setup_screened_kernel(mpi_env, input, info_unit, h5file, h5path, g_grid, omega, transitions, wscr)
      call bsh%set(wscr)
    end if

    call diagonalize(mpi_env, input, info_unit, h5file, h5path, bsh, dipole_matrix_elements)

    deallocate(transition_energies)
    call bsh%finalize()
  end subroutine fastBSE_main


!================================================================================================================================


  !> Setup the transitions to be considered for fastBSE.
  subroutine setup_transitions(mpi_env, input, info_unit, h5file, h5path, transition_energies, matrix_elements, transitions)
    !> MPI environment.      
    type(mpiinfo), intent(inout) :: mpi_env
    !> Input file container.
    type(input_type), intent(in) :: input 
    !> Unit of info file.
    integer(i32), intent(in) :: info_unit
    !> Name of the HDF5 file to read transitions from.
    character(*), intent(in) :: h5file
    !> Name of the HDF5 group to read transitions from.
    character(*), intent(in) :: h5path
    !> Transition energies of the electrons, either from DFT or GW.
    real(dp), intent(out), allocatable :: transition_energies(:)
    !> Renormalization of the matrix elements. 
    complex(dp), intent(out), allocatable :: matrix_elements(:, :)
    !> Transition container.
    type(transition_type), intent(out) :: transitions

   
    integer(i32), allocatable :: dims(:), band_indices(:, :), transition_mask(:)
    
    type(xhdf5_type) :: h5
    character(:), allocatable :: group

    call h5%initialize(h5file, mpi_env)
    call read_transitions_hdf5(mpi_env, h5, h5path, transition_energies, matrix_elements, band_indices, transition_mask)
    call h5%finalize()

    call transitions%initialize(band_indices, transition_mask)
  end subroutine setup_transitions


  !> Setup decomposed \(V_{x}\) as needed for fastBSE.
  subroutine setup_exchange_kernel(mpi_env, info_unit, h5file, h5path, g_grid, omega, transitions, vexc)
    !> MPI environment.      
    type(mpiinfo), intent(inout) :: mpi_env
    !> Unit of info file.
    integer(i32), intent(in) :: info_unit
    !> Name of the HDF5 file to read ISDF from.
    character(*), intent(in) :: h5file
    !> Name of the HDF5 group to read ISDF from.
    character(*), intent(in) :: h5path
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

    integer(i32) :: i_isdf, n_isdf, i_G_first
    real(dp) :: t_start, t_end 

    !> \(\mathbf G\) vectors to be considered for calculating the exhanges kernel.
    real(dp), allocatable :: G_vecs(:, :)
    !> Coulomb inter action in reciprocal space.
    real(dp), allocatable ::v_coulomb_hat(:)
    !> Occupied and unoccupied wave functions, evaluated at the inderpolation points (ISDF).
    complex(dp), allocatable :: u_o_isdf(:, :), u_u_isdf(:, :)
    !> Interpolation coefficients (ISDF).
    complex(dp), allocatable :: zeta(:, :)

    call timesec(t_start)

    ! Read dataset shapes
    call h5%initialize(h5file, mpi_env)
    call read_isdf_hdf5(mpi_env, h5, h5path, zeta, u_o_isdf=u_o_isdf, u_u_isdf=u_u_isdf)
    call h5%finalize() 
    
    n_isdf = size(zeta, 2)
    
    ! Fourier transform the interpolation coefficients
    call fft%initialize(g_grid%sampling, FFTW_FORWARD, zeta(:, 1))
    do i_isdf=1, n_isdf
      call fft%execute(zeta(:, i_isdf))
    end do 
    call fft%delete()

    ! Setup coulomb interaction in Fourier space
    call read_grid_coordinates(g_grid_file, G_vecs) ! Read in G vectors, ordered by length and |G + k| <= gqmax
    call terminate_if_false(mpi_env, size(G_vecs, 2) > 1, &
            'fastBSE_main: G grid only contains the Gamma point. fastBSE needs more G vectors. Increase gqmax.')
    
    i_G_first = 2 ! ignore long range part and ommit the Gamma point
    v_coulomb_hat = omega / transitions%n_k * 4._dp * pi / sum(G_vecs(:, i_G_first:) ** 2, 1) ! Test OMP reduction here

    ! Choose all G with |G + k| <= gqmax
    zeta = zeta(g_grid%composite_index(G_vecs(:, i_G_first:)), :)
   
    ! Setup exchange kernel
    call vexc%initialize(transitions, u_o_isdf, u_u_isdf, zeta, v_coulomb_hat)

    call timesec(t_end)

    write(info_unit, '(A)')        'Setup exchange interaction kernel.'
    write(info_unit, '(A, F15.6)') 'Time (s): ', t_end-t_start
    write(info_unit, *)
  end subroutine


  !> Setup decomposed screened interaction kernel for fastBSE.
  subroutine setup_screened_kernel(mpi_env, input, info_unit, h5file, h5path, g_grid, omega, transitions, wscr)
    !> MPI environment.      
    type(mpiinfo), intent(inout) :: mpi_env
    !> Input file container.
    type(input_type), intent(in) :: input 
    !> Unit of info file.
    integer(i32), intent(in) :: info_unit
    !> Name of the HDF5 file to read ISDF from.
    character(*), intent(in) :: h5file
    !> Name of the HDF5 group to read ISDF from.
    character(*), intent(in) :: h5path
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

    integer(i32) :: k_sampling(3), i_isdf, n_k, n_G_per_q_max, n_isdf_o, n_isdf_u, iq
    real(dp) :: t_start, t_end

    

    !> \(\mathbf G\), \(\mathbf q\) and \(\mathbf{G+q}\) vectors to be considered for calculating the screened kernel.
    real(dp), allocatable :: G_vecs(:, :), q(:, :), G_plus_q(:, :, :)
    !> Occupied and unoccupied wave functions, evaluated at the inderpolation points (ISDF).
    complex(dp), allocatable :: u_o_isdf(:, :), u_u_isdf(:, :)
    !> Interpolation coefficients (ISDF) for occupied and unoccupied pairing.
    complex(dp), allocatable :: zeta_o(:, :), zeta_u(:, :)
    !> Screened coulomb interaction
    complex(dp), allocatable :: w(:, :, :) ! n_G x n_G x n_q

    integer(i32), allocatable :: n_G_per_q(:), index_matrix(:, :)
    
    call timesec(t_start)
 
    call h5%initialize(h5file, mpi_env)
    call read_isdf_hdf5(mpi_env, h5, h5path, zeta_o, u_o_isdf = u_o_isdf)
    call read_isdf_hdf5(mpi_env, h5, h5path, zeta_u, u_u_isdf = u_u_isdf)
    call read_screened_coulomb_hdf5(mpi_env, h5, h5path, w, q, G_plus_q, n_G_per_q, k_sampling)
    call h5%finalize()
 
    n_k = product(k_sampling)
    n_G_per_q_max = maxval(n_G_per_q)
    n_isdf_o = size(zeta_o, 2)
    n_isdf_u = size(zeta_u, 2)
 
    ! Setup G+q grid
    G_vecs = reshape(G_plus_q - spread(q, 2, n_G_per_q_max), [3, n_G_per_q_max * n_k])
    index_matrix = reshape(g_grid%composite_index(G_vecs), [n_G_per_q_max, n_k])
    call gq_indices%init(index_matrix, n_G_per_q)
 
    ! Fourier transform zetas
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
            w, gq_indices, g_grid, k_sampling)
 
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
  subroutine diagonalize(mpi_env, input, info_unit, h5file, h5path, bsh, dipole_matrix_elements)
    !> MPI environment.      
    type(mpiinfo), intent(inout) :: mpi_env
    !> Input file container.
    type(input_type), intent(in) :: input 
    !> Unit of info file.
    integer(i32), intent(in) :: info_unit
    !> Name of the HDF5 file to write results to.
    character(*), intent(in) :: h5file
    !> Name of the HDF5 group to write results to.
    character(*), intent(in) :: h5path
    !> fastBSE Bethe-Salpeter Hamiltonian.
    type(bsh_type), intent(inout) :: bsh
    !> Optical absorption vector.
    complex(dp), intent(in) :: dipole_matrix_elements(:, :)

    integer(i32) :: nlanczos
    real(dp) :: scaling, broadening

    type(xhdf5_type) :: h5
    character(:), allocatable :: group
    real(dp) :: t_start, t_end

    logical :: save_exc_evecs
    logical, allocatable :: eval_mask(:), evec_mask(:), Q_lanczos_mask(:)
    character(:), allocatable :: bse_type, text_file
    integer(i32) :: i_dim, n_its, n_exc, n_its_save(3), n_exc_save(3), n_k, n_o, n_u, n_transitions, n_omega, ierr, i_exc
    real(dp) :: omega_vol, lattice(3, 3), omega_intervall(2), energy_conversion, symmetric_matrix(3, 3)

    real(dp), allocatable :: omega(:), absspec(:, :), exc_evals_save(:, :), oscstr(:), oscstr_save(:, :)
    real(dp), allocatable :: alpha(:), beta(:), evals(:), tridiag_vec(:, :)
    real(dp), allocatable, target :: evecs_tridiag(:, :)
    real(dp), pointer :: real_ptr(:)


    complex(dp), allocatable :: Q_lanczos(:, :), exc_evecs_save(:, :, :), exc_evecs(:, :)
    complex(dp), allocatable, target :: Q_lanczos_gq(:, :)
    complex(dp), pointer :: complex_ptr(:)
    

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
    n_omega  = input%xs%energywindow%points
    bse_type        = trim(adjustl(input%xs%BSE%bsetype))
    nlanczos        = select_nlanczos(input, n_transitions)
    save_exc_evecs  = input%xs%fastBSE%saveQ

    ! Setup omega
    omega = linspace(omega_intervall, n_omega)

    ! Initialize arrays
    allocate(exc_evals_save(2 * nlanczos - 1, 3), source = 0._dp)
    allocate(oscstr_save(2 * nlanczos - 1, 3), source = 0._dp)
    allocate(absspec(n_omega, 3), source = 0._dp)
    if (save_exc_evecs) allocate(exc_evecs_save(n_transitions, 2 * nlanczos - 1, 3), source = zzero)

    do i_dim=1, 3

      ! Run Lanczos iteration
      if(save_exc_evecs) then
        call lanczos(nlanczos, bsh_times_vector, dipole_matrix_elements(:, i_dim), alpha, beta, Q_lanczos)
      else
        call lanczos(nlanczos, bsh_times_vector, dipole_matrix_elements(:, i_dim), alpha, beta)
      end if
      call terminate_if_false(allocated(alpha), 'Error(fastBSE/diagonalize): For i_dim = ' // to_char(i_dim) &
              // 'Lanczos broke down in the first iteration. This means that the compressed BSH is extremly linear ' &
              // 'dependent. Thus, either the ISDF does not work well and you need to increase input%xs%fastBSEnisdf ' &
              // 'or fastBSE is not appropriate for your problem. This might be the case for problems with only few ' &
              // 'transitions.')
      n_its = size(alpha)

      ! Prepare tridiagonal matrix for gauss quadrature
      alpha = [alpha(1 : n_its), alpha(mesh_1d(n_its - 1, 1))]
      beta = [beta(1 : n_its), beta(mesh_1d(n_its - 2, 1))]

      ! Diagonalize
      call diagonalize_symtridiag(alpha, beta, evals, evecs_tridiag)
      
      ! Filter out elements with eigen values < 0
      eval_mask = evals >= 0._dp 
      evals = pack(evals, eval_mask)
      oscstr = pack(evecs_tridiag(1, :) ** 2, eval_mask)
      n_exc = size(evals)

      ! Calculate spectrum
      call calculate_absorption_spectrum(evals, oscstr, omega, broadening, lorentzian, absspec(:, i_dim))
      absspec(:, i_dim) = scaling * norm(dipole_matrix_elements(:, i_dim))**2 * absspec(:, i_dim)
      
      ! Save results
      n_its_save(i_dim) = n_its 
      n_exc_save(i_dim) = n_exc
      exc_evals_save(:n_exc, i_dim) = evals
      oscstr_save(:n_exc, i_dim) = oscstr

      ! Calculate approximation to eigen vectors if wished
      if (save_exc_evecs) then

        evec_mask = reshape(spread(eval_mask, dim=2, ncopies=size(eval_mask)) &
                      .and. spread(eval_mask, dim=1, ncopies=size(eval_mask)), [size(eval_mask)**2])
        real_ptr(1 : size(evecs_tridiag)) => evecs_tridiag(:, :)
        evecs_tridiag = reshape(pack(real_ptr, evec_mask), [n_exc, n_exc])

        allocate(Q_lanczos_gq(n_transitions, 2 * n_its - 1))
        Q_lanczos_gq(:, : n_its) = Q_lanczos
        Q_lanczos_gq(:, n_its + 1: ) = Q_lanczos(:, mesh_1d(n_its - 1, 1))
        Q_lanczos_mask = reshape(spread(eval_mask, dim=1, ncopies=n_transitions), [size(Q_lanczos_gq)])
        complex_ptr(1 : size(Q_lanczos_gq)) => Q_lanczos_gq
        Q_lanczos_gq = reshape(pack(complex_ptr, Q_lanczos_mask), [n_transitions, n_exc])

        allocate(exc_evecs(n_transitions, n_exc))
        call matrix_multiply(Q_lanczos_gq, evecs_tridiag, exc_evecs)
        exc_evecs_save(:, :n_exc, i_dim) = exc_evecs
        deallocate(Q_lanczos, Q_lanczos_mask, Q_lanczos_gq, exc_evecs)
      end if


      deallocate(alpha, beta, evals, oscstr, evecs_tridiag)
      if(save_exc_evecs) deallocate(eval_mask)
      
    end do

    ! Symmetrize quantities
    if(.not. input%xs%BSE%nosymspec) then
      call setup_symmetric_matrix(symmetric_matrix)
      symmetric_matrix = transpose(symmetric_matrix) ! Transpose such that the larger arrays do not need to be transposed
      
      absspec = matmul(absspec, symmetric_matrix)
      exc_evals_save = matmul(exc_evals_save, symmetric_matrix)
      oscstr_save = matmul(oscstr_save, symmetric_matrix)

      if(save_exc_evecs) then
        do i_exc=1, nlanczos
          exc_evecs = exc_evecs_save(:, i_exc, :)
          exc_evecs = matmul(exc_evecs, symmetric_matrix)
          exc_evecs_save(:, i_exc, :) = exc_evecs
        end do
      end if
    end if

    ! Write results to drive
    call h5%initialize(h5file, mpi_env)
    call h5%initialize_group(h5path, result_group)
    group = join_paths(h5path, result_group)

    call h5%initialize_group(group, bse_type)
    group = join_paths(group, bse_type)

    call h5%write(group, ip_bandgap_dataset, bsh%ip_gap())
    call h5%write(group, n_its_dataset, n_its_save)
    call h5%write(group, n_exc_dataset, n_exc_save)
    call h5%write(group, omega_dataset, omega)
    call h5%write(group, absspec_dataset, absspec)
    call h5%write(group, exc_evals_dataset, exc_evals_save)
    call h5%write(group, oscstr_dataset, oscstr_save)
    if (save_exc_evecs) call h5%write(group, exc_evecs_dataset, exc_evecs_save)
    call h5%finalize()

    deallocate(omega, absspec, exc_evals_save, oscstr_save)
    if (save_exc_evecs) deallocate(exc_evecs_save)

    call timesec(t_end)

    write(info_unit, '(A)')        'Diagonalization done.'
    write(info_unit, '(A, I8)')    'Lanczos iterations done: ', n_its 
    write(info_unit, '(A, F15.6)') 'Time (s): ', t_end - t_start
    write(info_unit, *)    

    contains 

    !> Select the number of lanzcos n_its from the input file. There are two parameters that allow to define that number.
    !> One way is to set it directly (`input%xs%fastBSE%nlanczos`), the other way is to define it as a portion of the number of
    !> transitions (`input%xs%fastBSE%clanczos`). 
    !>
    !> If `input%xs%fastBSE%nlanczos > 0` this value is taken. The default is defined as `0`. 
    !> Else `int(input%xs%fastBSE%clanczos * n_transitions)` is taken.
    integer(i32) function select_nlanczos(input, n_transitions)
      type(input_type), intent(in) :: input 
      integer(i32), intent(in) :: n_transitions

      call warn_if_false(input%xs%fastBSE%nlanczos > n_transitions, 'fastBSE: diagonalize: ' &
              // 'input%xs%fastBSE%nlanczos > n_transitions: Parameter will be set to the maximum limit of ' &
              // 'n_transitions for further execution!')

      select_nlanczos = min(input%xs%fastBSE%nlanczos, n_transitions)
    end function


    !> Wrapper for applying the Bethe Salpeter Hamiltonian to a vector that has the
    !> correct signiture as input for [[lanczos]].
    subroutine bsh_times_vector(vector_in, vector_out)
      !> Input vector.
      complex(dp), intent(in) :: vector_in(:)
      !> Output vector.
      complex(dp), intent(out) :: vector_out(:)

      call bsh%times_vector(vector_in, vector_out)
    end subroutine

  end subroutine


  !> Read fastBSE results from HDF5 fike and write formatted output files for the spectrum, 
  !> the eigen energies, and the oscillator strengths.
  subroutine fastBSE_human_readable_output(mpi_env, input, h5file, h5path)
    !> MPI environment.
    type(mpiinfo), intent(inout) :: mpi_env
    !> Input file container.
    type(input_type) :: input
    !> Path to the HDF5 file. It contains input data as wave functions and
    !> is used to store the results of the run.
    character(*), intent(in) :: h5file
    !> Group in the HDF5 file with the input data and to store the results to.
    character(*), intent(in) :: h5path

    type(xhdf5_type) :: h5 
    character(:), allocatable :: group, bse_type

    integer(i32) :: n_exc(3)
    integer(i32), allocatable :: dset_shape(:)
    real(dp) :: ip_gap, energy_conversion
    real(dp), allocatable :: omega(:), absspec(:, :), exc_evals(:, :), oscstr(:, :)

    energy_conversion = 1._dp
    if(input%xs%tevout) energy_conversion = hartree_to_ev

    bse_type = trim(adjustl(input%xs%BSE%bsetype))
    group = join_paths(h5path, result_group)
    group = join_paths(group, bse_type)

    ! Read data
    call h5%initialize(h5file, mpi_env)

    call h5%read(group, ip_bandgap_dataset, ip_gap)
    
    call h5%dataset_shape(group, omega_dataset, dset_shape)
    allocate(omega(dset_shape(1)))
    call h5%read(group, omega_dataset, omega)

    call h5%dataset_shape(group, absspec_dataset, dset_shape)
    allocate(absspec(dset_shape(1), dset_shape(2)))
    call h5%read(group, absspec_dataset, absspec)

    call h5%dataset_shape(group, exc_evals_dataset, dset_shape)
    allocate(exc_evals(dset_shape(1), dset_shape(2)))
    call h5%read(group, exc_evals_dataset, exc_evals)

    call h5%dataset_shape(group, oscstr_dataset, dset_shape)
    allocate(oscstr(dset_shape(1), dset_shape(2)))
    call h5%read(group, oscstr_dataset, oscstr)

    call h5%finalize()


    ! Write Human readable output
    call write_absorption_spectrum_textfile(fname_calculate_absorption_spectrum , omega, absspec, input%xs%broad, energy_conversion)
    call write_exciton_energy_textfile(fname_exciton_energies, exc_evals, energy_conversion, ip_gap)
    call write_oscillator_strength_textfile(fname_oscillator_strengths, oscstr)

    deallocate(dset_shape, omega, absspec, exc_evals, oscstr)

    contains


    !> Write the imaginary part of the diagonal of \(\epsilon^M\) to a text file.
    !> The first column contains \(\omega\), the following \(\epsilon^M_{11}), \(\epsilon^M_{22}) and \(\epsilon^M_{33}).
    subroutine write_absorption_spectrum_textfile(fname, omega, absspec, broadening, energy_conversion)
      !> Name of the file
      character(*), intent(in) :: fname
      !> Energy grid for \(\omega\)
      real(dp), intent(in) :: omega(:)
      !> Imaginary part of the diagonal of \(\epsilon^M\) with dimensions stored columnwise.
      real(dp), intent(in) :: absspec(:, :)
      !> Broadening used to generate the spectrum.
      real(dp), intent(in) :: broadening
      !> Energy scaling. As usual 1.0 means the energy is in Hartree.
      real(dp), intent(in) :: energy_conversion

      integer(i32) :: n_omega, i_omega, unit 

      n_omega = size(omega)

      open(newunit=unit, file=fname, form='formatted', action='write', status='replace')

      write(unit, '(A)') '# fastBSE imaginary macroscopic dielectric function'
      write(unit, '(A)') '# '
      write(unit, '(A, E23.16, A)') '# Energy unit: ', 1 / energy_conversion, ' Hartree'
      write(unit, '(A, E23.16, A)') '# Broadening:  ', broadening * energy_conversion, ' energy unit'
      write(unit, '(A)') '#'
      write(unit, '(A, A22, 1x, A23, 1x, A23, 1x, A23)') '#',  'omega', 'oc11', 'oc22', 'oc33'
      write(unit, '(SP, E23.16, 1x, E23.16, 1x, E23.16, 1x, E23.16)') &
              (omega(i_omega) * energy_conversion, absspec(i_omega, 1), absspec(i_omega, 2), absspec(i_omega, 3), i_omega=1, n_omega)
      
      close(unit)
    end subroutine


    !> Write the exciton eigen energies to a text file.
    !> The three columns contain the eigen energies corresponding to the three lanczos runs
    !> for the directions of \(\langle \mathbf p \rangle\) as starting point.
    subroutine write_exciton_energy_textfile(fname, exc_evals_save, energy_conversion, ip_gap)
      !> Name of the file
      character(*), intent(in) :: fname
      !> Exciton eigen energies from Lanczos runs.
      real(dp), intent(in) :: exc_evals_save(:, :)
      !> Energy scaling. As usual 1.0 means the energy is in Hartree.
      real(dp), intent(in) :: energy_conversion
      !> Independent particle band gap
      real(dp), intent(in) :: ip_gap

      integer(i32) :: n_exciton, i_exciton, unit 

      n_exciton = size(exc_evals_save, dim=1)

      open(newunit=unit, file=fname, form='formatted', action='write', status='replace')

      write(unit, '(A)') '# fastBSE exciton eigen energies'
      write(unit, '(A)') '# The three columns correspond to the results of the three Lanczos runs, each for one of the'
      write(unit, '(A)') '# directions of <p> as starting point.'
      write(unit, '(A)') '# '
      write(unit, '(A, E23.16, A)') '# Energy unit: ', 1 / energy_conversion, ' Hartree'
      write(unit, '(A, E23.16, A)') '# IP band gap: ', ip_gap * energy_conversion, ' energy unit'
      write(unit, '(A)') '# '
      write(unit, '(A, A22, 1x, A23, 1x, A23)') '#', 'E -> <p_1>', 'E -> <p_2>', 'E -> <p_3>'
      write(unit, '(SP, E23.16, 1x, E23.16, 1x, E23.16)') &
              (exc_evals_save(i_exciton, 1) * energy_conversion, exc_evals_save(i_exciton, 2) * energy_conversion, exc_evals_save(i_exciton, 3) * energy_conversion, i_exciton=1, n_exciton)
      
      close(unit)
    end subroutine


    !> Write the exciton eigen energies to a text file.
    !> The three columns contain the eigen energies corresponding to the three lanczos runs
    !> for the directions of \(\langle \mathbf p \rangle\) as input for the Lanczos algorithm..
    subroutine write_oscillator_strength_textfile(fname, oscillator_strength)
      !> Name of the file
      character(*), intent(in) :: fname
      !> Exciton eigen energies from Lanczos runs.
      real(dp), intent(in) :: oscillator_strength(:, :)

      integer(i32) :: n_exciton, i_exciton, unit 

      n_exciton = size(oscillator_strength, dim=1)

      open(newunit=unit, file=fname, form='formatted', action='write', status='replace')

      write(unit, '(A)') '# fastBSE oscillator strength'
      write(unit, '(A)') '# The three columns correspond to the results of the three Lanczos runs, each for one of the'
      write(unit, '(A)') '# directions of <p> as input for the lanczos algorithm.'
      write(unit, '(A)') '# '
      write(unit, '(A, A22, 1x, A23, 1x, A23)') '#', 'osc. str. -> <p_1>', 'osc. str. -> <p_2>', 'osc. str. -> <p_3>'
      write(unit, '(SP, E23.16, 1x, E23.16, 1x, E23.16)') &
              (oscillator_strength(i_exciton, 1), oscillator_strength(i_exciton, 2), oscillator_strength(i_exciton, 3), i_exciton=1, n_exciton)
      close(unit)
    end subroutine

  end subroutine 

  
  !> Read excitons as calculated with `fastBSE`.
  subroutine read_fastBSE_excitons_hdf5(input, h5file, h5path, exc_evals_save, exc_evecs_save)
    !> Input file container.
    type(input_type) :: input
    !> Path to the HDF5 file. It contains input data as wave functions and
    !> is used to store the results of the run.
    type(xhdf5_type), intent(inout) :: h5file
    !> Group in the HDF5 file with the input data and to store the results to.
    character(*), intent(in) :: h5path
    !> Exciton eigen values
    real(dp), allocatable :: exc_evals_save(:, :)
    !> Exciton eigen vectors
    complex(dp), allocatable :: exc_evecs_save(:, :, :)

    character(:), allocatable :: group, bse_type
    integer(i32) :: n_excitons, n_transitions
    integer(i32), allocatable :: dset_shape(:)

    bse_type = trim(adjustl(input%xs%BSE%bsetype))
    group = join_paths(h5path, result_group)
    group = join_paths(group, bse_type)

    call h5file%dataset_shape(group, exc_evecs_dataset, dset_shape, complex_dataset=.true.)

    n_transitions = dset_shape(1)
    n_excitons = dset_shape(2)

    if (allocated(exc_evals_save)) deallocate(exc_evals_save)
    allocate(exc_evals_save(n_excitons, 3))
    call h5file%read(group, exc_evals_dataset, exc_evals_save, [1, 1])

    if (allocated(exc_evecs_save)) deallocate(exc_evecs_save)
    allocate(exc_evecs_save(n_transitions, n_excitons, 3))
    call h5file%read(group, exc_evecs_dataset, exc_evecs_save, [1, 1, 1])
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
      call terminate_mpi_env(mpi_env, &
              'Error(fastBSE): fastBSE only supports TDA.')
    end if 

    if(input%xs%BSE%xas) then
      call terminate_mpi_env(mpi_env, &
              'Error(fastBSE): fastBSE only supports valence excitations.')
    end if

    if(input%xs%fastBSE%nlanczos <= 0) then
      call terminate_mpi_env(mpi_env, &
              'Error(fastBSE): input%xs%fastBSE%nlanczos < 0. Choose a value > 0.')
    end if

  end subroutine

end module fastBSE   
