! Created by  on 21/02/2022.

module wscr_isdf_kernel
  use precision, only: dp
  use constants, only: zone, zzero
  use asserts, only: assert

  use math_utils, only: mod1
  use grid_utils, only: n_grid_diff, flattened_map
  use xlapack, only: matrix_multiply
  use lapack_f95_interfaces, only: zgemm, zhemm
  use xfftw, only: fft_type, FFTW_FORWARD, FFTW_BACKWARD
  use multi_index_conversion, only: indices_to_composite_index, composite_index_to_indices

  use regular_grid, only: regular_grid_type
  use dynamic_indices, only: dynindex_type

  use interaction_kernel, only: interaction_kernel_type

  use bse_transitions, only: transition_type

  use pointer_group, only: matrix_pointer_group_cmplx_dp

  use modmpi, only: mpiinfo


  private
  public :: wscr_isdf_kernel_type

  type, extends(interaction_kernel_type) :: wscr_isdf_kernel_type
    !> Sizes
    integer, private :: n_occupied_bands, n_q, n_isdf_o, n_isdf_u, q_sampling(3), k_sampling(3)
    !> Fourier transformation over q-space of fastBSE screened interaction kernel
    complex(dp), allocatable :: W_tilde_hat(:, :, :)
    !> Periodic part of real space wave functions on the interpolation grid
    complex(dp), private, allocatable :: u_o_mu(:, :), u_u_mu(:, :)
    !> Interger map to zero padded extension
    integer, allocatable :: map_to_zero_padded(:)
    !> FFTW plans for q grid FFT
    type(fft_type) :: q_fft_fw, q_fft_bw
    !> Workspace for matrix vector multiplication
    complex(dp), allocatable :: A(:, :), B(:, :, :), D(:)

    type(matrix_pointer_group_cmplx_dp), allocatable, private :: A_(:), u_o_mu_(:), u_u_mu_(:), D_(:)

    contains
                    
    procedure :: intialize
    procedure :: finalize
    procedure :: times_vector
  end type wscr_isdf_kernel_type

  contains

  subroutine intialize(this, mpi_env, transitions, u_o_mu, u_u_mu, zeta_o_hat, zeta_u_hat, scaling, w_hat, Gq_indices, g_grid, k_sampling)
    class(wscr_isdf_kernel_type) :: this
    type(mpiinfo), intent(inout) :: mpi_env
    type(transition_type), intent(in) :: transitions
    !> Periodic part of the occupied/unoccupied states wavefunctions.
    complex(dp), intent(in) :: u_o_mu(:, :),  u_u_mu(:, :)
    !> Fourier transform of the interpolation vectors as obtained by [[isdf]] for the left wavefunction set.
    complex(dp), intent(in), contiguous :: zeta_u_hat(:, :), zeta_o_hat(:, :)
    !> Scaling of the screened Coulomb potential
    real(dp), intent(in) :: scaling
    !> Fourier transtransform of the screened Coulomb potential
    complex(dp), intent(in) :: w_hat(:, :, :)
    !>  \(\mathbf{G} + \mathbf{k}\)-grid on which the screened interaction is defined
    type(dynindex_type) :: Gq_indices
    !> \(\mathbf{G} + \mathbf{k}\)-grid on which the screened interaction is defined
    type(regular_grid_type), intent(in) :: g_grid
    !> Number of \(\mathbf{k}\)-points per dimension
    integer, intent(in) :: k_sampling(3)

    integer :: i, imu, inu
    type(fft_type) :: fft
    complex(dp), allocatable :: fft_space(:)

    this%k_sampling = k_sampling
    this%n_k = product(k_sampling)
    this%q_sampling = n_grid_diff(k_sampling)
    this%n_q = product(this%q_sampling)
    
    this%map_to_zero_padded = flattened_map(this%q_sampling, this%k_sampling)

    this%n_isdf_u = size(u_u_mu, dim=1)
    this%n_isdf_o = size(u_o_mu, dim=1)

    this%transitions = transitions
    
    this%u_o_mu = u_o_mu
    this%u_u_mu = u_u_mu

    ! Calculate fastBSE kernel
    this%W_tilde_hat = scaling * W_curl(w_hat, k_sampling, g_grid, Gq_indices, zeta_o_hat, zeta_u_hat)
    call fft%initialize(this%q_sampling, FFTW_FORWARD)
    do inu=1, this%n_isdf_o
      do imu=1, this%n_isdf_u
        call fft%execute(this%W_tilde_hat(:, imu, inu))
      end do
    end do
    call fft%delete()

    call initialize_workspace(this)

  end subroutine intialize

  subroutine initialize_workspace(this)
    type(wscr_isdf_kernel_type), target, intent(inout) :: this

    integer :: ik 

    ! Setup workspace
    allocate(this%A(this%n_isdf_o, this%transitions%n_u))
    allocate(this%B(this%n_isdf_u, this%n_isdf_o, this%n_k))
    allocate(this%D(this%transitions%n_allowed))

    allocate(this%u_o_mu_(this%transitions%n_k))
    allocate(this%u_u_mu_(this%transitions%n_k))
    allocate(this%A_(this%transitions%n_k))
    allocate(this%D_(this%transitions%n_k))

    do ik=1, this%transitions%n_k 
      
      this%u_o_mu_(ik)%p(1 : this%n_isdf_o, 1 : this%transitions%n_o_per_k(ik)) &
              => this%u_o_mu(:, this%transitions%o_first(ik) : this%transitions%o_last(ik))

      this%u_u_mu_(ik)%p(1 : this%n_isdf_u, 1 : this%transitions%n_u_per_k(ik)) &
              => this%u_u_mu(:, this%transitions%u_first(ik) : this%transitions%u_last(ik))

      this%A_(ik)%p(1 : this%n_isdf_o, 1 : this%transitions%n_u_per_k(ik)) &
              => this%A(:, this%transitions%u_first(ik) : this%transitions%u_last(ik))

      this%D_(ik)%p(1 : this%transitions%n_u_per_k(ik), 1 : this%transitions%n_o_per_k(ik)) &
              => this%D(this%transitions%first(ik) : this%transitions%last(ik))
    
    end do 


    call this%q_fft_fw%initialize(this%q_sampling, FFTW_FORWARD)
    call this%q_fft_bw%initialize(this%q_sampling, FFTW_BACKWARD)

  end subroutine initialize_workspace

  subroutine finalize(this)
    class(wscr_isdf_kernel_type), intent(inout) :: this

    integer :: ik

    if(allocated(this%u_o_mu)) deallocate(this%u_o_mu)
    if(allocated(this%u_u_mu)) deallocate(this%u_u_mu)
    if(allocated(this%W_tilde_hat)) deallocate(this%W_tilde_hat)
    if(allocated(this%map_to_zero_padded)) deallocate(this%map_to_zero_padded)
    if(allocated(this%A)) deallocate(this%A)
    if(allocated(this%B)) deallocate(this%B)
    if(allocated(this%D)) deallocate(this%D)

    do ik=1, this%transitions%n_k
      call this%u_o_mu_(ik)%finalize()
      call this%u_u_mu_(ik)%finalize()
      call this%A_(ik)%finalize()
      call this%D_(ik)%finalize()
    end do

  end subroutine finalize


  !> Action of the resonant-resonant screened interaction kernel part \( W_A \) on a vector \( X \)
  !> \[
  !>   [W_A\cdot X](\mathbf k i_o i_u) = \sum_{\mathbf{k'}, j_o, j_u} W_A(\mathbf{k} i_o i_u, \mathbf{k'} j_o j_u)
  !>       X(\mathbf{k'} j_o j_u)
  !> \]
  !> where \( i_u, i_o, j_o \) and \( j_u \) are the indices of the occupied and unoccupied bands and \( \mathbf{k} \)
  !> and \( \mathbf{k'} \) the k-points for which the exchange interaction kernel is obtained. With
  !> \[
  !>   W_A(\mathbf{k} i_o i_u, \mathbf{k'} j_o j_u) = \frac{1}{N_k} \sum_{\mu,\nu=1}^{N_\mu^W}
  !>       \bar{u}_{\mathbf{k}i_u}(\mathbf{\hat{r}}_\mu^W) u_{\mathbf{k'}j_u}(\mathbf{\hat{r}}_\mu^W)
  !>           \tilde{W}_{\mathbf{k-k'}, \mu\nu} \bar{u}_{\mathbf{k'}j_o}(\mathbf{\hat{r}}_\nu^W)
  !>               u_{\mathbf{k}i_o}(\mathbf{\hat{r}}_\nu^W)
  !> \]
  !> where \( u_{\mathbf{k} i} ( \mathbf{\hat{r}}_\mu) \) is the wave function of the band `i` at k-point
  !> \( \mathbf{k} \) evaluated at the ISDF interpolation point \( \mathbf{\hat{r}}_\mu \), \( N_\mu^W \) is the
  !> number of ISDF interpolation points and \( \tilde{W}_{A,\mu\nu} = \mathcal_W(\zeta^{uu}_\mu, \zeta^{oo}_\nu) \)
  !> is the screened interaction kernel with ISDF compression as obtained by [assemble_W_curl]. Putting these
  !> equations together leads to
  !> \[
  !>   [W_A\cdot X](\mathbf k i_o i_u) = \frac{1}{N_k}\sum_{} \sum_{\mathbf{k'}} \sum_{j_o=1}^{N_v} \sum_{j_u=1}^{N_c}
  !>       \sum_{\mu=1}^{N_\mu^W} \sum_{\nu=1}^{N_\mu^W} \bar{u}_{\mathbf{k}i_u}(\mathbf{\hat{r}}_\mu^W)
  !>           u_{\mathbf{k'}j_u}(\mathbf{\hat{r}}_\mu^W) \tilde{W}_{\mathbf{k-k'}, \mu\nu}
  !>               \bar{u}_{\mathbf{k'}j_o}(\mathbf{\hat{r}}_\nu^W) u_{\mathbf{k}i_o}(\mathbf{\hat{r}}_\nu^W)
  !>                   X(\mathbf{k'} j_o j_u)
  !> \]
  !> where \( n_occupied_bands, n_u \) are the number of occupied and unoccupied bands considered for the solution of
  !> the BSE respectively. Regrouping the sum reveals the efficient action on the vector and leads to the
  !> expression that is evaluated in this routine:
  !> \[
  !>   \frac{1} {N_k} \sum_{\nu=1}^{N_\mu^W} u_{\mathbf{k} i_o}(\mathbf{\hat{r}}_\nu^W) \left\{\sum_{\mu=1}^{N_\mu^W}
  !>       \bar{u}_{\mathbf{k}i_u}(\mathbf{\hat{r}}_\mu^W) \left[\sum_{\mathbf{k'}} \tilde{W}_{\mathbf{k-k'}, \mu\nu}
  !>           \left(\sum_{j_u=1}^{N_c} u_{\mathbf{k'}j_u}(\mathbf{\hat{r}}_\mu^W) \left(\sum_{j_o=1}^{N_v}
  !>               \bar{u}_{\mathbf{k'}j_o}(\mathbf{\hat{r}}_\nu^W) X(\mathbf{k'} j_o j_u)\right)\right)\right]
  !>                   \right\}
  !> \]
  subroutine times_vector(this, vector)
    class(wscr_isdf_kernel_type), target, intent(inout) :: this
    complex(dp), contiguous, target, intent(inout) :: vector(:)

    integer :: ik, imu, inu
    complex(dp), allocatable :: fft_space(:), C(:, :)

    allocate(fft_space(this%n_q))

    this%D = vector 
    
    do ik=1, this%n_k
      call matrix_multiply(conjg(this%u_o_mu_(ik)%p), this%D_(ik)%p, this%A_(ik)%p, trans_B='T')
      call matrix_multiply(this%u_u_mu_(ik)%p, this%A_(ik)%p, this%B(:, :, ik), trans_B='T')
    end do

    !$omp parallel do default(shared) private(imu, inu, fft_space) collapse(2)
    do inu=1, this%n_isdf_o
      do imu=1, this%n_isdf_u
        ! convolution via fast fft over the k grid
        fft_space = zzero
        fft_space(this%map_to_zero_padded) = this%B(imu, inu, :)
        call this%q_fft_fw%execute(fft_space, rescale_forward=.false.)
        fft_space(:) = this%W_tilde_hat(:, imu, inu) * fft_space(:)
        call this%q_fft_bw%execute(fft_space)
        this%B(imu, inu, :) = fft_space(this%map_to_zero_padded)
      end do 
    end do
    !$omp end parallel do

    do ik=1, this%n_k
      allocate(C(this%transitions%n_u_per_k(ik), this%n_isdf_o))
      call matrix_multiply(this%u_u_mu_(ik)%p, this%B(:, :, ik), C, trans_A = 'C')
      call matrix_multiply(C, this%u_o_mu_(ik)%p, this%D_(ik)%p)
      deallocate(C)
    end do

    vector = this%D

  end subroutine times_vector

  !> Calculate the screened interaction kernel with ISDF:
  !> \[
  !>  \mathcal{W}_\vec{q} (\zeta_{l,mu}, \zeta_{r,\nu}) = |\Omega| \sum_{\vec{G},\vec{G}'}
  !>          \bar{\hat{\zeta}}_{l,\mu}(\vec{G}) \hat{W}_\vec{q}(\vec{G},\vec{G}') \hat{\zeta}_{r,\nu}(\vec{G'}).
  !> \]
  function W_curl(w_hat, k_sampling, g_grid, Gq_indices, zeta_o_hat, zeta_u_hat) result(W_curl_out)
    !> Screened interaction \( \hat{W} \)
    complex(dp), intent(in) :: w_hat(:, :, :)
    !> \(\mathbf{q}\)-grid on which the screened interaction is defined
    integer, intent(in) :: k_sampling(3)
    !> \(\mathbf{G} + \mathbf{k}\)-grid on which the screened interaction is defined
    type(regular_grid_type), intent(in) :: g_grid
    !> Indices of the \(\mathbf{G}\)-vectors per \(\mathbf{q}\)-vector on the FFT grid.
    type(dynindex_type) :: Gq_indices
    !> ISDF interpolation coefficients \( \zeta_{o/u, mu} \)
    complex(dp), intent(in) :: zeta_o_hat(:, :), zeta_u_hat(:, :)

    complex(dp), allocatable :: W_curl_out(:, :, :)

    integer :: ik, iq, n_q, N_G, n_isdf_o, n_isdf_u, q_sampling(3), g_shift(3), q_int_coord(3), k_int_coord(3)
    integer, allocatable :: G_map(:), G_vecs(:, :)

    complex(dp), allocatable :: w_hat_times_zeta_o(:, :)

    q_sampling = n_grid_diff(k_sampling)
    n_q = product(q_sampling)

    n_isdf_o = size(zeta_o_hat, dim=2)
    n_isdf_u = size(zeta_u_hat, dim=2)

    allocate(W_curl_out(n_q, n_isdf_u, n_isdf_o))
    
    do iq = 1, n_q
      ! Calculate the mapping from G+q to the G+k grid.
      call composite_index_to_indices(iq, q_sampling, q_int_coord)
      k_int_coord = mod1(q_int_coord, k_sampling)
      ik = indices_to_composite_index(k_int_coord, k_sampling)
      g_shift = (q_int_coord - k_int_coord) / k_sampling

      N_G = Gq_indices%n_per_col(ik)
      G_vecs = g_grid%multi_index(Gq_indices%col(ik)) + spread(g_shift, 2, N_G)
      G_map = g_grid%composite_index(G_vecs)

      ! zeta_u_hat x w_hat x zeta_o_hat
      allocate(w_hat_times_zeta_o(N_G, n_isdf_o))
      call zhemm('l', 'u', N_G, n_isdf_o, zone, &
              w_hat(1 : N_G, 1 : N_G, ik), N_G, &
              zeta_o_hat(G_map, :), N_G, zzero, &
              w_hat_times_zeta_o, N_G)

      call zgemm('c', 'n', n_isdf_u, n_isdf_o, N_G, zone, &
              zeta_u_hat(G_map, :), N_G, &
              w_hat_times_zeta_o, N_G, zzero, &
              W_curl_out(iq, :, :), n_isdf_u)
      
      deallocate(w_hat_times_zeta_o)
    end do

    deallocate(G_map, G_vecs)
  end function W_curl

end module wscr_isdf_kernel