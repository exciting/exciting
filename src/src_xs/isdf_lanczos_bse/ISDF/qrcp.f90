!> Routines for calculating periodic interpolative seperable density fitting (ISDF). 
!> ISDF approximates wavefunction products on a interpolation grid.
!> For periodic systems, ISDF decomposition has the following general form:
!> \[ 
!>   Z_{i\boldsymbol{k},j\boldsymbol{k}'}(\boldsymbol{r}) := 
!>   u_{i\boldsymbol{k}}(\boldsymbol{r}) u_{\boldsymbol{k}'}(\boldsymbol{r}) 
!> \]
!> \[ 
!>   \approx \sum_{\mu=1}^{N_\mu} \zeta_\mu(\boldsymbol{r}) 
!>   u_{i\boldsymbol{k}}(\boldsymbol{r}_\mu) u_{\boldsymbol{k}'}(\boldsymbol{r}_\mu) 
!> \]
!> where (\ u_{i\boldsymbol{k}}(\boldsymbol{r}) \) is the wavefunction for state and k-point 
!> \( i, \boldsymbol{k} \) at the real space grid point \( \boldsymbol{r} \),
!> \( zeta_\mu(\boldsymbol{r}) /) is the interpolation vector for interpolation point 
!> \( \boldsymbol{r}_\mu \) and at the real space grid point \( \boldsymbol{r} \).
!> Seen \( \{Z_{i \boldsymbol{k}, j \boldsymbol{k}'}(\boldsymbol{r})\} \) as a \( N_r \times N_kN_iN_j \) 
!> matrix, choosing sufficiant interpolation points can be undestood as selecting  rows such 
!> that the other rows can be represented by a linear combination of these rows. 
!> This can be achieved by calculating a QR factorization with column pivoting of \( Z \).
!> Calculating the interpolation vactors can be achieved by a leat squares approximation.
!> Included routines for calculating interpolation points (qrcp) and interpolation vectors (isdf).
module qr_randomized_column_pivoting
  use precision, only: dp
  use constants, only: pi, zi, zzero
  use modmpi, only: terminate_if_false
  use xlapack, only: qr_column_pivot, outer_product
  use math_utils, only: shuffle_vector
  use multi_index_conversion, only: ind, composit_index_to_indices
  use seed_generation, only: set_seed
  use fftw_wrapper, only: fftw_plan_type, setup_fftw_plan, execute_fft

  implicit none

  private
  public :: qrcp

  !> Default seed source for random number generation
  character(1) :: seed_source_default = 'N'

  !> Calculate the indices of the interpolation points for the ISDF compression of the product between 
  !> wavefunctions.
  interface qrcp
    module procedure qrcp_ii, qrcp_ij
  end interface qrcp
    
contains
    
  !> Calculate the indices of the interpolation points for the ISDF compression of the product between 
  !> wavefunctions from a single set, 'e.g.' wavefunctions of the occupied states form such a set.
  !> The indices are calculated with randomized QR factorization with culomn pivoting (QRCP) as proposed in
  !> Lu, Ying; 2016; Fast Allgorithm for periodic Fitting for Bloch Waves.
  function qrcp_ii(N_mu, N_sub, u, seed_source_input) result(r_mu)
    !> Number of interpolation points
    integer, intent(in) :: N_mu
    !> Number of subsampling points
    integer, intent(in) :: N_sub
    !> Wave function set
    complex(dp), intent(in) :: u(:, :, :)
    !> Source for seed for random number generation.
    !> See [[set_seed]]
    character(1), intent(in), optional :: seed_source_input

    integer :: r_mu(N_mu)

    integer :: N_r, N_k, N, ir, isub1, isub2
    character(1) :: seed_source
    integer, allocatable :: rand_order(:), r_work(:)
    real(dp), allocatable :: random_phase(:)
    complex(dp), allocatable :: random_unit(:), u_work(:), M(:,:), M_r(:, :)
    type(fftw_plan_type) :: fft_plan

    N_r = size(u, dim=1)
    N   = size(u, dim=2)
    N_k = size(u, dim=3)

    ! Setup random complex numbers with lenght 1.0
    seed_source = seed_source_default
    if (present(seed_source_input)) seed_source = seed_source_input
    call set_seed(seed_source_input)

    allocate(random_phase(N * N_k))
    call random_number(random_phase)
    random_unit = exp(2 * zi * pi * random_phase)

    rand_order = shuffle_vector(N * N_k, seed_source)

    ! Setup fftw plan
    fft_plan = setup_fftw_plan([N * N_k], sign=-1, flag_effort="FFTW_EXHAUSTIVE")

    ! Build up matrix for qr factorization with column pivoting
    allocate(u_work(N * N_k))
    allocate(M_r(N_sub, N_sub), source = zzero)
    allocate(M(N_sub**2, N_r))
    do ir=1, N_r
      ! Reshape Wavefunction and multiply random phase
      u_work = random_unit * reshape(u(ir, :, :), [N * N_k])

      ! Calculate fft
      call execute_fft(fft_plan, u_work)

      ! Reorder the FT of reshaped u
      u_work = u_work(rand_order(1 : N_sub))
      
      ! Setup column of wave function product.
      call outer_product(conjg(u_work), u_work, M_r)

      M(:, ir) = reshape(M_r, [N_sub**2])

    end do

    allocate(r_work(N_r))
    call qr_column_pivot(M, r_work)

    ! choose first N_mu numbers from column pivoting permutation
    r_mu = r_work(1 : N_mu)
  end function qrcp_ii
    

  !> Calculate the indices of the interpolation point for the ISDF compression of the product between wave functions
  !> from two sets of wave functions, 'e.g.' wave functions of the occupied and unoccupied bands.
  !> The indices are calculated randomized QR factorization with culomn pivoting (QRCP) as proposed in
  !> Lu, Ying; 2016; Fast Allgorithm for periodic Fitting for Bloch Waves.
  function qrcp_ij(N_mu, N_sub, u_i, u_j, seed_source_input) result(r_mu)
    !> Number of interpolation points
    integer, intent(in) :: N_mu
    !> Number of subsampling points
    integer, intent(in) :: N_sub
    !> Wave function sets
    complex(dp), intent(in) :: u_i(:,:,:), u_j(:,:,:)
    !> Source for seed for random number generation.
    !> See [[set_seed]]
    character(1), intent(in), optional :: seed_source_input

    integer :: r_mu(N_mu)
       
    integer :: N_r, N_k, N_i, N_j, ir, i, j
    integer, allocatable :: rand_order(:), r_work(:)
    real(dp), allocatable :: random_phase(:)
    complex(dp), allocatable :: random_unit(:), u_ij_reshaped(:, :), u_ij_work(:), M(:,:)
    type(fftw_plan_type) :: fft_plan 
    character(1) :: seed_source
        
    call terminate_if_false(size(u_i, dim=1) == size(u_j, dim=1), &
             "Error(qrcp_ij): The number of real space grid points for u_i and u_j must be the same.")

    call terminate_if_false(size(u_i, dim=3) == size(u_j, dim=3), &
             "Error(qrcp_ij): The number of k-points for u_i and u_j must be the same.")
        
    N_r = size(u_i, dim=1)
    N_k = size(u_i, dim=3)
    N_i = size(u_i, dim=2)
    N_j = size(u_j, dim=2)
        
    allocate(random_phase(N_i * N_j * N_k))
    allocate(u_ij_reshaped(N_i * N_j, N_k))
    allocate(u_ij_work(N_i * N_j * N_k))
    allocate(M(N_sub, N_r))
        
    ! Initialize randomization
    seed_source = seed_source_default
    if (present(seed_source_input)) seed_source = seed_source_input
    call set_seed(seed_source)

    ! Calculate random unit
    call random_number(random_phase)
    random_unit = exp(2 * zi * pi * random_phase)

    ! Setup random permutation
    rand_order = shuffle_vector(N_i * N_j * N_k, seed_source)

    ! Calculate FFT plan
    u_ij_work = zzero
    fft_plan = setup_fftw_plan([N_i * N_j * N_k], u_ij_work, sign=-1, flag_effort="FFTW_EXHAUSTIVE")

    ! Build up matrix for qr factorization with column pivoting
    do ir=1, N_r

      do i=1, N_i
        do j=1, N_j
          u_ij_reshaped(ind([i, j], [N_i, N_j]), :) = conjg(u_i(ir, i, :)) * u_j(ir, j, :)
        end do
      end do

      u_ij_work = random_unit * reshape(u_ij_reshaped, [N_i * N_j * N_k])
          
      ! calculate execute_fft
      call execute_fft(fft_plan, u_ij_work)

      M(:, ir) = u_ij_work(rand_order(1 : N_sub)) ! step 3 & 4
    end do

    allocate(r_work(N_r))
    call qr_column_pivot(M, r_work)
    !write(*, *) shape(M)

    r_mu = r_work(1 : N_mu)      
  end function qrcp_ij
end module qr_randomized_column_pivoting