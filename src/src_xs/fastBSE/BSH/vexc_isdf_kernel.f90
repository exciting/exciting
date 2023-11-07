! Created by  on 18.02.22.
!> Module for setting up the exchange interaction kernel with ISDF compression (see [[isdf]]) and applying it
!> effitiently to the optical transition vector.
module vexc_isdf_kernel
  use precision, only: dp
  use constants, only: pi, zzero
  use grid_utils, only: mesh_1d
  use xlapack, only: matrix_multiply
  use regular_grid, only: regular_grid_type
  use interaction_kernel, only: interaction_kernel_type
  use modmpi, only: terminate_if_false
  use bse_transitions, only: transition_type
  use pointer_group, only: matrix_pointer_group_cmplx_dp

  implicit none

  private
  public :: vexc_isdf_kernel_type

  !> Type for ISDF exchange interaction
  type, extends(interaction_kernel_type) :: vexc_isdf_kernel_type
    !> Numbers of states, \(\mathbf{k}\)-points and interpolation points
    integer, private :: n_isdf
    !> fastBSE interaction kernel
    complex(dp), allocatable, private :: V_tilde(:, :)
    !> Periodic part of real space wave functions on the interpolation grid
    complex(dp), allocatable, private :: u_o_mu(:, :), u_u_mu(:, :)
    !> Workspace for matrix vector multiplication
    complex(dp), private, allocatable :: A(:, :), B(:), C(:), D(:)

    type(matrix_pointer_group_cmplx_dp), allocatable, private ::  A_(:), D_(:), u_o_mu_(:), u_u_mu_(:)

    contains

    procedure :: initialize 
    procedure :: finalize
    procedure :: times_vector
  end type vexc_isdf_kernel_type

  contains

  !> Constructor
  subroutine initialize(this, transitions, u_o_mu, u_u_mu, zeta_hat, v_coulomb_hat)
    !> fastBSE exchange interaction kernel type
    class(vexc_isdf_kernel_type) :: this
    !>
    type(transition_type), intent(in) :: transitions
    !> Periodic part of the wavefunctions for the occupied states, evaluated on the interpolation grid.
    complex(dp), intent(in), contiguous, target :: u_o_mu(:, :) ! n_isdf x n_o x N_k
    !> Periodic part of the wavefunctions for the unoccupied states, evaluated on the interpolation grid.
    complex(dp), intent(in), contiguous, target :: u_u_mu(:, :) ! n_isdf x n_u x N_k
    !> Fourier transfrom of the interpolation vectors as obtained by [[ISDF]].
    complex(dp), intent(in), contiguous, target :: zeta_hat(:, :) ! N_G x n_isdf
    !> Fourier transform of the bare Coulomb potential
    real(dp), intent(in), contiguous :: v_coulomb_hat(:) ! N_G
    
    this%transitions = transitions
    this%u_o_mu = u_o_mu
    this%u_u_mu = u_u_mu
    this%n_isdf = size(zeta_hat, dim=2)
    this%V_tilde = V_tilde_A(v_coulomb_hat, zeta_hat)
    
    call intialize_workspace(this)
  end subroutine initialize

  subroutine intialize_workspace(this)
    !> fastBSE exchange interaction kernel type
    type(vexc_isdf_kernel_type), target, intent(inout) :: this
    
    integer :: ik

    allocate(this%A(this%n_isdf, this%transitions%n_u))
    allocate(this%B(this%n_isdf))
    allocate(this%C(this%n_isdf))
    allocate(this%D(this%transitions%n_allowed), source=zzero)

    allocate(this%u_o_mu_(this%transitions%n_k))
    allocate(this%u_u_mu_(this%transitions%n_k))
    allocate(this%A_(this%transitions%n_k))
    allocate(this%D_(this%transitions%n_k))

    do ik=1, this%transitions%n_k

      this%u_o_mu_(ik)%p(1 : this%n_isdf, 1 : this%transitions%n_o_per_k(ik)) &
              => this%u_o_mu(:, this%transitions%o_first(ik) : this%transitions%o_last(ik))

      this%u_u_mu_(ik)%p(1 : this%n_isdf, 1 : this%transitions%n_u_per_k(ik)) &
              => this%u_u_mu(:, this%transitions%u_first(ik) : this%transitions%u_last(ik))

      this%A_(ik)%p(1 : this%n_isdf, 1 : this%transitions%n_u_per_k(ik)) &
              => this%A(:, this%transitions%u_first(ik) : this%transitions%u_last(ik))

      this%D_(ik)%p(1 : this%transitions%n_u_per_k(ik), 1 : this%transitions%n_o_per_k(ik)) &
              => this%D(this%transitions%first(ik) : this%transitions%last(ik))

    end do  
  end subroutine


  !> Set up resonant-resonant block of exchange interaction kernel with ISDF compression \(\tilde{W}_{\text{ex}}_{A}\):
  !> \[
  !>   \tilde{W}_{\text{ex}}_{A,\mu\nu} = \frac{|\Omega|} {N_k} \sum_{\mathbf{G} \neq 0}
  !>       \bar{\hat{\zeta}}^{ou}_\mu(\mathbf{G}) \hat{V}(\mathbf{G}) \hat{\zeta}^{ou}_\nu(\mathbf{G}).
  !> \]
  function V_tilde_A(v_hat, zeta_hat) result(W_tilde_A_out)
    !> FT of the bare coulomb interaction in
    real(dp), intent(in), contiguous :: v_hat(:)
    !> FT of ISDF interpolation vectors \( \zeta \) along the first dimension
    complex(dp), intent(in), contiguous :: zeta_hat(:, :)

    complex(dp), allocatable :: W_tilde_A_out(:, :)

    integer :: n_isdf, N_G
    complex(dp), allocatable :: V_times_zeta_hat(:, :)

    N_G = size(zeta_hat, dim=1)
    n_isdf = size(zeta_hat, dim=2)

    V_times_zeta_hat = spread(v_hat, 2, n_isdf) * zeta_hat
    allocate(W_tilde_A_out(n_isdf, n_isdf))
    call matrix_multiply(zeta_hat, V_times_zeta_hat, W_tilde_A_out, trans_A='C')
  end function V_tilde_A


  !> Set up resonant-antiresonant block of exchange interaction kernel with ISDF compression
  !> \( \tilde{W}_{\text{ex}}_{A} \):
  !> \[
  !>   \tilde{W}_{\text{ex}}_{B,\mu\nu} = \frac{|\Omega|} {N_k} \sum_{\mathbf{G} \neq 0}
  !>       \bar{\hat{\zeta}}^{ou}_\mu(\mathbf{G}) \hat{V}(\mathbf{G}) \bar{\hat{\zeta}}^{ou}_\nu(\mathbf{G}).
  !> \]
  function V_tilde_B(v_hat, zeta_hat) result(W_tilde_B_out)
    !> FT of the bare coulomb interaction in
    real(dp), intent(in) :: v_hat(:)
    !> FT of ISDF interpolation vectors \( \zeta \)
    complex(dp), intent(in) :: zeta_hat(:, :)

    complex(dp), allocatable :: W_tilde_B_out(:, :)

    integer :: n_isdf, N_G
    complex(dp), allocatable :: V_times_zeta_hat(:, :)

    N_G = size(zeta_hat, dim=1)
    n_isdf = size(zeta_hat, dim=2)

    V_times_zeta_hat = spread(v_hat, 2, n_isdf) * conjg(zeta_hat)

    allocate(W_tilde_B_out(n_isdf, n_isdf))
    call matrix_multiply(zeta_hat, zeta_hat, W_tilde_B_out, trans_A='C')
  end function V_tilde_B


  !> Destructor
  subroutine finalize(this)
    class(vexc_isdf_kernel_type), intent(inout) :: this

    integer :: ik

    if(allocated(this%u_o_mu)) deallocate(this%u_o_mu)
    if(allocated(this%u_u_mu)) deallocate(this%u_u_mu)
    if(allocated(this%V_tilde)) deallocate(this%V_tilde)
    if(allocated(this%A)) deallocate(this%A)
    if(allocated(this%B)) deallocate(this%B)
    if(allocated(this%C)) deallocate(this%C)
    if(allocated(this%D)) deallocate(this%D)

    do ik=1, this%transitions%n_k
      call this%u_o_mu_(ik)%finalize()
      call this%u_u_mu_(ik)%finalize()
      call this%A_(ik)%finalize()
      call this%D_(ik)%finalize()
    end do

    call this%transitions%finalize()
  end subroutine finalize

  !> Efficient action of exchange interaction kernel \( V_{A/B} \) on a vector \( X \) with ISDF compression
  !> \[
  !>   \left[ V_{A/B} \cdot X \right](\mathbf{k} i_o i_u) = \sum_{\mathbf{k'},j_o,j_u}
  !>       V_{A/B}(\mathbf{k} i_o i_u, \mathbf{k'} j_o j_u)
  !>           X(\mathbf{k'} j_o j_u),
  !> \]
  !> where \( i_u, i_o, j_o \) and \( j_u \) are the indices of the occupied and unoccupied bands and \( \mathbf{k} \)
  !> and \( \mathbf{k'} \) the k-points for which the exchange interaction kernel is obtained. With
  !> \[
  !>   V_{A/B} (\mathbf{k} i_o i_u, \mathbf{k'} j_o j_u) = \frac{1} {N_k}\sum_{\mu,\nu=1}^{N_\mu^V}
  !>       \bar{u}_{\mathbf{k} i_u} (\mathbf{\hat{r}}_\mu)
  !>           u_{\mathbf{k} i_o} (\mathbf{\hat{r}}_\mu)
  !>               \tilde{W}_{\text{ex}}_{A/B,\mu\nu}
  !>                   \bar{u}_{\mathbf{k'} j_o}(\mathbf{\hat{r}}_\nu)
  !>                       u_{\mathbf{k'} j_u}(\mathbf{\hat{r}}_\nu),
  !> \]
  !> where \( u_{\mathbf{k} i} ( \mathbf{\hat{r}}_\mu) \) is the wave function of the band `i` at k-point
  !> \( \mathbf{k} \) evaluated at the ISDF interpolation point \( \mathbf{\hat{r}}_\mu \), \( N_\mu^V \) is the number
  !> of ISDF interpolation points and \( \tilde{W}_{\text{ex}}_{A/B,\mu\nu} \) is the exchange interaction kernel with
  !> ISDF compression (see [assemble_VA_tilde] and [assemble_VB_tilde]). Putting these equations together leads to
  !> \[
  !>   [V_{A/B}\cdot X](\mathbf{k} i_o i_u)=
  !>       \frac{1} {N_k}
  !>           \sum_{\mathbf{k'}} \sum_{j_o=1}^{n_o} \sum_{j_u=1}^{n_u} \sum_{\mu=1}^{N_\mu^V} \sum_{\nu=1}^{N_\mu^V}
  !>               \bar{u}_{\mathbf{k} i_u}(\mathbf{\hat{r}}_\mu)
  !>                   u_{\mathbf{k} i_o}(\mathbf{\hat{r}}_\mu)
  !>                       \tilde{W}_{\text{ex}}_{A/B,\mu\nu}
  !>                           \bar{u}_{\mathbf{k'} j_o}(\mathbf{\hat{r}}_\nu)
  !>                               u_{\mathbf{k'} j_u}(\mathbf{\hat{r}}_\nu)
  !>                                   X(\mathbf{k'} j_o j_u),
  !> \]
  !> where \( n_o, n_u \) are the number of occupied and unoccupied bands considered for the solution of
  !> the BSE respectively. Regrouping the sum reveals the efficient action on the vector and leads to the
  !> expression that is evaluated in this routine:
  !> \[
  !>   \frac{1} {N_k}
  !>       \sum_{\mu=1}^{N_\mu^V}
  !>           \bar{u}_{\mathbf{k} i_u}(\mathbf{\hat{r}}_\mu) u_{\mathbf{k} i_o}(\mathbf{\hat{r}}_\mu)
  !>               \left\{ \sum_{\nu=1}^{N_\mu^V} \tilde{W}_{\text{ex}}_{A/B,\mu\nu} \left[ \sum_{\mathbf{k'}}
  !>                   \left( \sum_{j_u=1}^{n_u} u_{\mathbf{k'} j_u}(\mathbf{\hat{r}}_\nu) \left(\sum_{j_o=1}^{n_o}
  !>                       \bar{u}_{\mathbf{k'} j_o}(\mathbf{\hat{r}}_\nu)
  !>                           X(\mathbf{k'} j_o j_u) \right) \right) \right] \right\}.
  !> \]
  subroutine times_vector(this, vector)
    class(vexc_isdf_kernel_type), target, intent(inout) :: this
    complex(dp), contiguous, target, intent(inout) :: vector(:)

    integer :: ik

    this%D = vector

    ! Evaluate the innermost round bracket
    do ik=1, this%transitions%n_k
      call matrix_multiply(conjg(this%u_o_mu_(ik)%p), this%D_(ik)%p, this%A_(ik)%p, trans_B='T')
    end do

    ! Evaluate squared bracket (sum over unoccpied states and sum over k-points)
    this%B = sum(this%u_u_mu * this%A, dim=2)

    ! Evaluate curly bracket
    ! TODO: Hermetian matrix multiplication?
    call matrix_multiply(this%V_tilde, this%B, this%C)

    ! Evaluate result
    
    do ik=1, this%transitions%n_k
      call matrix_multiply(this%u_u_mu_(ik)%p, spread(this%C, 2, this%transitions%n_o_per_k(ik)) * this%u_o_mu_(ik)%p, this%D_(ik)%p, trans_A='C')
    end do


    vector = this%D

  end subroutine times_vector


  

end module vexc_isdf_kernel