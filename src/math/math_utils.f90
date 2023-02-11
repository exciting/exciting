!> Math utilities and functions
module math_utils
  use iso_fortran_env, only: error_unit
  use, intrinsic :: ISO_C_BINDING

  use precision, only: sp, dp
  use constants, only: pi, zzero, zone, zi, fourpi
  use asserts, only: assert
  use seed_generation, only: set_seed

  implicit none

  private
  public :: identity_integer, &
            identity_real_dp, &
            identity_complex_dp, &
            all_close, &
            all_zero, &
            diag, &
            is_square, &
            is_hermitian, &
            is_unitary, &
            kronecker_product, &
            determinant, &
            permanent, &
            mod1, &
            random_order, &
            round_down, &
            calculate_all_vector_distances, &
            calculate_all_vector_differences, &
            outer_sum, &
            get_subinterval_indices, &
            boundary_mask, &
            is_positive_definite, &
            fractional_part, &
            integer_part, &
            plane_wave_in_spherical_harmonics, &
            get_degeneracies

  !> Default tolerance
  real(dp), parameter :: default_tol = 1e-10

  !>  Return the diagonal of a 2D array
  interface diag
    module procedure diag_int_sp, diag_real_dp, diag_complex_dp
  end interface diag

  !>  Check if a matrix is square
  interface is_square
    module procedure is_square_real_dp, is_square_complex_dp, is_square_integer
  end interface is_square

  !> Check if a matrix is hermitian (symmetric)
  interface is_hermitian
    module procedure is_symmetric_real_dp, is_hermitian_complex_dp
  end interface is_hermitian

  !> Check if a matrix is unitary (orthogonal), such that
  !> \[
  !>    \mathbf{A}^{-1} = \mathbf{A}^\dagger.
  !> \]
  interface is_unitary
    module procedure is_orthogonal_integer, is_orthogonal_real_dp, is_unitary_complex_dp
  end interface is_unitary

  !>  Check if two arrays are close, within an absolute tolerance
  interface all_close
    module procedure all_close_rank0_real_dp, all_close_rank1_real_dp, &
                   & all_close_rank2_real_dp, all_close_rank3_real_dp,&
                     all_close_rank0_complex_dp,&
                   & all_close_rank1_complex_dp, all_close_rank2_complex_dp,&
                   all_close_rank3_complex_dp, all_close_rank4_complex_dp
  end interface all_close

  !>  Check if an array is zero, to within an absolute tolerance
  interface all_zero
    module procedure all_zero_rank0_real_dp, all_zero_rank1_real_dp, &
                   & all_zero_rank2_real_dp, all_zero_rank0_complex_dp,&
                   & all_zero_rank1_complex_dp, all_zero_rank2_complex_dp
  end interface all_zero

  !> Calculate the Kronecker product for three vectors
  !>
  !> The Kronecker product for two vectors is defined as:
  !> \[
  !>    \boldsymbol{a} \bigotimes \boldsymbol{b} =
  !>      \begin{pmatrix}
  !>         a_1\boldsymbol{b}\\ \vdots \\ a_n\boldsymbol{b}
  !>      \end{pmatrix}
  !> \]
  interface kronecker_product
     module procedure kronecker_product_integer, kronecker_product_real_dp, kronecker_product_complex_dp
  end interface kronecker_product

  !> Calculate the determinant of a matrix using Laplace extension.
  !> \[
  !>     \det A = \sum_{i=1}^n (-1)^{i+j}a_ij \cdot \det A_{ij}
  !> \]
  !> where \(A_{ij}\) is derived by canceling the \(i\)'th row and the \(j\)'th collumn of \(A\).
  !> This implementation is based on the reference found at
  !> [rossetta code](http://rosettacode.org/wiki/Determinant_and_permanent#Fortran)
  interface determinant
    module procedure integer_determinant, real_determinant_dp, complex_determinant_dp
  end interface determinant

  !> Calculate the permanent of a matrix using Laplace extension.
  !> \[
  !>     \text{per} A = \sum_{i=1}^n a_ij \cdot \text{per} A_{ij}
  !> \]
  !> where \(A_{ij}\) is derived by canceling the \(i\)'th row and the \(j\)'th collumn of \(A\).
  !> This implementation is based on the reference found at
  !> [rossetta code](http://rosettacode.org/wiki/Determinant_and_permanent#Fortran)
  interface permanent
    module procedure integer_permanent, real_permanent_dp, complex_permanent_dp
  end interface permanent

  !> Modulus after floor division, returning in the range \((0,N]\). Works like the modulus function but instead of
  !> \( \text{mod}(M, N) = 0 \), it returns \( \text{mod1}(M, N) = N \).
  interface mod1
    module procedure mod1_without_offset, mod1_with_offset
  end interface mod1

  !> For a vector of dimensions \( (N_1, ..., N_d) \), setup a \(d\)-rank logical mask with those dimension, such that either 
  !> the upper or the lower boundaries are set to false. The interface supports \(d=3\).
  !> 
  !> The lower boundaries are defined as the first subslabs for each dimension:
  !> `boundary_mask(1, :, ..., :) = .false., boundary_mask(:, 1, :, ..., :) = .false., boundary_mask(:, ..., :, 1) = .false.`
  !>
  !> The upper boundaries are defined as the last subslabs for each dimension:
  !> `boundary_mask(N(1), :, ..., :) = .false., boundary_mask(:, N(2), :, ..., :) = .false., boundary_mask(:, ..., :, N(d)) = .false.`
  !>
  !> By default, the upper boundaries are set to false.
  interface boundary_mask
    module procedure boundary_mask_3d
  end interface boundary_mask

  !> Calculate the fractional part of \(x\), where negative numbers are treated the same as negative
  !> Numbers:
  !> \[
  !>    x - \lfloor x \rfloor.
  !> \]
  !> The routine allows for defining an offset \( c \). If specified the function returns
  !> \[
  !>   x - \lfloor x - c \rfloor
  !> \]
  !> Supported interfaces for scalar vector and matrix input. For matrix input the offset is expected
  !> as vector. For the \(i\)'th row of the input matrix, the \(i\)'th element of the offset vector is
  !> taken as offset.
  interface fractional_part
    module procedure :: fractional_part_scalar, fractional_part_matrix
  end interface fractional_part


  !> Calculate the integer part of of a real number \( x \):
  !> \[
  !>   \lfloor x \rfloor.
  !> \]
  !> The routine allows for defining an offset \( c \). Then it returns
  !> \[
  !>    \lfloor x - c \rfloor.
  !> \]
  !> If the result is close to one, with respect to a tolerance that can be defined,
  !> the result as given by the intrinsic `floor` function is increased by 1.
  interface integer_part
    module procedure :: integer_part_scalar, integer_part_matrix
  end interface integer_part

contains

! identity_real_dp, identity_complex_dp
!
! Setup identity matrix

  !> Setup integer identity matrix.
  function identity_integer(N) result(identity)
    !> Dimension of the identity
    integer, intent(in) :: N

    integer, allocatable :: identity(:, :)

    integer :: i

    allocate(identity(N, N), source=0)

    do i=1, N
      identity(i, i) = 1
    end do
  end function identity_integer

  !> Setup real identity matrix.
  pure function identity_real_dp(N) result(identity)
    !> Dimension of the identity
    integer, intent(in) :: N

    real(dp), allocatable :: identity(:, :)

    integer :: i

    allocate(identity(N, N), source=0.0_dp)

    do i=1, N
      identity(i, i) = 1.0_dp
    end do
  end function identity_real_dp


  !> Setup complex identity matrix.
  pure function identity_complex_dp(N) result(identity)
    !> Dimension of the identity
    integer, intent(in) :: N

    complex(dp), allocatable :: identity(:, :)

    integer :: i

    allocate(identity(N, N), source=zzero)

    do i=1, N
      identity(i, i) = zone
    end do
  end function identity_complex_dp


! diag
!
! Extract the diagonal of a N x N matrix.

  !> Extract the diagonal of an integer \(N \times N\) matrix as vector of size N.
  function diag_int_sp(a) result(diagonal)

    integer(sp), intent(in) :: a(:, :)
    integer(sp), allocatable :: diagonal(:)

    integer(sp) :: i

   call assert(size(a, dim=1) == size(a, dim=2), &
     & 'diag_int_sp: Input needs to be a square matrix.')

   allocate (diagonal(size(a, dim=1)))

   do i = 1, size(diagonal)
     diagonal(i) = a(i, i)
   end do

  end function diag_int_sp


  !> Extract the diagonal of a real \(N \times N\) matrix as vector of size N.
  function diag_real_dp(a) result(diagonal)

    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable :: diagonal(:)

    integer(sp) :: i

    call assert(size(a, dim=1) == size(a, dim=2), &
      & 'diag_real_dp: Input needs to be a square matrix.')

    allocate (diagonal(size(a, dim=1)))

    do i = 1, size(diagonal)
      diagonal(i) = a(i, i)
    end do

  end function diag_real_dp


  !> Extract the diagonal of a complex \(N \times N\) matrix as vector of size N.
  function diag_complex_dp(a) result(diagonal)

    complex(dp), intent(in) :: a(:, :)
    complex(dp), allocatable :: diagonal(:)
    integer(sp) :: i

    call assert(size(a, dim=1) == size(a, dim=2), &
      & 'diag_complex_dp: Input needs to be a square matrix.')

    allocate (diagonal(size(a, dim=1)))

    do i = 1, size(diagonal)
      diagonal(i) = a(i, i)
    end do

  end function diag_complex_dp


! is_square
!
! Check if a matrix is suqare

  !> Check if a real matrix is square
  pure function is_square_real_dp(a) result(square)
    !> Input matrix
    real(dp), intent(in) :: a(:, :)
    !> Is the matrix square
    logical :: square
    square = size(a, dim=1) == size(a, dim=2)
  end function


  !> Check if a complex matrix is square
  pure function is_square_complex_dp(a) result(square)
    !> Input matrix
    complex(dp), intent(in) :: a(:, :)
    !> Is the matrix square
    logical :: square
    square = size(a, dim=1) == size(a, dim=2)
  end function


  !> Check if a integer matrix is square
  pure function is_square_integer(a) result(square)
    !> Input matrix
    integer, intent(in) :: a(:, :)
    !> Is the matrix square
    logical :: square
    square = size(a, dim=1) == size(a, dim=2)
  end function


! is_hermitian
!
! Check if a matrix is hermitian or symmetric respectively.

  !> Check if a real matrix is symmetric.
  logical function is_symmetric_real_dp(A, tol)
    !> Matrix to check for
    real(dp), intent(in) :: A(:, :)
    !> Tolerance for defining equiv
    real(dp), intent(in), optional :: tol

    real(dp) tol_

    tol_ = default_tol
    if(present(tol)) tol_ = tol

    if (.not. is_square(A)) then
      is_symmetric_real_dp = .false.
    else
      is_symmetric_real_dp = all_close(A, transpose(A), tol_)
    end if
  end function is_symmetric_real_dp


  !> Check if a complex matrix is hermitian.
  logical function is_hermitian_complex_dp(A, tol)
    !> Matrix to check for
    complex(dp), intent(in) :: A(:, :)
    !> Tolerance for defining equiv
    real(dp), intent(in), optional :: tol

    real(dp) tol_

    tol_ = default_tol
    if(present(tol)) tol_ = tol

    if (.not. is_square(A)) then
      is_hermitian_complex_dp = .false.
    else
      is_hermitian_complex_dp = all_close(A, transpose(conjg(A)), tol_)
    end if
  end function is_hermitian_complex_dp

! is_unitary
!
! Check if a matrix is unitary or orthogonal respectively.

  !> Check if an integer matrix is orthogonal
  !> The definition of orthogonality is extended to rectangular matrices such that
  !> a matrix \(\mathbf{A} \in \mathbb{Z}^{m \times n}\), if \(m \le n\), is understood to be orthogonal if
  !> \[
  !>   \mathbf{A} \cdot \mathbf{A}^T = mathbf{I}_m
  !> \]
  !> and if \(n > m\)
  !> \[
  !>   \mathbf{A}^T \cdot \mathbf{A} = mathbf{I}_n,
  !> \]
  !> where \(\mathbf{I}_m\) is the \(m \times m\) dimensional identity matrix.
  logical function is_orthogonal_integer(A)
    !> Matrix to check for
    integer, intent(in) :: A(:, :)

    integer, allocatable :: A_times_A_T(:, :)

    integer :: m, n

    m = size(A, dim=1)
    n = size(A, dim=2)

    if(m <= n) then
      A_times_A_T = matmul(A, transpose(A))
    else
      A_times_A_T = matmul(transpose(A), A)
    end if

    is_orthogonal_integer = all(A_times_A_T == identity_integer(size(A_times_A_T, dim=1)))
  end function is_orthogonal_integer

  !> Check if an real matrix is orthogonal
  !> The definition of orthogonality is extended to rectangular matrices such that
  !> a matrix \(\mathbf{A} \in \mathbb{R}^{m \times n}\), if \(m \le n\), is understood to be orthogonal if
  !> \[
  !>   \mathbf{A} \cdot \mathbf{A}^T = mathbf{I}_m
  !> \]
  !> and if \(n > m\)
  !> \[
  !>   \mathbf{A}^T \cdot \mathbf{A} = mathbf{I}_n,
  !> \]
  !> where \(\mathbf{I}_m\) is the \(m \times m\) dimensional identity matrix.
  logical function is_orthogonal_real_dp(A, tol)
    !> Matrix to check for
    real(dp), intent(in) :: A(:, :)
    !> Tolerance for defining equiv
    real(dp), intent(in), optional :: tol

    real(dp) tol_
    real(dp), allocatable :: A_times_A_T(:, :)
    integer :: m, n

    tol_ = default_tol
    if(present(tol)) tol_ = tol

    m = size(A, dim=1)
    n = size(A, dim=2)

    if(m <= n) then
      A_times_A_T = matmul(A, transpose(A))
    else
      A_times_A_T = matmul(transpose(A), A)
    end if

    is_orthogonal_real_dp = all_close(A_times_A_T, identity_real_dp(size(A_times_A_T, dim=1)), tol_)
  end function is_orthogonal_real_dp


  !> Check if a complex matrix is unitary
  !> The definition of unitarity is extended to rectangular matrices such that
  !> a matrix \(\mathbf{A} \in \mathbb{C}^{m \times n}\), if \(m \le n\), is understood to be orthogonal if
  !> \[
  !>   \mathbf{A} \cdot \mathbf{A}^T = mathbf{I}_m
  !> \]
  !> and if \(n > m\)
  !> \[
  !>   \mathbf{A}^T \cdot \mathbf{A} = mathbf{I}_n,
  !> \]
  !> where \(\mathbf{I}_m\) is the \(m \times m\) dimensional identity matrix.
  logical function is_unitary_complex_dp(A, tol)
    !> Matrix to check for
    complex(dp), intent(in) :: A(:, :)
    !> Tolerance for defining equiv
    real(dp), intent(in), optional :: tol

    real(dp) tol_
    complex(dp), allocatable :: A_times_A_C(:, :)
    integer :: m, n

    tol_ = default_tol
    if(present(tol)) tol_ = tol

    m = size(A, dim=1)
    n = size(A, dim=2)

    if(m <= n) then
      A_times_A_C = matmul(A, transpose(conjg(A)))
    else
      A_times_A_C = matmul(transpose(conjg(A)), A)
    end if

    is_unitary_complex_dp = all_close(A_times_A_C, identity_complex_dp(size(A_times_A_C, dim=1)), tol_)
  end function is_unitary_complex_dp

! all_close
!
! Check if two scalars, vectors or matrices are close to each other element wise
! with respect to a certain tolerance
  
  
  !> Check if two real scalars \( a \) and \( b \) are equal,
  !> where equal is defined as \( |a - b| \leq abs\_tol \).
  logical function all_close_rank0_real_dp(a, b, tol)
    !> Input array
    real(dp), intent(in) :: a
    !> Reference array
    real(dp), intent(in) :: b
    !> Absolute tolerance for input and reference to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_close_rank0_real_dp = abs(a - b) <= tol_
  end function all_close_rank0_real_dp


  !> Check if two real rank-1 arrays \( \mathbf{a} \) and \( \mathbf{b} \)
  !> are equal, where equal is defined as
  !> \( |a_i - b_i| \leq abs\_tol,  \forall i \).
  !> As such, the tolerance is checked elementwise.
  logical function all_close_rank1_real_dp(a, b, tol)
    !> Input array
    real(dp), intent(in) :: a(:)
    !> Reference array
    real(dp), intent(in) :: b(:)
    !> Absolute tolerance for input and reference to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    call assert(size(a) == size(b), &
      & 'all_close_rank1_real_dp: size of input arrays differs.')

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_close_rank1_real_dp = all(abs(a - b) <= tol_)
  end function all_close_rank1_real_dp


  !> Check if two real rank-2 arrays \( \mathbf{a} \) and \( \mathbf{b} \)
  !> are equal, where equal is defined as
  !> \( |a_{ij} - b_{ij}| \leq abs\_tol,  \forall i,j \).
  !> As such, the tolerance is checked elementwise.
  logical function all_close_rank2_real_dp(a, b, tol)
    !> Input array
    real(dp), intent(in) :: a(:,:)
    !> Reference array
    real(dp), intent(in) :: b(:,:)
    !> Absolute tolerance for input and reference to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    call assert(size(a) == size(b), &
      & 'all_close_rank2_real_dp: size of input arrays differs.')

    call assert(all(shape(a) == shape(b)), &
      & 'all_close_rank2_real_dp: shape of input arrays differs.')

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_close_rank2_real_dp = all(abs(a - b) <= tol_)
  end function all_close_rank2_real_dp


  !> Check if two real rank-3 arrays \( \mathbf{a} \) and \( \mathbf{b} \)
  !> are equal, where equal is defined as
  !> \( |a_{ijk} - b_{ijk}| \leq abs\_tol,  \forall i,j,k \).
  !> As such, the tolerance is checked elementwise.
  logical function all_close_rank3_real_dp(a, b, tol)

    !> Input array
    real(dp), intent(in) :: a(:, :, :)
    !> Reference array
    real(dp), intent(in) :: b(:, :, :)
    !> Absolute tolerance for input and reference to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    call assert(size(a) == size(b), &
      & 'all_close_rank3_real_dp: size of input arrays differs.')

    call assert(all(shape(a) == shape(b)), &
      & 'all_close_rank3_real_dp: shape of input arrays differs.')

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_close_rank3_real_dp = all(abs(a - b) <= tol_)

  end function all_close_rank3_real_dp


  !> Check if two complex scalars \( a \) and \( b \)
  !> are equal, where equal is defined as \( |a - b| \leq abs\_tol \).
  !> The abs of a complex scalar \( a \) is defined as:
  !> \[
  !>    |a_i| = \sqrt{ a_i \cdot {a_i}^*}.
  !> \]
  !> As such, the tolerance is a real value.
  logical function all_close_rank0_complex_dp(a, b, tol)
    !> Input array
    complex(dp), intent(in) :: a
    !> Reference array
    complex(dp), intent(in) :: b
    !> Absolute tolerance for input and reference to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_close_rank0_complex_dp = abs(a - b) <= tol_
  end function all_close_rank0_complex_dp


  !> Check if two complex rank-1 arrays \( \mathbf{a} \) and \( \mathbf{b} \)
  !>  are equal, where equal is defined as
  !>  \( |a_i - b_i| \leq abs\_tol,  \forall i \).
  !> The abs of a complex scalar \( a_i \) is defined as:
  !> \[
  !>    |a_i| = \sqrt{ a_i \cdot {a_i}^*}.
  !> \]
  !> As such, the tolerance is a real value and is
  !> checked elementwise.
  logical function all_close_rank1_complex_dp(a, b, tol)
    !> Input array
    complex(dp), intent(in) :: a(:)
    !> Reference array
    complex(dp), intent(in) :: b(:)
    !> Absolute tolerance for input and reference to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    call assert(size(a) == size(b), &
      & 'all_close_rank1_complex_dp: size of input arrays differs.')

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_close_rank1_complex_dp = all(abs(a - b) <= tol_)
  end function all_close_rank1_complex_dp


  !> Check if two rank-2 arrays \( \mathbf{a} \) and \( \mathbf{b} \)
  !> of are equal, where equal is defined as
  !>  \( |a_{ij} - b_{ij}| \leq abs\_tol,  \forall i,j \).
  !> The abs of a complex scalar \( a_i \) is defined as:
  !> \[
  !>    |a_i| = \sqrt{ a_i \cdot {a_i}^*}.
  !> \]
  !> As such, the tolerance is a real value and is
  !> checked elementwise.
  logical function all_close_rank2_complex_dp(a, b, tol)
    !> Input array
    complex(dp), intent(in) :: a(:,:)
    !> Reference array
    complex(dp), intent(in) :: b(:,:)
    !> Absolute tolerance for input and reference to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    call assert(size(a) == size(b), &
      & 'all_close_rank2_complex_dp: size of input arrays differs.')

    call assert(all(shape(a) == shape(b)), &
      & 'all_close_rank2_complex_dp: shape of input arrays differs.')

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_close_rank2_complex_dp = all(abs(a - b) <= tol_)
  end function all_close_rank2_complex_dp

  !> Check if two rank-3 arrays \( \mathbf{a} \) and \( \mathbf{b} \)
  !> of are equal, where equal is defined as
  !>  \( |a_{ijk} - b_{ijk}| \leq abs\_tol,  \forall i,j,k \).
  !> As such, the tolerance is a real value and is
  !> checked elementwise.
  logical function all_close_rank3_complex_dp(a, b, tol)

    !> Input array
    complex(dp), intent(in) :: a(:, :, :)
    !> Reference array
    complex(dp), intent(in) :: b(:, :, :)
    !> Absolute tolerance for input and reference to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    call assert(size(a) == size(b), &
      & 'all_close_rank3_complex_dp: size of input arrays differs.')

    call assert(all(shape(a) == shape(b)), &
      & 'all_close_rank3_complex_dp: shape of input arrays differs.')

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_close_rank3_complex_dp = all(abs(a - b) <= tol_)

  end function all_close_rank3_complex_dp

  !> Check if two rank-4 arrays \( \mathbf{a} \) and \( \mathbf{b} \)
  !> of are equal, where equal is defined as
  !>  \( |a_{ijkl} - b_{ijkl}| \leq abs\_tol,  \forall i,j,k,l \).
  !> As such, the tolerance is a real value and is
  !> checked elementwise.
  logical function all_close_rank4_complex_dp(a, b, tol)

    !> Input array
    complex(dp), intent(in) :: a(:, :, :, :)
    !> Reference array
    complex(dp), intent(in) :: b(:, :, :, :)
    !> Absolute tolerance for input and reference to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    call assert(size(a) == size(b), &
      & 'all_close_rank4_complex_dp: size of input arrays differs.')

    call assert(all(shape(a) == shape(b)), &
      & 'all_close_rank4_complex_dp: shape of input arrays differs.')

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_close_rank4_complex_dp = all(abs(a - b) <= tol_)

  end function all_close_rank4_complex_dp

  
! all_zero
!
! Check if a sclalar, a vector or a matrix is element wise close to zero.

  !> Check if a complex scalar \( a \) is zero,
  !> where zero is defined as \( |a - b| \leq abs\_tol \).
  !> The abs of a complex scalar \( a \) is defined as:
  !> \[
  !>    |a| = \sqrt{ a \cdot {a}^*}.
  !> \]
  !> As such, the tolerance is a real value.
  pure logical function all_zero_rank0_complex_dp(a, tol)
    !> Input array
    complex(dp), intent(in) :: a
    !> Absolute tolerance for input and reference to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_zero_rank0_complex_dp = abs(a ) <= tol_
  end function all_zero_rank0_complex_dp


  !> Check if a complex rank-1 array \( \mathbf{a} \) is zero,
  !>  where zero is defined as \( |a_i| \leq abs\_tol,  \forall i \).
  !> The abs of a complex scalar \( a_i \) is defined as:
  !> \[
  !>    |a_i| = \sqrt{ a_i \cdot {a_i}^*}.
  !> \]
  !> As such, the tolerance is a real value and is
  !> checked elementwise.
  pure logical function all_zero_rank1_complex_dp(a, tol)
    !> Input array
    complex(dp), intent(in) :: a(:)
    !> Absolute tolerance for input and reference to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_zero_rank1_complex_dp = all(abs(a) <= tol_)
  end function all_zero_rank1_complex_dp


  !> Check if a complex rank-2 array \( \mathbf{a} \) is zero,
  !> where zero  is defined as \( |a_{ij}| \leq abs\_tol,  \forall i,j \).
  !> The abs of a complex scalar \( a_i \) is defined as:
  !> \[
  !>    |a_i| = \sqrt{ a_i \cdot {a_i}^*}.
  !> \]
  !> As such, the tolerance is a real value and is
  !> checked elementwise.
  pure logical function all_zero_rank2_complex_dp(a, tol)
    !> Input array
    complex(dp), intent(in) :: a(:,:)
    !> Absolute tolerance for input and reference to be considered equal

    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_zero_rank2_complex_dp = all(abs(a) <= tol_)
  end function all_zero_rank2_complex_dp


  !> Check if a real scalar \(a \) is zero, where zero
  !> is defined as \( |a| \leq abs\_tol \).
  pure logical function all_zero_rank0_real_dp(a, tol)
    !> Input array
    real(dp), intent(in) :: a
    !> Absolute tolerance for input to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_zero_rank0_real_dp = abs(a) <= tol_
  end function all_zero_rank0_real_dp


  !> Check if a real rank-1 array \( \mathbf{a} \) is zero,
  !> where zero  is defined as \( |a_i| \leq abs\_tol,  \forall i \).
  !> As such, the tolerance is checked elementwise.
  pure logical function all_zero_rank1_real_dp(a, tol)
    !> Input array
    real(dp), intent(in) :: a(:)
    !> Absolute tolerance for input to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_zero_rank1_real_dp = all(abs(a) <= tol_)
  end function all_zero_rank1_real_dp


  !> Check if a real rank-2 array \( \mathbf{a} \) is zero,
  !> where zero  is defined as \( |a_{ij}| \leq abs\_tol,  \forall i,j \).
  !> As such, the tolerance and is checked elementwise.
  pure logical function all_zero_rank2_real_dp(a, tol)
    !> Input array
    real(dp), intent(in) :: a(:,:)
    !> Absolute tolerance for input to be considered equal
    real(dp), intent(in), optional :: tol

    !> Local absolute tolerance
    real(dp) :: tol_

    tol_ = default_tol
    if (present(tol)) tol_ = tol

    all_zero_rank2_real_dp = all(abs(a) <= tol_)
  end function all_zero_rank2_real_dp


! kronecker_product
!
! Calculate the Kronecker product for three vectors.

  !> Calculate the Kronecker product for three integer vectors with various dimensions:
  !> \( \mathbf{A} \bigotimes \mathbf{B} \bigotimes \mathbf{C} \)
  !>
  !> The Kronecker product for two vectors is defined as:
  !> \[
  !>    \mathbf{A} \bigotimes \mathbf{B} =
  !>           \begin{pmatrix}
  !>               a_1\mathbf{B}\\ \vdots \\ a_n\mathbf{B}
  !>           \end{pmatrix}
  !> \]
  function kronecker_product_integer(A, B, C) result(output_vector)
    !> first input vector A
    integer, intent(in) :: A(:), &
      !> second input vector B
      B(:), &
      !> third input vector C
      C(:)
    !> result
    integer, allocatable :: output_vector(:)
    ! local variables
    integer :: i, NA, NB, NC, ilow, iup
    integer, allocatable :: work(:)

    NA = size(A)
    NB = size(B)
    NC = size(C)

    allocate (work(NA * NB))
    do i = 1, NA
      ilow = (i - 1) * NB + 1
      iup = ilow - 1 + NB
      work(ilow: iup) = A(i) * B
    end do

    allocate (output_vector(NA * NB * NC))
    do i = 1, size(A)*size(B)
      ilow = (i - 1) * NC + 1
      iup = ilow - 1 + NC
      output_vector(ilow: iup) = work(i) * C
    end do
  end function kronecker_product_integer


  !> Calculate the Kronecker product for three real vectors with various dimensions:
  !> \( \mathbf{A} \bigotimes \mathbf{B} \bigotimes \mathbf{C} \)
  !>
  !> The Kronecker product for two vectors is defined as:
  !> \[
  !>    \mathbf{A} \bigotimes \mathbf{B} =
  !>           \begin{pmatrix}
  !>               a_1\mathbf{B}\\ \vdots \\ a_n\mathbf{B}
  !>           \end{pmatrix}
  !> \]
  function kronecker_product_real_dp(A, B, C) result(output_vector)
    !> first input vector A
    real(dp), intent(in) :: A(:), &
                            !> second input vector B
                            B(:), &
                            !> third input vector A
                            C(:)
    !> result
    real(dp), allocatable :: output_vector(:)
    ! local variables
    integer :: i, NA, NB, NC, ilow, iup
    real(dp), allocatable :: work(:)

    NA = size(A)
    NB = size(B)
    NC = size(C)

    allocate (work(NA * NB))
    do i = 1, NA
      ilow = (i - 1) * NB + 1
      iup = ilow - 1 + NB
      work(ilow: iup) = A(i) * B
    end do

    allocate (output_vector(NA * NB * NC))
    do i = 1, size(A)*size(B)
      ilow = (i - 1) * NC + 1
      iup = ilow - 1 + NC
      output_vector(ilow: iup) = work(i) * C
    end do
  end function kronecker_product_real_dp


  !> Calculate the Kronecker product for three complex vectors with various dimensions,
  !> \( \mathbf{A} \bigotimes \mathbf{B} \bigotimes \mathbf{C} \).
  !>
  !> The Kronecker product for two vectors is defined as:
  !> \[
  !>     \mathbf{A} \bigotimes \mathbf{B} =
  !>        \begin{pmatrix} a_1\mathbf{B}\\ \vdots \\ a_n\mathbf{B}
  !>     \end{pmatrix}
  !> \]
  function kronecker_product_complex_dp(A, B, C) result(output_vector)
    !> first input vector A
    complex(dp), intent(in) :: A(:), &
                               !> second input vector B
                               B(:), &
                               !> third input vector A
                               C(:)
    !> result
    complex(dp), allocatable :: output_vector(:)
    ! local variables
    integer :: i, NA, NB, NC, ilow, iup
    complex(dp), allocatable :: work(:)

    NA = size(A)
    NB = size(B)
    NC = size(C)

    allocate (work(NA * NB))
    do i = 1, NA
      ilow = (i - 1) * NB + 1
      iup = ilow - 1 + NB
      work(ilow: iup) = A(i) * B
    end do

    allocate (output_vector(NA * NB * NC))
    do i = 1, size(A)*size(B)
      ilow = (i - 1) * NC + 1
      iup = ilow - 1 + NC
      output_vector(ilow: iup) = work(i) * C
    end do
  end function kronecker_product_complex_dp


! determinant
!
! Calculate the determinant of a matrix. See determinent_laplace for more information.

  !> Calculates the determinant of a real matrix using Laplace extension.
  real(dp) function real_determinant_dp(A)
    !> Matrix for which the determinant is calculated
    real(dp), dimension(:, :), intent(in) :: A

    call assert(is_square(A), &
      & message='real_determinant_dp: Input needs to be a square matrix.')
    real_determinant_dp = real_determinant_laplace_dp(A, -1)
  end function real_determinant_dp


  !> Calculates the determinant of a complex matrix using Laplace extension.
  complex(dp) function complex_determinant_dp(A)

    !> Matrix for which the determinant is calculated
    complex(dp), dimension(:, :), intent(in) :: A

    call assert(is_square(A), &
      & message='complex_determinant_dp: Input needs to be a square matrix.')
    complex_determinant_dp = complex_determinant_laplace_dp(A, -1)
  end function complex_determinant_dp


  !> Calculates the determinant of a integer matrix using Laplace extension.
  integer function integer_determinant(A)

    !> Matrix for which the determinant is calculated
    integer, dimension(:, :), intent(in) :: A

    call assert(is_square(A), &
      & message='integer_determinant: Input needs to be a square matrix.')
    integer_determinant = integer_determinant_laplace(A, -1)
  end function integer_determinant


! permanent
!
! Calculate the determinant of a matrix. See determinent_laplace for more information.

  !> Calculates the permanent of a real matrix using Laplace extension.
  real(dp) function real_permanent_dp(A)
    !> Matrix for which the permanent is calculated
    real(dp), dimension(:, :), intent(in) :: A

    call assert(is_square(A), &
      &  message='real_permanent_dp: Input needs to be a square matrix.')
    real_permanent_dp = real_determinant_laplace_dp(A, 1)
  end function real_permanent_dp



  !> Calculates the permanent of a complex matrix using Laplace extension.
  complex(dp) function complex_permanent_dp(A)
    !> Matrix for which the permanent is calculated
    complex(dp), dimension(:, :), intent(in) :: A

    call assert(is_square(A), &
      &  message='complex_permanent_dp: Input needs to be a square matrix.')
    complex_permanent_dp = complex_determinant_laplace_dp(A, 1)
  end function complex_permanent_dp


  !> Calculates the permanent of a integer matrix using Laplace extension.
  integer function integer_permanent(A)
    !> Matrix for which the permanent is calculated
    integer, dimension(:, :), intent(in) :: A

    call assert(is_square(A), &
      &  message='integer_permanent: Input needs to be a square matrix.')
    integer_permanent = integer_determinant_laplace(A, 1)
  end function integer_permanent


! determinant_laplace
!
! Calculate the determinant of a real matrix using Laplace extension.

  !> Calculate the determinant of a real matrix using Laplace extension.
  !> \[
  !>     \det A = \sum_{i=1}^n (-1)^{i+j}a_ij \cdot \det A_{ij}
  !> \]
  !> where $A_{ij}$ is derived by canceling the $i$'th row and the $j$'th collumn of $A$.
  !>
  !> This implementation is based on the reference found at
  !> [rossetta code](http://rosettacode.org/wiki/Determinant_and_permanent#Fortran)
  recursive function real_determinant_laplace_dp(a, permanent) result(accumulation)
    !> Matrix for which the determinant is calculated
    real(dp), dimension(:, :), intent(in) :: a

    !> Setting permanent to 1 computes the permanent.
    !> Setting permanent to -1 computes the determinant.
    integer, intent(in), optional :: permanent
    real(dp), allocatable :: b(:,:)
    integer :: permanent_, i, sgn
    real(dp) :: accumulation
    !> Dimension of input matrix
    integer(sp) :: n

    call assert(is_square(a), &
      &  message='real_determinant_laplace_dp: Input needs to be a square matrix.')

    n = size(a, dim=2)

    allocate(b(n - 1, n - 1))
    permanent_ = -1
    if (present(permanent)) then
      call assert( abs(permanent) == 1, &
             & message="real_determinant_laplace_dp: permanent needs to be 1 or -1.")
      permanent_ = permanent
    end if

    if (n == 1) then
      accumulation = a(1, 1)
    else
      accumulation = 0
      sgn = 1
      do i = 1, n
        b(:, :(i - 1)) = a(2:, :i - 1)
        b(:, i:) = a(2:, i + 1:)
        accumulation = accumulation + sgn * a(1, i) * real_determinant_laplace_dp(b, permanent_)
        sgn = sgn * permanent_
      end do
    end if
  end function real_determinant_laplace_dp


  !> Calculate complex determinant or permanent respectively.
  !> See [[real_determinant_laplace_dp(subroutine)]] for details.
  recursive function complex_determinant_laplace_dp(a, permanent) result(accumulation)
    !> Matrix for which the determinant is calculated
    complex(dp), dimension(:, :), intent(in) :: a
    !> Setting permanent to 1 computes the permanent.
    !> Setting permanent to -1 computes the determinant.
    integer, intent(in), optional :: permanent

    complex(dp), allocatable :: b(:,:)
    integer :: permanent_, i, sgn
    complex(dp) :: accumulation
    !> Dimension of input matrix
    integer(sp) :: n

    call assert(is_square(a), &
      &  message='complex_determinant_laplace_dp: Input needs to be a square matrix.')

    n = size(a, dim=2)

    allocate(b(n - 1, n - 1))
    permanent_ = -1
    if (present(permanent)) then
      call assert( abs(permanent) == 1, &
             & message="complex_determinant_laplace_dp: permanent needs to be 1 or -1.")
      permanent_ = permanent
    end if

    if (n == 1) then
      accumulation = a(1, 1)
    else
      accumulation = 0
      sgn = 1
      do i = 1, n
        b(:, :(i - 1)) = a(2:, :i - 1)
        b(:, i:) = a(2:, i + 1:)
        accumulation = accumulation + sgn * a(1, i) * complex_determinant_laplace_dp(b, permanent_)
        sgn = sgn * permanent_
      end do
    end if
  end function complex_determinant_laplace_dp


  !> Calculate integer determinant or permanent respectively.
  !> See [[real_determinant_laplace_dp(subroutine)]] for details.
  recursive function integer_determinant_laplace(a, permanent) result(accumulation)
    !> Matrix for which the determinant is calculated
    integer, dimension(:, :), intent(in) :: a
    !> Setting permanent to 1 computes the permanent.
    !> Setting permanent to -1 computes the determinant.
    integer, intent(in), optional :: permanent

    integer, allocatable :: b(:,:)
    integer :: permanent_, i, sgn
    integer :: accumulation
    !> Dimension of input matrix
    integer(sp) :: n

    call assert(is_square(a), &
      &  message='integer_determinant_laplace: Input needs to be a square matrix.')

    n = size(a, dim=2)

    allocate(b(n - 1, n - 1))
    permanent_ = -1
    if (present(permanent)) then
      call assert( abs(permanent) == 1, &
             & message="integer_determinant_laplace: permanent needs to be 1 or -1.")
      permanent_ = permanent
    end if

    if (n == 1) then
      accumulation = a(1, 1)
    else
      accumulation = 0
      sgn = 1
      do i = 1, n
        b(:, :(i - 1)) = a(2:, :i - 1)
        b(:, i:) = a(2:, i + 1:)
        accumulation = accumulation + sgn * a(1, i) * integer_determinant_laplace(b, permanent_)
        sgn = sgn * permanent_
      end do
    end if
  end function integer_determinant_laplace

! mod1

  !> Modulus after floor division. Maps an integer number to an interval \((0,N]\).
  !> \[
  !>   \begin{split}
  !>     \text{mod1}(1, 4) &= 1 \\\
  !>     \text{mod1}(4, 4) &= 4 \\\
  !>     \text{mod1}(6, 4) &= 2 \\\
  !>     \text{mod1}(12, 4) &= 4 \\\
  !>     \text{mod1}(0, 4) &= 4 \\\
  !>     \text{mod1}(-1, 4) &= 3 \\\
  !>     \text{mod1}(-1, -4) &= -1 \\\
  !>     \text{mod1}(1, -4) &= -3
  !>   \end{split}
  !> \]
  integer elemental function mod1_without_offset(M, N)
    !> integer to translate
    integer, intent(in) :: M
    !> length of the cycle
    integer, intent(in) :: N

    if (mod(M, N) == 0) then
      mod1_without_offset = N
    elseif(real(M) / N > 0.0) then
      mod1_without_offset = mod(M, N)
    elseif(real(M) / N < 0.0) then
      mod1_without_offset = mod(M, N) + N
    end if
  end function mod1_without_offset

  !> Modulus after floor division with integer offset, returning in the range \((\text{offset}, N + \text{offset}]\):
  !> `mod1_offset(M, N, offset) = mod1(M - offset, N) + offset. Sess [[mod1_]].
  integer elemental function mod1_with_offset(M, N, offset)
    !> integer to translate
    integer, intent(in) :: M
    !> length of the cycle
    integer, intent(in) :: N
    !> Offset of the cycle
    integer, intent(in) :: offset

    mod1_with_offset = mod1_without_offset(M - offset, N) + offset
  end function mod1_with_offset

! random_order

  !> Return an integer vector of length N with all numbers between 1 and N randomly ordered.
  function random_order(N) result (p)
    !> length of the permutation
    integer, intent(in) :: N

    integer :: p(N)

    integer :: j, k
    real(sp) :: u

    p = 0
    do j = 1, N
      call random_number(u)
      k = floor(j * u) + 1
      p(j) = p(k)
      p(k) = j
    end do
  end function random_order

! round down

  !> Rounds a real number down to the nth decimal place, e.g.
  !>
  !> If n is greater than 0, then number is rounded down to the specified number of decimal places.
  !> 1123.234523 -> 1123.234 for n = 3.
  !>
  !> If n is 0, then number is rounded down to the nearest integer.
  !> -2123.77963 -> -2123 for n = 0
  !>
  !> If n is less than 0, then number is rounded down to the left of the decimal point.
  !> 2123.77963 -> 2100 for n = -2.
  elemental function round_down(x, n) result(x_rounded)
    !> Real number to round down
    real(dp), intent(in) :: x
    !> Number of decimal digits to round down to
    integer, intent(in) :: n

    real(dp) :: x_rounded
    x_rounded = dble(int(x*(10.0_dp**n)))/(10.0_dp**n)
  end function round_down

  !> Given two sets of vectors with the same dimension, calculate the Eucledian distances between all vectors of the sets.
  subroutine calculate_all_vector_distances(vector_set1, vector_set2, all_vector_distances)
    !> Vector sets.
    real(dp), intent(in), contiguous :: vector_set1(:, :), vector_set2(:, :)
    !> Matrix that holds the differences, such that `all_vector_distances(i, j) = norm2(vector_set1(:, i) - vector_set2(:, j))`
    real(dp), intent(out), contiguous :: all_vector_distances(:, :)

    integer :: i, j, k, l, m, n
    real(dp) :: vector2(3)

    k = size(vector_set1, 1)
    l = size(vector_set1, 2)
    m = size(vector_set2, 1)
    n = size(vector_set2, 2) 

    call assert(k == m, 'vector_set1 and vector_set2 have not the same number of rows.')
    call assert(all(shape(all_vector_distances) == [l, n]), "all_vector_distances has the wrong shape.")

    do j = 1, n
      vector2 = vector_set2(:, j)
      do i = 1, l
        all_vector_distances(i, j) = norm2(vector_set1(:, i) - vector2)
      end do
    end do

  end subroutine calculate_all_vector_distances

  !> Given two sets of vectors with the same dimension, calculate the differences between all vectors of the sets.
  subroutine calculate_all_vector_differences(vector_set1, vector_set2, all_vector_differences)
    !> Vector sets
    real(dp), intent(in), contiguous :: vector_set1(:, :), vector_set2(:, :)
    !> Differences between the vectors such that `all_vector_differences(:, i, j) = vector_set1(:, i) - vector_set2(:, j)`
    real(dp), intent(out), contiguous :: all_vector_differences(:, :, :)

    integer :: i, j, k, l, m, n
    real(dp) :: vector2(3)

    k = size(vector_set1, 1)
    l = size(vector_set1, 2)
    m = size(vector_set2, 1)
    n = size(vector_set2, 2)

    call assert(k == m, 'vector_set1 and vector_set2 have not the same number of rows.')
    call assert(all(shape(all_vector_differences) == [k, l, n]), 'all_vector_differences has not the correct shape.')

    do j = 1, n
      vector2 = vector_set2(:, j)
      do i = 1, l
        all_vector_differences(:, i, j) = vector_set1(:, i) - vector2
      end do
    end do

  end subroutine calculate_all_vector_differences

  !> Calculate the outer sum of two vectors \( \mathbf{a} \in \mathbb{R}^n \) and \( \mathbf{b} \in \mathbb{R}^m \)
  !> such that the result is a matrix \( \mathbf{C} \in \mathbb{R}^{n \times m} \), given element wise by
  !> \[
  !>   \text{C}_{ij} = \text{a}_{i} + \text{b}_{j}.
  !> \]
  subroutine outer_sum(a, b, C)
    !> Input vector a
    real(dp), intent(in), contiguous :: a(:)
    !> Input vector b
    real(dp), intent(in), contiguous :: b(:)
    !> Output matrix C
    real(dp), intent(out), contiguous :: C(:, :)

    integer :: i, j

    call assert(all(shape(C) == [size(a), size(b)]), "First dimension of output matrix C &
                has to equal the size of input vector a and second dimension of output matrix C &
                has to equal the size of input vector b.")

    do j = 1, size(b)
       do i = 1, size(a)
          C(i, j) = a(i) + b(j)
       end do
    end do

  end subroutine outer_sum


  !> Given the size of an interval \( i \in \mathbb{N} \) and a subinterval size \( s \in \mathbb{N} \),
  !> where \(i\) must be a multiple of \(s\), such that \( n = \frac{i}{s} \in \mathbb{N} \), splits
  !> the interval into \(n\) subintervals of size \(s\) and returns the lower and upper boundary
  !> indice of each subinterval, saved columnwise in the output matrix.
  !>
  !> \[ \text{Example: } i = 6, s = 3, \text{ returns: }
  !> \mathbf{indices} = \left( \begin{smallmatrix}
  !>  1 & 4\\
  !>  3 & 6
  !> \end{smallmatrix} \right)
  !> \]
  function get_subinterval_indices(interval_size, subinterval_size) result(indices)
    !> Interval size
    integer, intent(in) :: interval_size
    !> Subinterval size
    integer, intent(in) :: subinterval_size
    !> Number of multiples
    integer ::  n_multiple

    integer :: i, i1, i2
    !> Output indices matrix
    integer, allocatable :: indices(:, :)

    call assert(interval_size >= subinterval_size, message = "Interval size has to be &
                larger than or equal to the subinterval size.")

    n_multiple = interval_size/ subinterval_size
    allocate(indices(2, n_multiple))

    call assert(n_multiple * subinterval_size == interval_size, message = "n_multiple has to be an &
                integer. Input sizes have to be multiples of each other.")

    do i = 1, n_multiple
      i1 = 1 + (i-1) * subinterval_size
      i2 = i1 + subinterval_size - 1
      indices(:, i) = [i1, i2]
    end do

  end function get_subinterval_indices

! boundary_mask

  !> For a vector of dimensions \( (N_1, N_2, N_3) \), setup a \(3\)-rank logical mask with those dimension, such that either
  !> the upper or the lower boundaries are set to false.
  !>
  !> The lower boundaries are defined as the first subslabs for each dimension:
  !> `boundary_mask(1, :, :) = .false., boundary_mask(:, 1, :) = .false., boundary_mask(:, :, 1) = .false.`
  !>
  !> The upper boundaries are defined as the last subslabs for each dimension:
  !> `boundary_mask(N(1), :, :) = .false., boundary_mask(:, N(2), :) = .false., boundary_mask(:, :, N(3)) = .false.`
  !>
  !> By default, the upper boundaries are set to false.
  function boundary_mask_3d(upper_bounds, use_lower_bounds) result(boundary_mask)
    !> Number of grid points per dimension
    integer :: upper_bounds(3)
    !> Set the lower boundaries to .false. if .true.
    !> If false (default), the upper boundaries are set to .false.
    logical, optional :: use_lower_bounds

    logical, allocatable :: boundary_mask(:, :, :)

    integer :: bounds(3)
    logical :: use_lower_bounds_local

    use_lower_bounds_local = .false.
    if(present(use_lower_bounds)) use_lower_bounds_local = use_lower_bounds

    bounds = upper_bounds
    if(use_lower_bounds_local) bounds = [1, 1, 1]

    allocate(boundary_mask(upper_bounds(1), upper_bounds(2), upper_bounds(3)))
    boundary_mask = .true.

    boundary_mask(bounds(1), :, :) = .false.
    boundary_mask(:, bounds(2), :) = .false.
    boundary_mask(:, :, bounds(3)) = .false.
  end function boundary_mask_3d


  !> Check if a hermitian matrix is positive-definite. This is done by checking
  !> if all the eigenvalues are positive
  logical function is_positive_definite( A )
    !> Matrix to be checked
    complex(dp), intent(in)   :: A(:, :)

    integer                   :: dim, info, lwork
    real(dp), allocatable     :: rwork(:), eigenvalues(:)
    complex(dp), allocatable  :: A_copy(:, :), work(:)
    character(200)            :: error_msg

    call assert( is_hermitian(A), 'A is not hermitian' )
    ! TODO Issue #25: Lapack wrapper for ZHEEV needed
    dim = size( A, 1 )
    allocate( A_copy(dim, dim), eigenvalues(dim), rwork(3*dim-2) )
    A_copy = A
    ! Obtain the optimum lwork
    allocate( work(2) )
    lwork = -1
    call ZHEEV( 'N', 'U', dim, A_copy, dim, eigenvalues, work, lwork, rwork, &
      & info )
    lwork = work(1)
    deallocate( work )
    allocate( work(lwork) )
    ! Obtain the eigenvalues of A
    call ZHEEV( 'N', 'U', dim, A_copy, dim, eigenvalues, work, lwork, rwork, info )
    write(error_msg,*) 'Error(is_positive_definite): ZHEEV returned info = ', info
    call assert( info==0, error_msg )
    is_positive_definite = all( eigenvalues > 0._dp )
  end function


  !> Calculate the fractional part of a real number \(x\):
  !> \[
  !>    x - \lfloor x \rfloor
  !>  \]
  !> Optionally an integer offset \( c \) for the intervall can be defined, such that \( x \) is mapped to
  !> \( [c, c + 1) \).
  !> If the actual result is close to zero or one (or \(c\) and \(c+1\) respectively) with respect to a tolerance that
  !> can be specified, the routine returns zero as fractional part.
  elemental function fractional_part_scalar(x_in, c_in, tol) result(frac)
    !> Number to map to \([c, c + 1)\)
    real(dp), intent(in) :: x_in
    !> Offset of the intervall
    integer, intent(in), optional :: c_in
    !> Tolerance for real numbers to be defined as zero.
    real(dp), intent(in), optional :: tol

    real(dp) :: frac

    real(dp) :: tol_, c
    real(dp), parameter :: default_c = 0._dp

    c = default_c
    if(present(c_in)) c = real(c_in, dp)

    tol_ = default_tol
    if(present(tol)) tol_ = tol

    frac = x_in - real(floor(x_in - c), dp)
    if(abs(frac - c) < tol_) frac = c
    if(abs(frac - c - 1._dp) < tol_) frac = c
  end function fractional_part_scalar


  !> Calculate the fractional part of a matrix. See [[fractional_part_scalar(function)]].
  pure function fractional_part_matrix(X_in, C_in, tol) result(X_out)
    !> Matrix to map to the intervall
    real(dp), intent(in) :: X_in(:, :)
    !> Offset. It is expected that it has at least the same dimension as rank 1 of [[X]].
    !> Only the elements 1 to the dimension of rank 1 of [[X]] are taken into account. The \(i\)'th element of
    !> this array is taken as the offset for the \(i\)'th row of [[X]].
    integer, intent(in), optional :: C_in(:)
    !> Tolerance for real numbers to be defined as zero.
    real(dp), intent(in), optional :: tol

    real(dp), allocatable :: X_out(:, :)

    real(dp) :: tol_
    integer, parameter :: default_c = 0
    real(dp), allocatable :: X(:)
    integer, allocatable :: C(:)
    integer :: M, N

    tol_ = default_tol
    if(present(tol)) tol_ = tol

    M = size(X_in, 1)
    N = size(X_in, 2)

    if (present(C_in)) then
      if (size(C_in) /= M) error stop 'C_in has not the same dimension as rank 1 of X_in.'
      C = reshape(spread(C_in, 2, N), [M * N])
    else
      C = spread(default_c, 1, M * N)
    end if

    X = reshape(X_in, [M * N])
    X_out = reshape(fractional_part_scalar(X, C, tol_), [M, N])
  end function fractional_part_matrix


  !> Calculate the integer part of of a real number \( x \):
  !> \[
  !>   \lfloor x \rfloor.
  !> \]
  !> The routine allows for defining an offset \( c \). Then it returns
  !> \[
  !>    \lfloor x - c \rfloor.
  !> \]
  !> If the result is close to one, with respect to a tolerance that can be defined,
  !> the result as given by the intrinsic `floor` function is increased by 1.
  elemental function integer_part_scalar(x_in, c_in, tol) result(int)
    !> Number to map to \([c, c + 1)\)
    real(dp), intent(in) :: x_in
    !> Offset of the intervall
    integer, intent(in), optional :: c_in
    !> Tolerance for real numbers to be defined as zero. For default see [[default_tol]].
    real(dp), intent(in), optional :: tol

    integer :: int

    real(dp) :: tol_, c, frac
    real(dp), parameter :: default_c = 0._dp

    c = default_c
    if(present(c_in)) c = real(c_in, dp)

    tol_ = default_tol
    if(present(tol)) tol_ = tol

    int = floor(x_in - c)
    frac = x_in - c - real(int, dp)
    if(abs(frac - 1._dp) < tol_) int = int + 1
  end function integer_part_scalar


  !> Calculate the integer part of a matrix. See [[integer_part_scalar(function)]].
  pure function integer_part_matrix(X_in, C_in, tol) result(I_out)
    !> Matrix to map to the intervall
    real(dp), intent(in) :: X_in(:, :)
    !> Offset. It is expected that it has at least the same dimension as rank 1 of [[X]].
    !> Only the elements 1 to the dimension of rank 1 of [[X]] are taken into account. The \(i\)'th element of
    !> this array is taken as the offset for the \(i\)'th row of [[X]].
    integer, intent(in), optional :: C_in(:)
    !> Tolerance for real numbers to be defined as zero.
    real(dp), intent(in), optional :: tol

    integer, allocatable :: I_out(:, :)

    real(dp) :: tol_
    integer, parameter :: default_c = 0
    real(dp), allocatable :: X(:)
    integer, allocatable :: C(:)
    integer :: M, N

    tol_ = default_tol
    if(present(tol)) tol_ = tol

    M = size(X_in, 1)
    N = size(X_in, 2)

    if (present(C_in)) then
      if (size(C_in) /= M) error stop 'C_in has not the same dimension as rank 1 of X_in.'
      C = reshape(spread(C_in, 2, N), [M * N])
    else
      C = spread(default_c, 1, M * N)
    end if

    X = reshape(X_in, [M * N])
    I_out = reshape(integer_part_scalar(X, C, tol_), [M, N])
  end function integer_part_matrix

  !> Get the spherical harmonics expansion of a plane wave.
  !>
  !> \[ {\rm e}^{{\rm i} {\bf p} \cdot {\bf r}}
  !>    = 4\pi \sum_{l,m} {\rm i}^l\, j_l(pr)\, Y_{lm}^\ast(\hat{\bf p})\, Y_{lm}(\hat{\bf r}) \]
  subroutine plane_wave_in_spherical_harmonics( p, radial_grid, lmax, plane_wave_sh )
    !> wavevector in Cartesian coordinates
    real(dp), intent(in) :: p(3)
    !> radial grid
    real(dp), intent(in) :: radial_grid(:)
    !> maximum angular momentum \(l\) in expansion
    integer, intent(in) :: lmax
    !> radial functions of spherical harmonics expansion
    !> (\((l,m)\) index in first dimension, radial grid point in second)
    complex(dp), allocatable, intent(out) :: plane_wave_sh(:,:)
  
    integer :: nr, lmmax, l, m, lm, ir
    real(dp) :: p_length, p_angles(2)
    complex(dp) :: fourpi_il
  
    real(dp), allocatable :: besselj(:)
    complex(dp), allocatable :: ylm(:)
  
    nr = size( radial_grid )
    lmmax = (lmax + 1)**2
    allocate( besselj(0:lmax) )
    allocate( ylm(lmmax) )
  
    if( allocated( plane_wave_sh ) ) deallocate( plane_wave_sh )
    allocate( plane_wave_sh(lmmax, nr) )
  
    ! decompose wavevector into length and angles
    call sphcrd( p, p_length, p_angles )
    ! generate spherical harmonics of wavevector
    call genylm( lmax, p_angles, ylm )
  
    do ir = 1, nr
      ! generate spherical bessel functions
      call sbessel( lmax, p_length*radial_grid(ir), besselj )
      lm = 1
      fourpi_il = cmplx( fourpi, 0.0_dp, dp ) ! 4*pi*i^l
      ! synthesize radial functions
      do l = 0, lmax
        do m = -l, l
          plane_wave_sh(lm, ir) = fourpi_il * besselj(l) * conjg( ylm(lm) )
          lm = lm + 1
        end do
        ! multiply with i to increment exponent
        fourpi_il = cmplx( -aimag( fourpi_il ), dble( fourpi_il ), dp )
      end do
    end do
  
    deallocate( ylm, besselj )
  end subroutine plane_wave_in_spherical_harmonics

  !> For a list of non-decreasing eigenvalues \(e_{i+1} \geq e_i\) with \(i=1,\dots,N\),
  !> find blocks of degenerate eigenvalues \((l,m)_k\) with \(l \leq m\) and
  !> \(e_m - e_l < \epsilon\) for some given tolerance \(\epsilon\).
  !> The results will be returned as triples \((l,m,n)_k\), where \(l\) (\(m\))
  !> is the index of the first (last) degenerate state in the block and
  !> \(n = m-l+1\) is the degree of the degeneracy.
  pure function get_degeneracies( eval, tol ) result( deg )
    !> list of increasing eigenvalues
    real(dp), intent(in) :: eval(:)
    !> tolerance for degeneracies
    real(dp), intent(in) :: tol
    !> list of degenerate blocks
    integer, allocatable :: deg(:,:)

    real(dp), parameter :: eps0 = 1e-128_dp ! effective 0

    integer :: n, i, j
    real(dp) :: e
    integer, allocatable :: tmp(:,:)

    n = size( eval )

    do i = 2, n
      if( eval(i) + eps0 < eval(i-1) ) error stop 'Eigenvalues are not increasing.'
    end do

    if( allocated( deg ) ) deallocate( deg )
    allocate( tmp(3, n) )

    j = 0; e = -huge( 1.0_dp )
    do i = 1, n
      if( eval(i) - e < tol ) then
        tmp(2, j) = i
        tmp(3, j) = tmp(3, j) + 1
      else
        j = j + 1
        e = eval(i)
        tmp(:, j) = [i, i, 1]
      end if
    end do

    allocate( deg, source=tmp(:, 1:j) )
    deallocate( tmp )
  end function get_degeneracies

end module math_utils


