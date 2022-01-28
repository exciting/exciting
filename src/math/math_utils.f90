!> Math utilities and functions
module math_utils
  use iso_fortran_env, only: error_unit
  use, intrinsic :: ISO_C_BINDING

  use precision, only: sp, dp
  use constants, only: pi, zzero, zone, zi
  use asserts, only: assert

  implicit none

  private
  public :: identity_real_dp, &
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
            shuffle_vector, &
            mask_vector, &
            round_down, &
            is_positive_definite

  !> default tolerance
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
    module procedure is_orthogonal_real_dp, is_unitary_complex_dp
  end interface is_unitary

  !>  Check if two arrays are close, within an absolute tolerance
  interface all_close
    module procedure all_close_rank0_real_dp, all_close_rank1_real_dp, &
                   & all_close_rank2_real_dp, all_close_rank0_complex_dp,&
                   & all_close_rank1_complex_dp, all_close_rank2_complex_dp
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
     module procedure kronecker_product_real_dp, kronecker_product_complex_dp
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


contains


! identity_real_dp, identity_complex_dp
!
! Setup identity matrix

  !> Setup real identity matrix.
  function identity_real_dp(N) result(identity)
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
  function identity_complex_dp(N) result(identity)
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
  function is_square_real_dp(a) result(square)
    !> Input matrix
    real(dp), intent(in) :: a(:, :)
    !> Is the matrix square
    logical :: square
    square = size(a, dim=1) == size(a, dim=2)
  end function


  !> Check if a complex matrix is square
  function is_square_complex_dp(a) result(square)
    !> Input matrix
    complex(dp), intent(in) :: a(:, :)
    !> Is the matrix square
    logical :: square
    square = size(a, dim=1) == size(a, dim=2)
  end function


  !> Check if a integer matrix is square
  function is_square_integer(a) result(square)
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

  !> Check if a real matrix is orthognal
  logical function is_orthogonal_real_dp(A, tol)
    !> Matrix to check for
    real(dp), intent(in) :: A(:, :)
    !> Tolerance for defining equiv
    real(dp), intent(in), optional :: tol

    real(dp) tol_
    real(dp), allocatable :: A_A_T(:, :)

    tol_ = default_tol
    if(present(tol)) tol_ = tol

    if(size(A, dim=1) <= size(A, dim=2)) then
      A_A_T = matmul(A, transpose(A))
    else
      A_A_T = matmul(transpose(A), A)
    end if

    is_orthogonal_real_dp = all_close(A_A_T, identity_real_dp(size(A_A_T, dim=1)), tol_)
  end function is_orthogonal_real_dp


  !> Check if a complex matrix is unitary
  logical function is_unitary_complex_dp(A, tol)
    !> Matrix to check for
    complex(dp), intent(in) :: A(:, :)
    !> Tolerance for defining equiv
    real(dp), intent(in), optional :: tol

    real(dp) tol_
    real(dp), allocatable :: A_A_T(:, :)

    tol_ = default_tol
    if(present(tol)) tol_ = tol

    if(size(A, dim=1) <= size(A, dim=2)) then
      A_A_T = matmul(A, transpose(conjg(A)))
    else
      A_A_T = matmul(transpose(conjg(A)), A)
    end if

    is_unitary_complex_dp = all_close(A_A_T, identity_real_dp(size(A_A_T, dim=1)), tol_)
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
  logical function all_zero_rank0_complex_dp(a, tol)

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
  logical function all_zero_rank1_complex_dp(a, tol)

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
  logical function all_zero_rank2_complex_dp(a, tol)

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
  logical function all_zero_rank0_real_dp(a, tol)

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
  logical function all_zero_rank1_real_dp(a, tol)

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
  logical function all_zero_rank2_real_dp(a, tol)

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
    integer :: i
    real(dp), allocatable :: work(:)

    allocate (output_vector(size(A)*size(B)*size(C)))
    allocate (work(size(A)*size(B)))

    do i = 1, size(A)
      work((i - 1)*size(B) + 1:i*size(B)) = A(i)*B
    end do

    do i = 1, size(A)*size(B)
      output_vector((i - 1)*size(C) + 1:i*size(C)) = work(i)*C
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
    integer :: i
    complex(dp), allocatable :: work(:)

    allocate (output_vector(size(A)*size(B)*size(C)))
    allocate (work(size(A)*size(B)))

    do i = 1, size(A)
      work((i - 1)*size(B) + 1:i*size(B)) = A(i)*B
    end do

    do i = 1, size(A)*size(B)
      output_vector((i - 1)*size(C) + 1:i*size(C)) = work(i)*C
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

  !> Modulus after floor division, returning in the range (0,N]
  !> For example:
  !>   mod1(3,2) = 1
  !>   mod1(4,2) = 2
  integer function mod1(M, N)
    !> integer to translate
    integer, intent(in) :: M
    !> length of the cycle
    integer, intent(in) :: N

    if (M > 0) then
     if (mod(M, N) == 0) then
       mod1 = N
     else
       mod1 = mod(M, N)
     end if

    else
      if (mod(M + N, N) == 0) then
        mod1 = N
      else
        mod1 = mod(M + N, N)
      end if
    end if
  end function


! shuffle_vector

  !> Return an integer vector of length N with all numbers between 1 and N randomly ordered.
  !> The result can be used to shuffle a vector.
  function shuffle_vector(N) result (p)
    ! input/output
    !> length of the permutation
    integer, intent(in) :: N
    !> random permutation
    integer :: p(N)
    ! local variables
    integer :: j, k
    real(sp) :: u

    p = 0
    do j = 1, N
      call random_number(u)
      k = floor(j*u) + 1
      p(j) = p(k)
      p(k) = j
    end do
  end function shuffle_vector



! mask_vector

  !> Filter a vector with a mask vector.
  !>
  !> The mask vector is a vector of logicals, which is used to determine
  !> which elements to retain in 'vector_in'
  !>
  !> For example:
  !>   [1, 3] = mask_vector([1, 2, 3, 4], [.true., .false., .true., .false.])
  function mask_vector(vector_in, mask) result(vector_out)
    !> Input vector
    complex(dp), intent(in) :: vector_in(:)
    !> Mask vector
    logical, intent(in) :: mask(:)
    !> Output vector
    complex(dp), allocatable :: vector_out(:)
    ! local variables
    integer :: i, j, dim_out

    call assert(size(vector_in) == size(mask), &
        message = "mask_vector: vector_in and mask_vector differ in size.")

    dim_out = count(mask)
    allocate(vector_out(dim_out))

    j = 1
    do i = 1, size(vector_in)
      if(mask(i)) then
        vector_out(j) = vector_in(i)
        j = j + 1
      end if
    end do
  end function mask_vector

  !> !> Rounds a real number down to the nth decimal place, e.g.
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


  !> Check if a hermitian matrix is positive-definite. This is done by checking
  !> if all the eigenvalues are positive
  logical function is_positive_definite( A )
    !> Matrix to be checked
    complex(dp), intent(in)   :: A(:, :)

    integer                   :: dim, info, lwork
    real(dp), allocatable     :: rwork(:), eigenvalues(:)
    complex(dp), allocatable  :: A_copy(:, :), work(:)
    character(200)            :: error_msg

    call assert( is_hermitian(A), 'Error(is_positive_definite): A is not &
      & hermitian' )
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
    call ZHEEV( 'N', 'U', dim, A_copy, dim, eigenvalues, work, lwork, rwork, &
      & info )
    write(error_msg,*) &
      & 'Error(is_positive_definite): ZHEEV returned info = ', info
    call assert( info==0, error_msg )
    is_positive_definite = all( eigenvalues > 0._dp )
  end function

  


  

end module math_utils
