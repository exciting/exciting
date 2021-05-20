!> Math utilities and functions
module math_utils
  use iso_fortran_env, only: error_unit
  use, intrinsic :: ISO_C_BINDING

  use precision, only: sp, dp
  use constants, only: pi, zi
  use asserts, only: assert

  implicit none

  private
  public :: all_close, &
            all_zero, &
            diag, &
            kronecker_product, &
            determinant, &
            permanent, &
            mod1, &
            shuffle_vector, &
            mask_vector


  !>  Return the diagonal of a 2D array 
  Interface diag
    Module Procedure diag_int_sp, diag_real_dp, diag_complex_dp
  End Interface diag

  !>  Check if two arrays are close, within an absolute tolerance 
  Interface all_close
    Module Procedure all_close_rank0_real_dp, all_close_rank1_real_dp, &
                   & all_close_rank2_real_dp, all_close_rank0_complex_dp,&
                   & all_close_rank1_complex_dp, all_close_rank2_complex_dp
  End Interface all_close

  !>  Check if an array is zero, to within an absolute tolerance 
  Interface all_zero
  Module Procedure all_zero_rank0_real_dp, all_zero_rank1_real_dp, &
                 & all_zero_rank2_real_dp, all_zero_rank0_complex_dp,&
                 & all_zero_rank1_complex_dp, all_zero_rank2_complex_dp
  End Interface all_zero


  !> Calculate the Kronecker product for three vectors
  !>          
  !> The Kronecker product for two vectors is defined as:
  !> [\
  !>    \boldsymbol{a} \bigotimes \boldsymbol{b} =
  !>      \begin{pmatrix} 
  !>         a_1\boldsymbol{b}\\ \vdots \\ a_n\boldsymbol{b}
  !>      \end{pmatrix}
  !> \]
  interface kronecker_product
     module procedure kronecker_product_real_dp, kronecker_product_complex_dp
  end interface kronecker_product
 

contains

  !> Extracts the diagonal of an integer NxN matrix as vector of size N
  function diag_int_sp(a) result(diagonal)
   
    integer(sp), intent(in) :: a(:, :)
    integer(sp), allocatable :: diagonal(:)
   
    integer(sp) :: i
   
   call assert(size(a, dim=1) == size(a, dim=2), &
     & 'diag: Input needs to be a square matrix.')
   
   allocate (diagonal(size(a, dim=1)))
   
   do i = 1, size(diagonal)
     diagonal(i) = a(i, i)
   end do
   
  end function diag_int_sp


  !> Extracts the diagonal of a real NxN matrix as vector of size N
  function diag_real_dp(a) result(diagonal)

    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable :: diagonal(:)

    integer(sp) :: i

    call assert(size(a, dim=1) == size(a, dim=2), &
      & 'diag: Input needs to be a square matrix.')

    allocate (diagonal(size(a, dim=1)))

    do i = 1, size(diagonal)
      diagonal(i) = a(i, i)
    end do

  end function diag_real_dp


  !> Extracts the diagonal of a complex NxN matrix as vector of size N
  function diag_complex_dp(a) result(diagonal)

    complex(dp), intent(in) :: a(:, :)
    complex(dp), allocatable :: diagonal(:)
    integer(sp) :: i

    call assert(size(a, dim=1) == size(a, dim=2), &
      & 'diag: Input needs to be a square matrix.')

    allocate (diagonal(size(a, dim=1)))

    do i = 1, size(diagonal)
      diagonal(i) = a(i, i)
    end do

  end function diag_complex_dp


  !> Check if two real scalars \( a \) and \( b \) are equal,
  !> where equal is defined as \( |a - b| \leq abs\_tol \).
  logical function all_close_rank0_real_dp(a, b, abs_tolerance)

    !> Input array 
    real(dp), intent(in) :: a
    !> Reference array
    real(dp), intent(in) :: b
    !> Absolute tolerance for input and reference to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if

    all_close_rank0_real_dp = abs(a - b) <= abs_tol

  end function all_close_rank0_real_dp

  !> Check if two real rank-1 arrays \( \mathbf{a} \) and \( \mathbf{b} \)
  !> are equal, where equal is defined as
  !> \( |a_i - b_i| \leq abs\_tol,  \forall i \).
  !> As such, the tolerance is checked elementwise.
  logical function all_close_rank1_real_dp(a, b, abs_tolerance)

    !> Input array 
    real(dp), intent(in) :: a(:)
    !> Reference array
    real(dp), intent(in) :: b(:)
    !> Absolute tolerance for input and reference to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    call assert(size(a) == size(b), &
      & 'all_close_rank1_real_dp: size of input arrays differs.')

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if


    all_close_rank1_real_dp = all(abs(a - b) <= abs_tol)

  end function all_close_rank1_real_dp

  
  !> Check if two real rank-2 arrays \( \mathbf{a} \) and \( \mathbf{b} \)
  !> are equal, where equal is defined as
  !> \( |a_{ij} - b_{ij}| \leq abs\_tol,  \forall i,j \).
  !> As such, the tolerance is checked elementwise.
  logical function all_close_rank2_real_dp(a, b, abs_tolerance)

    !> Input array 
    real(dp), intent(in) :: a(:,:)
    !> Reference array
    real(dp), intent(in) :: b(:,:)
    !> Absolute tolerance for input and reference to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    call assert(size(a) == size(b), &
      & 'all_close_rank2_real_dp: size of input arrays differs.')

    call assert(all(shape(a) == shape(b)), &
      & 'all_close_rank2_real_dp: shape of input arrays differs.')

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if

    all_close_rank2_real_dp = all(abs(a - b) <= abs_tol)

  end function all_close_rank2_real_dp


  !> Check if two complex scalars \( a \) and \( b \)
  !> are equal, where equal is defined as \( |a - b| \leq abs\_tol \).
  !> The abs of a complex scalar \( a \) is defined as:
  !> [\
  !>    |a_i| = \sqrt{ a_i \cdot {a_i}^*}.
  !> ]\
  !> As such, the tolerance is a real value.
  logical function all_close_rank0_complex_dp(a, b, abs_tolerance)

    !> Input array 
    complex(dp), intent(in) :: a
    !> Reference array
    complex(dp), intent(in) :: b
    !> Absolute tolerance for input and reference to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if

    all_close_rank0_complex_dp = abs(a - b) <= abs_tol

  end function all_close_rank0_complex_dp

  !> Check if two complex rank-1 arrays \( \mathbf{a} \) and \( \mathbf{b} \)
  !>  are equal, where equal is defined as
  !>  \( |a_i - b_i| \leq abs\_tol,  \forall i \).
  !> The abs of a complex scalar \( a_i \) is defined as:
  !> [\
  !>    |a_i| = \sqrt{ a_i \cdot {a_i}^*}.
  !> ]\
  !> As such, the tolerance is a real value and is 
  !> checked elementwise.
  logical function all_close_rank1_complex_dp(a, b, abs_tolerance)

    !> Input array 
    complex(dp), intent(in) :: a(:)
    !> Reference array
    complex(dp), intent(in) :: b(:)
    !> Absolute tolerance for input and reference to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    call assert(size(a) == size(b), &
      & 'all_close_rank1_complex_dp: size of input arrays differs.')

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if

    all_close_rank1_complex_dp = all(abs(a - b) <= abs_tol)

  end function all_close_rank1_complex_dp


  !> Check if two rank-2 arrays \( \mathbf{a} \) and \( \mathbf{b} \)
  !> of are equal, where equal is defined as
  !>  \( |a_{ij} - b_{ij}| \leq abs\_tol,  \forall i,j \).
  !> The abs of a complex scalar \( a_i \) is defined as:
  !> [\
  !>    |a_i| = \sqrt{ a_i \cdot {a_i}^*}.
  !> ]\
  !> As such, the tolerance is a real value and is 
  !> checked elementwise.
  logical function all_close_rank2_complex_dp(a, b, abs_tolerance)

    !> Input array 
    complex(dp), intent(in) :: a(:,:)
    !> Reference array
    complex(dp), intent(in) :: b(:,:)
    !> Absolute tolerance for input and reference to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    call assert(size(a) == size(b), &
      & 'all_close_rank2_complex_dp: size of input arrays differs.')

    call assert(all(shape(a) == shape(b)), &
      & 'all_close_rank2_complex_dp: shape of input arrays differs.')

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if

    all_close_rank2_complex_dp = all(abs(a - b) <= abs_tol)

  end function all_close_rank2_complex_dp


  !> Check if a complex scalar \( a \) is zero, 
  !> where zero is defined as \( |a - b| \leq abs\_tol \).
  !> The abs of a complex scalar \( a \) is defined as:
  !> [\
  !>    |a| = \sqrt{ a \cdot {a}^*}.
  !> ]\
  !> As such, the tolerance is a real value.
  logical function all_zero_rank0_complex_dp(a, abs_tolerance)

    !> Input array 
    complex(dp), intent(in) :: a
    !> Absolute tolerance for input and reference to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if

    all_zero_rank0_complex_dp = abs(a ) <= abs_tol

  end function all_zero_rank0_complex_dp

  !> Check if a complex rank-1 array \( \mathbf{a} \) is zero,
  !>  where zero is defined as \( |a_i| \leq abs\_tol,  \forall i \).
  !> The abs of a complex scalar \( a_i \) is defined as:
  !> [\
  !>    |a_i| = \sqrt{ a_i \cdot {a_i}^*}.
  !> ]\
  !> As such, the tolerance is a real value and is 
  !> checked elementwise.
  logical function all_zero_rank1_complex_dp(a, abs_tolerance)

    !> Input array 
    complex(dp), intent(in) :: a(:)
    !> Absolute tolerance for input and reference to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if

    all_zero_rank1_complex_dp = all(abs(a) <= abs_tol)

  end function all_zero_rank1_complex_dp


  !> Check if a complex rank-2 array \( \mathbf{a} \) is zero, 
  !> where zero  is defined as \( |a_{ij}| \leq abs\_tol,  \forall i,j \).
  !> The abs of a complex scalar \( a_i \) is defined as:
  !> [\
  !>    |a_i| = \sqrt{ a_i \cdot {a_i}^*}.
  !> ]\
  !> As such, the tolerance is a real value and is 
  !> checked elementwise.
  logical function all_zero_rank2_complex_dp(a, abs_tolerance)

    !> Input array 
    complex(dp), intent(in) :: a(:,:)
    !> Absolute tolerance for input and reference to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if

    all_zero_rank2_complex_dp = all(abs(a) <= abs_tol)

  end function all_zero_rank2_complex_dp


   
  !> Check if a real scalar \(a \) is zero, where zero 
  !> is defined as \( |a| \leq abs\_tol \).
  logical function all_zero_rank0_real_dp(a, abs_tolerance)

    !> Input array 
    real(dp), intent(in) :: a
    !> Absolute tolerance for input to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if

    all_zero_rank0_real_dp = abs(a) <= abs_tol

  end function all_zero_rank0_real_dp


  !> Check if a real rank-1 array \( \mathbf{a} \) is zero, 
  !> where zero  is defined as \( |a_i| \leq abs\_tol,  \forall i \).
  !> As such, the tolerance is checked elementwise.
  logical function all_zero_rank1_real_dp(a, abs_tolerance)

    !> Input array 
    real(dp), intent(in) :: a(:)
    !> Absolute tolerance for input to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if

    all_zero_rank1_real_dp = all(abs(a) <= abs_tol)

  end function all_zero_rank1_real_dp


  !> Check if a real rank-2 array \( \mathbf{a} \) is zero, 
  !> where zero  is defined as \( |a_{ij}| \leq abs\_tol,  \forall i,j \).
  !> As such, the tolerance and is checked elementwise.
  logical function all_zero_rank2_real_dp(a, abs_tolerance)

    !> Input array 
    real(dp), intent(in) :: a(:,:)
    !> Absolute tolerance for input to be considered equal 
    real(dp), intent(in), optional :: abs_tolerance

    !> Local absolute tolerance 
    real(dp) :: abs_tol

    if (present(abs_tolerance)) then
      abs_tol = abs_tolerance
    else
      abs_tol = 1.e-10_dp
    end if

    all_zero_rank2_real_dp = all(abs(a) <= abs_tol)

  end function all_zero_rank2_real_dp


  !> Calculate the Kronecker product for three vectors dimensions:
  !> (\ \mathbf{A} \bigotimes \mathbf{B} \bigotimes \mathbf{C} \)
  !>
  !> The Kronecker product for two vectors is defined as:
  !> [\
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


  !> Calculate the Kronecker product for three complex double vectors with 
  !> arbitrary dimensions, (\ \mathbf{A} \bigotimes \mathbf{B} \bigotimes \mathbf{C} \).
  !>
  !> The Kronecker product for two vectors is defined as:
  !> [\
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


  !> Calculates the determinant of a real matrix using Laplace extension.
  !> See [[determinant_laplace(subroutine)]] for details.
  real(dp) function determinant(A)
    !> Matrix for which the determinant is calculated
    real(dp), dimension(:, :), intent(in) :: A

    call assert(size(A, 1) == size(A, 2), &
      & message='determinant: Input needs to be a square matrix.')
    determinant = determinant_laplace(A, -1)
      
  end function determinant


  !> Calculates the permanent of a real matrix using Laplace extension.
  !> See [[determinant_laplace(subroutine)]] for details.
  real(dp) function permanent(A)
    !> Matrix for which the permanent is calculated
    real(dp), dimension(:, :), intent(in) :: A

    call assert(size(A, 1) == size(A, 2), &
      &  message='permanent: Input needs to be a square matrix.')
    permanent = determinant_laplace(A, 1)
      
  end function permanent


  !> Calculate the determinant of a real matrix using Laplace extension.
  !>
  !> [\
  !>     \det A = \sum_{i=1}^n (-1)^{i+j}a_ij \cdot \det A_{ij}
  !> \]
  !> where $A_{ij}$ is derived by canceling the $i$'th row and the $j$'th collumn of $A$.
  !>
  !> This implementation is based on the reference found at 
  !> [rossetta code](http://rosettacode.org/wiki/Determinant_and_permanent#Fortran)
  recursive function determinant_laplace(a, permanent) result(accumulation)
    !> Matrix for which the determinant is calculated
    real(dp), dimension(:, :), intent(in) :: a

    !> Setting permanent to 1 computes the permanent.
    !> Setting permanent to -1 computes the determinant.
    integer, intent(in), optional :: permanent
      
    ! local variables:
    real(dp), allocatable :: b(:,:)
    integer :: permanent_, i, sgn
    real(dp) :: accumulation
    !> Dimension of input matrix
    integer(sp) :: n

    call assert(size(a, 1) == size(a, 2), &
      &  message='permanent: Input needs to be a square matrix.')
      
    n = size(a, dim=2)

    allocate(b(n - 1, n - 1))
    permanent_ = -1
    if (present(permanent)) then
      call assert( abs(permanent) == 1, & 
             & message="determinant_laplace: permanent needs to be 1 or -1.")
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
        accumulation = accumulation + sgn * a(1, i) * determinant_laplace(b, permanent_)
        sgn = sgn * permanent_
      end do
    end if

  end function determinant_laplace

   
  !> Modulus after division, returning in the range (0,N]
  !>
  !> For example:
  !>   mod1(3,2) = 1
  !>   mod1(4,2) = 2
  integer function mod1(M, N)
    implicit none
    !> integer to translate
    integer, intent(in) :: M, &
                           !> length of the cycle
                           N

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

  end function mod1

  !> Returns an integer vector of length N with all numbers between 1 and N randomly ordered. The result can be used to shuffle a vector.
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

  !> Filter a vector with a mask vector.
  !> 
  !> The mask vector is a vector of logicals, which is used to determine
  !> which elements to retain in 'vector_in' 
  !> 
  !> For example:
  !>   [1, 3] = mask_vector([1, 2, 3, 4], [.true., .false., .true., .false.])
  !>   
  function mask_vector(vector_in, mask) result(vector_out)
    implicit none
    !> Input vector
    complex(dp), intent(in) :: vector_in(:)
    !> Mask vector
    logical, intent(in) :: mask(:)
    !> Output vector
    complex(dp), allocatable :: vector_out(:)
    ! local variables
    integer :: i, j, dim_out
        
    call assert(size(vector_in) == size(mask), &
        message = "math_utils: mask: vector_in and mask_vector differ in size.")
    
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

end module math_utils
