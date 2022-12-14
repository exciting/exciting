!> Module for general tensor contractions.
module tensor_contractions
    use precision, only: dp
    use constants, only: zone, zzero, real_zero, real_one
    use lapack_f95_interfaces, only: zgemm
    use general_matrix_multiplication, only: gemm_parameters
    use asserts, only: assert

    implicit none

    private
    public :: complex_tensor_contraction_dp,&
              real_tensor_contraction_dp, shape_rank_2

    !> Default for **trans_A** and **trans_B** respectively which define
    !> if the tensor is used, its transpose or its
    !> conjugate transpose in case of complex tensors. Default is to use the matrix (`'N'`).
    character(len=1), parameter :: default_trans_char = 'N'



contains

    !> Calculates the contraction of two complex tensors \( \mathbf{A} \) 
    !> and \( \mathbf{B} \) of ranks > 2  over one or more indices. 
    !> The contractions can be performed according to one of the following 
    !> possibilites (i.e., the contraction is performed over all indices
    !> from \(d_i \) to \(d_f \) which are common to  
    !> \( \mathbf{A} \) and  \( \mathbf{B} \)):
    !>
    !>      <li> 1)  \[
    !>              C_{a_1,...,a_n,b_1,...,b_n} =
    !>                \sum_{d_i...d_f} A_{a_1,...,a_n,d_i,...,d_f}
    !>                           \cdot B_{d_i,...,d_f,b_1,...,b_n}
    !>                                                                      \],
    !>      </li>
    !>
    !>      <li> 2)  \[
    !>              C_{a_1,...,a_n,b_1,...,b_n} =
    !>                \sum_{d_i...d_f} A_{a_1,...,a_n,d_i,...,d_f}
    !>                           \cdot B_{b_1,...,b_n,d_i,...,d_f}
    !>                                                                      \],
    !>      </li>
    !>
    !>      <li> 3)  \[
    !>              C_{a_1,...,a_n,b_1,...,b_n} =
    !>                \sum_{d_i...d_f} A_{d_i,...,d_f,a_1,...,a_n}
    !>                           \cdot B_{d_i,...,d_f,b_1,...,b_n}
    !>                                                                      \],
    !>      </li>
    !>
    !>      <li> 4)  \[
    !>              C_{a_1,...,a_n,b_1,...,b_n} =
    !>                \sum_{d_i...d_f} A_{d_i,...,d_f,a_1,...,a_n}
    !>                           \cdot B_{b_1,...,b_n,d_i,...,d_f}
    !>                                                                      \],
    !>      </li>
    !>
    !>      <li> 5)  \[
    !>              C_{a_1,...,a_n,b_1,...,b_n} =
    !>                \sum_{d_i...d_f} A_{a_1,...,a_n,d_i,...,d_f}
    !>                           \cdot B^*_{b_1,...,b_n,d_i,...,d_f}
    !>                                                                      \],
    !>      <li> 6)  
    !>          \[
    !>              C_{a_1,...,a_n,b_1,...,b_n} =
    !>                \sum_{d_i...d_f} A^*_{d_i,...,d_f,a_1,...,a_n}
    !>                           \cdot B_{d_i,...,d_f,b_1,...,b_n}
    !>                                                                      \], !>      </li>
    !>
    !>      <li> 7)
    !>          \[
    !>              C_{a_1,...,a_n,b_1,...,b_n} =
    !>                \sum_{d_i...d_f} A^*_{d_i,...,d_f,a_1,...,a_n}
    !>                           \cdot B^*_{b_1,...,b_n,d_i,...,d_f}
    !>                                                                      \],
    !>     </li>
    !>
    !> which correspond to the following input parameters:
    !>
    !>   <li>    1) **trans_A** = `'n'`, **trans_B** = `'n'`,   </li>
    !>   <li>    2) **trans_A** = `'n'`, **trans_B** = `'t'`,   </li>
    !>   <li>    3) **trans_A** = `'t'`, **trans_B** = `'n'`,   </li>
    !>   <li>    4) **trans_A** = `'t'`, **trans_B** = `'t'`,   </li>
    !>   <li>    5) **trans_A** = `'n'`, **trans_B** = `'c'`,   </li>
    !>   <li>    6) **trans_A** = `'c'`, **trans_B** = `'n'`,   </li>
    !>   <li>    7) **trans_A** = `'c'`, **trans_B** = `'c'`.   </li>
    !>
    !> The contractions are performed by considering tensors of rank d as 
    !> matrices of rank 2 of the form
    !>
    !> \[
    !>    A_{a_1,...,a_n,d_i,...,d_f} \rightarrow
    !>    A_{a_1\cdot...\cdot a_n,d_i \cdot ... \cdot d_f}
    !>                                                            \]
    !>
    !> such that general LAPACK routines can be used.
    !> Contractions are only possible if the contracted dimensions are
    !>
    !>    <li>  a) contigous in memory,  </li>
    !>    <li>  b) are the leading or trailing dimensions 
    !>             of the contracted tensors.               </li>
    !>
    subroutine complex_tensor_contraction_dp(A, shape_A, B, shape_B,&
        C, shape_C,  trans_A, trans_B)
        
        !> input tensor A
        complex(dp), intent(in) :: A(*)
        !> shape of tensor A
        integer, intent(in) :: shape_A(:)
        !> input tensor B
        complex(dp), intent(in) :: B(*)
        !> shape of tensor B
        integer, intent(in) :: shape_B(:)
        !> shape of tensor C
        integer, intent(in) :: shape_C(:)
        !> output tensor C
        complex(dp), intent(out) :: C(*)
        !> transposition of A ('N' (default), 'T', or 'C')
        character(len=1), optional, intent(in) :: trans_A
        !> transposition of B ('N' (default), 'T', or 'C')
        character(len=1), optional, intent(in) :: trans_B

        !> number of dimensions to be contracted
        integer:: n_contracted
        !> LAPACK variables
        integer :: M, N, K, LDA, LDB, LDC
        !> Local copy of trans_A
        character(len=1) :: trans_A_char
        !> Local copy of trans_B
        character(len=1) :: trans_B_char
        !> True if A is transposed
        logical :: A_is_transposed
        !> True if B is transposed
        logical :: B_is_transposed
        !> Two-dimensional shapes of tensors A, B and C
        integer :: shape_A_2d(2), shape_B_2d(2), shape_C_2d(2)
        !> Ranks of tensors A, B and C
        integer :: rank_A, rank_B, rank_C

        rank_A = size(shape_A)
        rank_B = size(shape_B)
        rank_C = size(shape_C)

        ! Check: rank_C = rank_A + rank_B - 2*n_contracted
        call assert(mod(rank_A + rank_B - rank_C, 2) == 0,&
                        'Ranks of arrays are not correct: &
                         mod(rank_A + rank_B - rank_C, 2) /= 0.')

        n_contracted = (rank_A + rank_B - rank_C) / 2

        trans_A_char = default_trans_char
        if (present(trans_A)) trans_A_char = trans_A

        trans_B_char = default_trans_char
        if (present(trans_B)) trans_B_char = trans_B

        A_is_transposed = any(trans_A_char == ['T', 't', 'C', 'c'])
        B_is_transposed = any(trans_B_char == ['T', 't', 'C', 'c'])

        shape_A_2d = shape_rank_2(shape_A, n_contracted, .not. A_is_transposed)
        shape_B_2d = shape_rank_2(shape_B, n_contracted, B_is_transposed)
   
        call gemm_parameters(shape_A_2d, shape_B_2d, trans_A_char, trans_B_char, &
                          M, N, K, LDA, LDB, LDC)

        shape_C_2d = [M, N]

        call zgemm(trans_A_char, trans_B_char, M, N, K, zone, &
                   A, LDA, B, LDB, zzero, C, LDC)

    end subroutine complex_tensor_contraction_dp


    !> Calculates the contraction of two real tensors \( \mathbf{A} \) 
    !> and \( \mathbf{B} \) of ranks > 2  over one or more indices. 
    !> The contractions can be performed according to one of the following 
    !> possibilites (i.e., the contraction is performed over all indices
    !> from \(d_i \) to \(d_f \) which are common to A and B):
    !>
    !>      <li> 1)  \[
    !>              C_{a_1,...,a_n,b_1,...,b_n} =
    !>                \sum_{d_i...d_f} A_{a_1,...,a_n,d_i,...,d_f}
    !>                           \cdot B_{d_i,...,d_f,b_1,...,b_n}
    !>                                                                      \],
    !>      </li>
    !>
    !>      <li> 2)  \[
    !>              C_{a_1,...,a_n,b_1,...,b_n} =
    !>                \sum_{d_i...d_f} A_{a_1,...,a_n,d_i,...,d_f}
    !>                           \cdot B_{b_1,...,b_n,d_i,...,d_f}
    !>                                                                      \],
    !>      </li>
    !>
    !>      <li> 3)  \[
    !>              C_{a_1,...,a_n,b_1,...,b_n} =
    !>                \sum_{d_i...d_f} A_{d_i,...,d_f,a_1,...,a_n}
    !>                           \cdot B_{d_i,...,d_f,b_1,...,b_n}
    !>                                                                      \],
    !>      </li>
    !>
    !>      <li> 4)  \[
    !>              C_{a_1,...,a_n,b_1,...,b_n} =
    !>                \sum_{d_i...d_f} A_{d_i,...,d_f,a_1,...,a_n}
    !>                           \cdot B_{b_1,...,b_n,d_i,...,d_f}
    !>                                                                      \],
    !>      </li>
    !>
    !> which correspond to the following input parameters:
    !>
    !>   <li>    1) **trans_A** = `'n'`, **trans_B** = `'n'`,   </li>
    !>   <li>    2) **trans_A** = `'n'`, **trans_B** = `'t'`,   </li>
    !>   <li>    3) **trans_A** = `'t'`, **trans_B** = `'n'`,   </li>
    !>   <li>    4) **trans_A** = `'t'`, **trans_B** = `'t'`,   </li>
    !>
    !> The contractions are performed by considering tensors of rank d as 
    !> matrices of rank 2 of the form
    !>
    !> \[
    !>    A_{a_1,...,a_n,d_i,...,d_f} \rightarrow
    !>    A_{a_1\cdot...\cdot a_n,d_i \cdot ... \cdot d_f}
    !>                                                            \]
    !>
    !> such that general LAPACK routines can be used.
    !> Contractions are only possible if the contracted dimensions are
    !>
    !>    <li>  a) contigous in memory,  </li>
    !>    <li>  b) are the leading or trailing dimensions 
    !>             of the contracted tensors.               </li>
    !>
    subroutine real_tensor_contraction_dp(A, shape_A, B, shape_B,&
        C, shape_C,  trans_A, trans_B)
        
        !> input tensor A
        real(dp), intent(in) :: A(*)
        !> shape of tensor A
        integer, intent(in) :: shape_A(:)
        !> input tensor B
        real(dp), intent(in) :: B(*)
        !> shape of tensor B
        integer, intent(in) :: shape_B(:)
        !> shape of tensor C
        integer, intent(in) :: shape_C(:)
        !> output tensor C
        real(dp), intent(out) :: C(*)
        !> transposition of A ('N' (default), 'T', or 'C')
        character(len=1), optional, intent(in) :: trans_A
        !> transposition of B ('N' (default), 'T', or 'C')
        character(len=1), optional, intent(in) :: trans_B

        !> number of dimensions to be contracted
        integer:: n_contracted
        !> LAPACK variables
        integer :: M, N, K, LDA, LDB, LDC
        !> Local copy of trans_A
        character(len=1) :: trans_A_char
        !> Local copy of trans_B
        character(len=1) :: trans_B_char
        !> True if A is transposed
        logical :: A_is_transposed
        !> True if B is transposed
        logical :: B_is_transposed
        !> Two-dimensional shapes of tensors A, B and C
        integer :: shape_A_2d(2), shape_B_2d(2), shape_C_2d(2)
        !> Ranks of tensors A, B and C
        integer :: rank_A, rank_B, rank_C


        rank_A = size(shape_A)
        rank_B = size(shape_B)
        rank_C = size(shape_C)

        ! Check: rank_C = rank_A + rank_B - 2*n_contracted
        call assert(mod(rank_A + rank_B - rank_C, 2) == 0,&
                        'Ranks of arrays are not correct.')

        n_contracted = (rank_A + rank_B - rank_C) / 2

        trans_A_char = default_trans_char
        if (present(trans_A)) trans_A_char = trans_A

        trans_B_char = default_trans_char
        if (present(trans_B)) trans_B_char = trans_B

        A_is_transposed = any(trans_A_char == ['T', 't'])
        B_is_transposed = any(trans_B_char == ['T', 't'])

        shape_A_2d = shape_rank_2(shape_A, n_contracted, .not. A_is_transposed)
        shape_B_2d = shape_rank_2(shape_B, n_contracted, B_is_transposed)

        call gemm_parameters(shape_A_2d, shape_B_2d, trans_A_char, trans_B_char, &
                          M, N, K, LDA, LDB, LDC)

        shape_C_2d = [M, N]

        call dgemm(trans_A_char, trans_B_char, M, N, K, real_one, &
                   A, LDA, B, LDB, real_zero, C, LDC)

    end subroutine real_tensor_contraction_dp


    !> Returns the 2-dimensional shape of a complex array \( \mathbf{A} \)
    !> (which is of \( rank(\mathbf{A}) > 2) such that a tensor contraction  
    !> using the routines defined in this module can be performed.
    !> Given the number \( n_c \) of dimensions over which the contraction 
    !> is performed and the 'location' of these dimensions, i.e. if they
    !> are the leading or trailing dimensions, the 2-dimensional shape is 
    !> defined in the following way:
    !> 
    !> If the  contracted dimensions are the trailing dimensions,
    !> the  lower and upper limits of the combined leading dimensions
    !>  are defined by  \(  [1, rank(\mathbf{A}) - n_c] \)
    !>  and the lower and upper limits of the combined  trailing dimensions by 
    !> \( [rank(\mathbf{A}) - n_c + 1, rank(\mathbf{A})] \).
    !>
    !> If  the contracted dimensions are the leading dimensions,
    !> the  lower and upper limits of the combined leading dimensions
    !> are defined by \(  [1, n_c] \) and the lower and upper limits
    !> of the combined trailing dimensions by \( [n_c + 1, rank(\mathbf{A})] \).
    !> 
    !> Example: 
    !> The 2d shape of the 3 x 3 x 4  tensor \( \mathbf{A} \) 
    !> for a contraction over two indices is given by 
    !> (9, 4) if the contracted dimensions are the leading dimensions and by
    !> (3, 12) if the contracted dimensions are the trailing dimensions.
    function shape_rank_2(shape_A, n_contracted_dims, contracted_dims_are_trailing) result(shape_A_2d)

        use modmpi, only: terminate_mpi_env, mpiglobal

        !> Shape of A
        integer, intent(in) :: shape_A(:)
        !> Number of dimensions over which the contraction is performed
        integer, intent(in) :: n_contracted_dims
        !> True if the contracted dimensions are the trailing dimensions
        logical, intent(in) :: contracted_dims_are_trailing
        
        !> Rank of A
        integer :: rank_A
        !> 2d-shape A
        integer :: shape_A_2d(2)
        !> Number of leading dimension to be combined to a single dimension
        integer :: n_combined_leading
        !> Number of trailing dimension to be combined to a single dimension
        integer :: n_combined_trailing


        rank_A = size(shape_A)

        call assert(rank_A >= 2, 'rank_A < 2') 

       if (contracted_dims_are_trailing) then
            n_combined_leading = product(shape_A(1:rank_A - n_contracted_dims))
            n_combined_trailing = product(shape_A(rank_A - n_contracted_dims +1:rank_A))
        else   
            n_combined_leading =  product(shape_A(1:n_contracted_dims))
            n_combined_trailing = product(shape_A(n_contracted_dims +1:rank_A))
        end if
             
        shape_A_2d = [n_combined_leading, n_combined_trailing]

    end function

end module
