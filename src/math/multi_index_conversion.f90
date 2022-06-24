!> Module for converting between multi and composit indices where column-major order is assumed.
!> _E.g._, a 2 rank array of shape \( (4, 3) \) would be indexed as following (left is 
!> the multi index, right the composit index):
!> \[
!>    \begin{matrix}
!>       (1, 1) \rightarrow  1 & (1, 2) \rightarrow  5 & (1, 3) \rightarrow  9 \\\
!>       (2, 1) \rightarrow  2 & (2, 2) \rightarrow  6 & (2, 3) \rightarrow 10 \\\
!>       (3, 1) \rightarrow  3 & (3, 2) \rightarrow  7 & (3, 3) \rightarrow 11 \\\
!>       (4, 1) \rightarrow  4 & (4, 2) \rightarrow  8 & (4, 3) \rightarrow 12 
!>    \end{matrix}
!> \]
!>
!> The routines support fortrn 1-indexing.
module multi_index_conversion
  use asserts, only: assert
  use math_utils, only: mod1

  private
  public :: indices_to_composite_index, composite_index_to_indices
  
  !> Calculate for a given index \( i \) and the given shape of an array \( (N_1, \cdots, N_k) \) the corresponding
  !> multi index \( (i_1, \cdots, i_k) \) with \( 1 \le i_j \le N_j \).
  !>
  !> This routine can be used to unnest nested loops, _e.g._, a triple nested loop:
  !> 
  !><div class="hl"><pre>
  !><a></a></a><span class="k">do </span><span class="n">i = 1, N_1</span>
  !><a></a><span class="k">  do </span><span class="n">j = 1, N_2</span>
  !><a></a><span class="k">    do </span><span class="n">k = 1, N_3</span>
  !><a></a><span class="c">      ! Do something with i, j, k</span>
  !><a></a><span class="k">    end do </span>
  !><a></a><span class="k">  end do </span>
  !><a></a><span class="k">end do </span>
  !></pre></div>
  !> 
  !> can be rewritten as
  !>
  !><div class="hl"><pre>
  !><a></a><span class="k">do </span><span class="n">indices_to_composite_index = 1, N_1 * N_2 * N_3</span>
  !><a></a><span class="k">  call </span><span class="n">composite_index_to_indices(indices_to_composite_index, [N_1, N_2, N_3], i, j, k)</span>
  !><a></a><span class="c">  ! Do something with i, j, k </span>
  !><a></a><span class="k">end do </span>
  !></pre></div>
  !>
  !> The calculation of the indices `i, j, k` by calling `composite_index_to_indices`, which is necessary for unnesting the loop, causes a certain
  !> overhead which grows linearly with the rank of the array or the depth of the nested loop respectively.
  !> Therefore, use this routine only for unesting loops, if the overhead due to the `composite_index_to_indices` call is negiglible or you want to
  !> parallelize the nested loop.
  interface composite_index_to_indices
    module procedure :: calculate_inner_index, &
                        composite_index_to_multi_index_k, &
                        composite_index_to_multi_index, &
                        composite_index_to_double_index, &
                        composite_index_to_triple_index, &
                        composite_index_to_quartuple_index, &
                        composite_index_to_quintuple_index
  end interface composite_index_to_indices

contains

  !> Calculate for a multi index \( (i_1, \cdots, i_k) \), pointing to an element in a multi dimensional array
  !> of a given shape \( (N_1, \cdots, N_k) \), saved in column-major order,
  !> the corresponding single index \( i \):
  !> \[
  !>     (i_1, \cdots, i_k) \rightarrow i, \text{ } 1 \le i_j \le N_j,
  !> \]
  !> where the single index is calculated as
  !> \[
  !>     i = \sum_{j=k}^{2} (i_{j} - 1) \cdot \prod_{l = 1}^{j-1} N_l + i_1.
  !> \]
  integer function indices_to_composite_index(composite_index, N)
    !> Multi index \( (i_1, \cdots, i_k) \)
    integer, intent(in), contiguous :: composite_index(:)
    !> Shape of the corresponding array \( (N_1, \cdots, N_k) \)
    integer, intent(in), contiguous :: N(:)

    integer :: j

    call assert(size(composite_index) == size(N), 'composite_index and N have not the same size.')
    call assert(all(N > 0), 'Some elements of N are not greater than 0.')
    call assert(all(composite_index <= N), 'Some elements of composite_index are greater than the corresponding elements of N.')
    call assert(all(composite_index > 0), 'Some elements of composite_index are greater than the corresponding elements of N.')

    indices_to_composite_index = sum([( (composite_index(j) - 1) * product(N(: j-1)), j=size(composite_index), 2, -1 )]) + composite_index(1)
  end function indices_to_composite_index

  !> Calculate the inner and outer index from a composit index, _e.g._ from unrolling a nested loop. In other words,
  !> calculate \(i_1\) from a composit index \(i\) given by (\(i_j \in [1, \cdots, N_j]\):
  !> \[
  !>   \begin{split}
  !>     i &= \sum_{j=k}^{2} (i_{j} - 1) \cdot \prod_{l = 1}^{j-1} N_l + i_1 \\\
  !>       &= i_1 + N_1 \cdot \Biggl[(i_2 - 1)
  !>              + N_2 \cdot \biggl[(i_3 - 1)
  !>              + N_3 \cdot  \Bigl[(i_4 - 1)
  !>              + N_4 \cdot  \bigl[(i_5 - 1) + \cdots \bigl] \Bigl] \biggl]  \Biggl].
  !>   \end{split}
  !> \]
  !> Thus (see [[mod1]]):
  !> \[
  !>     i_1 = mod1(i, N_1),
  !> \]
  !> and the outer index
  !> \[
  !>     i_2 \rightarrow \frac{i - i_1} {N_1} + 1.
  !> \]
  subroutine calculate_inner_index(composite_index, inner_limit, i_1, i_2)
    !> composit index \(i\). It is updated to \( \frac{i - i_1} {N_1} + 1 \).
    integer, intent(in) :: composite_index
    !> Limit of the inner index \(N_1\).
    integer, intent(in) :: inner_limit
    !> Inner index
    integer, intent(out) :: i_1
    !> Outer index
    integer, intent(out) :: i_2

    i_1 = mod1(composite_index, inner_limit)
    i_2 = (composite_index - i_1) / inner_limit + 1
  end subroutine calculate_inner_index

  !> Calculate for a composit index \( i \), composed of \(k\) indices \( i_1, \cdots, i_k \) whith limits
  !> \( N_1, \cdots, N_k \), the indices \( i_1, \cdots, i_k \).
  subroutine composite_index_to_multi_index_k(index, k, N, composite_index)
    !> Single index \( i \)
    integer, intent(in) :: index
    !> Rank of the multi index (rank of the corresponding array).
    integer, intent(in) :: k
    !> Shape of the corresponding array \( (N_1, \cdots, N_k) \)
    integer, intent(in) :: N(k)
    !> Multi index \( (i_1, \cdots, i_k) \)
    integer, intent(out) :: composite_index(k)

    integer :: j, m, m_updated

    call assert(index > 0, 'index <= 0.')
    call assert(index <= product(N), 'composite_index exceeds the maximum possible element as determined by the product of loop limits.')
    call assert(all(N > 0), 'Some N <= 0.')

    m = index
    do j = 1, k
      call calculate_inner_index(m, N(j), composite_index(j), m_updated)
      m = m_updated
    end do
  end subroutine composite_index_to_multi_index_k

  !> Calculate for a composit index \( i \), composed of \(k\) indices \( i_1, \cdots, i_k \) whith limits
  !> \( N_1, \cdots, N_k \), the indices \( i_1, \cdots, i_k \).
  subroutine composite_index_to_multi_index(index, N, composite_index)
    !> Single index \( i \)
    integer, intent(in) :: index
    !> Shape of the corresponding array \( N_1, \cdots, N_k \)
    integer, intent(in), contiguous :: N(:)
    !> Multi index \( i_1, \cdots, i_k \)
    integer, intent(out), contiguous :: composite_index(:)

    integer :: k, j, m, m_updated

    call assert(size(composite_index) == size(N), 'composite_index and N do not have the same size.')
    call assert(index > 0, 'index <= 0.')
    call assert(index <= product(N), 'composite_index exceeds the maximum possible element as determined by the product of loop limits.')
    call assert(all(N > 0), 'Some N <= 0.')

    k = size(composite_index)
    m = index
    do j = 1, k
      call calculate_inner_index(m, N(j), composite_index(j), m_updated)
      m = m_updated
    end do
  end subroutine composite_index_to_multi_index

  !> Calculate for a composit index \( i \) , composed of  \( i_1, i_2 \) whith limits
  !> \( N_1, N_2 \), the indices \( i_1, i_2 \).
  subroutine composite_index_to_double_index(index, N, i_1, i_2)
    !> Single index \( i \)
    integer, intent(in) :: index
    !> Shape of the corresponding array \( (N_1, N_2) \)
    integer, intent(in) :: N(2)
    !> Elements of the double index
    integer, intent(out) :: i_1, i_2

    call assert(index > 0, 'index <= 0.')
    call assert(index <= product(N), 'composite_index exceeds the maximum possible element as determined by the product of loop limits.')
    call assert(all(N > 0), 'Some N <= 0.')

    call calculate_inner_index(index, N(1), i_1, i_2)
  end subroutine composite_index_to_double_index

  !> Calculate for a composit index \( i \) , composed of  \( i_1, i_2, i_3 \) whith limits
  !> \( N_1, N_2, N_3 \), the indices \( i_1, i_2, i_3 \).
  subroutine composite_index_to_triple_index(index, N, i_1, i_2, i_3)
    !> Single index \( i \)
    integer, intent(in) :: index
    !> Shape of the corresponding array \( (N_1, N_2, N_3) \)
    integer, intent(in) :: N(3)
    !> Elements of the triple index
    integer, intent(out) :: i_1, i_2, i_3

    integer :: m

    call assert(index > 0, 'index <= 0.')
    call assert(index <= product(N), 'composite_index exceeds the maximum possible element as determined by the product of loop limits.')
    call assert(all(N > 0), 'Some N <= 0.')

    m = index
    call calculate_inner_index(m, N(1), i_1, i_2)
    m = i_2
    call calculate_inner_index(m, N(2), i_2, i_3)
  end subroutine composite_index_to_triple_index

  !> Calculate for a composit index \( i \) , composed of  \( i_1, i_2, i_3, i_4 \) whith limits
  !> \( N_1, N_2, N_3, N_4 \), the indices \( i_1, i_2, i_3, i_4 \).
  subroutine composite_index_to_quartuple_index(index, N, i_1, i_2, i_3, i_4)
    !> Single index \( i \)
    integer, intent(in) :: index
    !> Shape of the corresponding array \( (N_1, N_2, N_3, N_4) \)
    integer, intent(in) :: N(4)
    !> Elements of the quartuple index
    integer, intent(out) :: i_1, i_2, i_3, i_4

    integer :: m

    call assert(index > 0, 'index <= 0.')
    call assert(index <= product(N), 'composite_index exceeds the maximum possible element as determined by the product of loop limits.')
    call assert(all(N > 0), 'Some N <= 0.')

    m = index
    call calculate_inner_index(m, N(1), i_1, i_2)
    m = i_2
    call calculate_inner_index(m, N(2), i_2, i_3)
    m = i_3
    call calculate_inner_index(m, N(3), i_3, i_4)
  end subroutine composite_index_to_quartuple_index

  !> Calculate for a composit index \( i \) , composed of  \( i_1, i_2, i_3, i_4, i_5 \) whith limits
  !> \( N_1, N_2, N_3, N_4, N_5 \), the indices \( i_1, i_2, i_3, i_4, i_5 \).
  subroutine composite_index_to_quintuple_index(index, N, i_1, i_2, i_3, i_4, i_5)
    !> Single index \( i \)
    integer, intent(in) :: index
    !> Shape of the corresponding array \( (N_1, N_2, N_3, N_4, N_5) \)
    integer, intent(in) :: N(5)
    !> Elements of the quintuple index
    integer, intent(out) :: i_1, i_2, i_3, i_4, i_5

    integer :: m

    call assert(index > 0, 'index <= 0.')
    call assert(index <= product(N), 'composite_index exceeds the maximum possible element as determined by the product of loop limits.')
    call assert(all(N > 0), 'Some N <= 0.')

    m = index
    call calculate_inner_index(m, N(1), i_1, i_2)
    m = i_2
    call calculate_inner_index(m, N(2), i_2, i_3)
    m = i_3
    call calculate_inner_index(m, N(3), i_3, i_4)
    m = i_4
    call calculate_inner_index(m, N(4), i_4, i_5)
  end subroutine composite_index_to_quintuple_index
end module multi_index_conversion