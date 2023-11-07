! Created by  on 22/09/2022.

module dynamic_indices
  use asserts, only: assert
  use iso_c_binding, only: c_ptr

  implicit none

  private
  public :: dynindex_type

  type dynindex_type
    integer, allocatable :: unrolled_indices(:)
    integer, allocatable :: n_per_col(:)
    integer :: n_cols
    integer :: n_indices

    contains

    procedure :: init_dynint, init_dynint_from_matrix
    generic :: init => init_dynint, init_dynint_from_matrix

    procedure :: finit => finit_dynint

    procedure :: col => get_col

    procedure :: get_index, get_indizes
    generic :: index => get_index, get_indizes
  end type dynindex_type

  contains

  subroutine init_dynint(this, unrolled_indices, n_per_col)
    class(dynindex_type), intent(inout) :: this
    integer, intent(in) :: unrolled_indices(:)
    integer, intent(in) :: n_per_col(:)

    call assert(size(unrolled_indices) == sum(n_per_col), 'size(unrolled_indices) /= sum(n_per_col).')

    this%unrolled_indices = unrolled_indices
    this%n_per_col = n_per_col
    this%n_cols = size(this%n_per_col)
    this%n_indices = size(this%unrolled_indices)
  end subroutine init_dynint


  subroutine init_dynint_from_matrix(this, index_matrix, n_per_col)
    class(dynindex_type), intent(inout) :: this
    integer, intent(in) :: index_matrix(:, :)
    integer, intent(in) :: n_per_col(:)

    integer :: col_index, upper_limit, lower_limit

    this%n_cols = size(n_per_col)
    this%n_per_col = n_per_col
    this%n_indices = sum(this%n_per_col)

    call assert(size(index_matrix, 2) == this%n_cols, 'size(index_matrix, 2) /= size(n_per_col).')
    call assert(size(index_matrix, 1) == maxval(n_per_col), 'size(index_matrix, 1) /= maxval(n_per_col).')

    allocate(this%unrolled_indices(this%n_indices))

    do col_index = 1, this%n_cols
      call unrolled_col_limits(this, col_index, lower_limit, upper_limit)
      this%unrolled_indices(lower_limit : upper_limit) = index_matrix(1 : this%n_per_col(col_index), col_index)
    end do
  end subroutine init_dynint_from_matrix

  subroutine finit_dynint(this)
    class(dynindex_type), intent(inout) :: this

    if (allocated(this%unrolled_indices)) deallocate(this%unrolled_indices)
    if (allocated(this%n_per_col)) deallocate(this%n_per_col)
  end subroutine finit_dynint

  function get_col(this, col_index) result(col)
    class(dynindex_type), intent(in) :: this
    integer, intent(in) :: col_index
    integer, allocatable :: col(:)

    integer :: lower_limit, upper_limit

    call assert(col_index >= 1, 'col_index < 1.')
    call assert(col_index <= this%n_cols, 'col_index > this%n_cols.')

    call unrolled_col_limits(this, col_index, lower_limit, upper_limit)

    col = this%unrolled_indices(lower_limit : upper_limit)
  end function get_col

  function get_index(this, row_index, col_index) result(index)
    class(dynindex_type), intent(in) :: this
    integer, intent(in) :: row_index, col_index
    integer :: index

    integer :: unrolled_index

    call assert(col_index >= 1, 'col_index < 1.')
    call assert(col_index <= this%n_cols, 'col_index > this%n_cols.')

    call assert(row_index >= 1, 'row_index < 1.')
    call assert(row_index <= this%n_per_col(col_index), 'row_index > this%n_per_col(col_index).')

    unrolled_index = sum(this%n_per_col(: col_index - 1)) + row_index

    index = this%unrolled_indices(unrolled_index)
  end function    

  function get_indizes(this, row_indizes, col_index) result(indices)
    class(dynindex_type), intent(in) :: this
    integer, intent(in) :: row_indizes(:), col_index
    integer, allocatable :: indices(:)

    integer, allocatable :: unrolled_indices(:)

    call assert(col_index >= 1, 'col_index < 1.')
    call assert(col_index <= this%n_cols, 'col_index > this%n_cols.')

    call assert(all(row_indizes >= 1), 'For some elements: row_index < 1.')
    call assert(all(row_indizes <= this%n_per_col(col_index)), 'For some elements: row_index > this%n_per_col(col_index).')

    unrolled_indices = spread(sum(this%n_per_col(: col_index - 1)), 1, size(row_indizes)) + row_indizes

    indices = this%unrolled_indices(unrolled_indices)
  end function

  subroutine unrolled_col_limits(this, col_index, lower_limit, upper_limit)
    type(dynindex_type), intent(in) :: This
    integer, intent(in) :: col_index 
    integer, intent(out) :: lower_limit, upper_limit

    call assert(col_index >= 1, 'col_index < 1.')
    call assert(col_index <= this%n_cols, 'col_index > this%n_cols.')

    if(col_index == 1) then
      lower_limit = 1
    else   
      lower_limit = sum(this%n_per_col(: col_index - 1)) + 1
    end if   
    upper_limit = sum(this%n_per_col(: col_index))
  end subroutine




end module dynamic_indices