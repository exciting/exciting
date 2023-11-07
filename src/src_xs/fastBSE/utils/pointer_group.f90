module pointer_group
  use precision, only: dp

  type matrix_pointer_group_cmplx_dp
    complex(dp), pointer :: p(:, :)

    contains 

    procedure :: finalize
  end type matrix_pointer_group_cmplx_dp

  contains 

  subroutine finalize(this)
    class(matrix_pointer_group_cmplx_dp) :: this

    this%p => null()
  end subroutine finalize


end module pointer_group