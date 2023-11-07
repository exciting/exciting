module isdf_test_utils
  use precision, only: dp 
  use asserts, only: assert
  use multi_index_conversion, only: composite_index_to_indices
  use xlapack, only: matrix_multiply

  private
  public :: calculate_wavefunction_product_samek, calculate_isdf_wavefunction_product_samek

  contains 


  subroutine calculate_wavefunction_product_samek(wf_o, wf_u, wf_product)
    complex(dp), intent(in), contiguous :: wf_o(:, :, :), wf_u(:, :, :)
    complex(dp), intent(out), contiguous :: wf_product(:, :, :, :)
    
    integer :: n_k, n_o, n_u, n_r
    integer :: ik, io, iu, i_composite

    n_r = size(wf_o, 1)
    n_o = size(wf_o, 2)
    n_u = size(wf_u, 2)
    n_k = size(wf_o, 3)

    call assert(size(wf_u, 1) == n_r, 'size(u_u, 1) /= n_r.')
    call assert(size(wf_u, 3) == n_k, 'size(u_u, 3) /= n_k.')
    call assert(size(wf_product, 1) == n_r, 'size(wf_product, 1) /= n_r.')
    call assert(size(wf_product, 2) == n_o, 'size(wf_product, 2) /= n_o.')
    call assert(size(wf_product, 3) == n_u, 'size(wf_product, 3) /= n_u.')
    call assert(size(wf_product, 4) == n_k, 'size(wf_product, 4) /= n_k.')

    do i_composite=1, n_k * n_u * n_o 
      call composite_index_to_indices(i_composite, [n_u, n_o, n_k], iu, io, ik)

      wf_product(:, io, iu, ik) = conjg(wf_o(:, io, ik)) * wf_u(:, iu, ik)
    end do 

  end subroutine calculate_wavefunction_product_samek


  subroutine calculate_isdf_wavefunction_product_samek(zeta, wf_o_isdf, wf_u_isdf, wf_product)
    complex(dp), intent(in), contiguous :: zeta(:, :), wf_u_isdf(:, :, :), wf_o_isdf(:, :, :)
    complex(dp), intent(out), contiguous :: wf_product(:, :, :, :)

    integer :: n_k, n_o, n_u, n_r, n_isdf 
    integer :: ik, io, iu, i_composite

    n_r = size(zeta, 1)
    n_isdf = size(zeta, 2)
    n_o = size(wf_o_isdf, 2)
    n_u = size(wf_u_isdf, 2)
    n_k = size(wf_o_isdf, 3)

    call assert(size(wf_o_isdf, 1) == n_isdf, 'size(u_u, 1) /= n_isdf.')
    call assert(size(wf_u_isdf, 1) == n_isdf, 'size(u_u, 1) /= n_isdf.')
    call assert(size(wf_u_isdf, 3) == n_k, 'size(u_u, 3) /= n_k.')
    call assert(size(wf_product, 1) == n_r, 'size(wf_product, 1) /= n_r.')
    call assert(size(wf_product, 2) == n_o, 'size(wf_product, 2) /= n_o.')
    call assert(size(wf_product, 3) == n_u, 'size(wf_product, 3) /= n_u.')
    call assert(size(wf_product, 4) == n_k, 'size(wf_product, 4) /= n_k.')

    do i_composite=1, n_k * n_u * n_o 
      call composite_index_to_indices(i_composite, [n_u, n_o, n_k], iu, io, ik)

      call matrix_multiply(zeta, conjg(wf_o_isdf(:, io, ik)) * wf_u_isdf(:, iu, ik), wf_product(:, io, iu, ik))
    end do 
  end subroutine 

end module isdf_test_utils