!> Compute the one-particle Greens function in the band representation. 
module greens_bands 
  use precision, only: wp 
  !use constants, only: i_cmplx => zi 
  use time_freq_grid, only: grid_type 
  implicit none 
  private 
  public :: g0_bands_img_time 
 
  contains 
 
  !
  !> G0 in bands representation with finite temperature (i.e., with beta)
  !> \[
  !>    G_{n n^{\prime}}^{0\bf k}(\tau) = - \delta_{n n^{\prime}} 
  !>      \frac{e^{\left(\epsilon_n^{\mathbf{k}}-\mu\right)(\beta-\tau)}}
  !>      {1+e^{\left(\epsilon_n^{\mathbf{k}}-\mu\right) \beta}}
  !> \]
  !> 
  !> Construct G0 in reciprocal space, in the bands representation. 
  !> Zero-temperature formalism. 
  !> \[
  !>    G_{n n^{\prime}}^{0\bf k}(\tau) =
  !>     - \delta_{n n^{\prime}} e^{- \left(\epsilon_n^{\mathbf{k}}-\mu\right)\tau}
  !>       \Theta(-\tau)\Theta(\mu-\epsilon_n^{\mathbf{k}})
  !>     - \delta_{n n^{\prime}} e^{- \left(\epsilon_n^{\mathbf{k}}-\mu\right)\tau}
  !>       \Theta(\tau)\Theta(\epsilon_n^{\mathbf{k}}-\mu)
  !> \]
  ! 
  subroutine g0_bands_img_time(eigenvalues, mu, img_time_grid, n_occupied, G0_bands) 
    !> KS eigenvalues 
    real(wp), intent(in) :: eigenvalues(:, :) 
    !> Fermi level/chemical potential 
    real(wp), intent(in) :: mu 
    !> Imagainary time grid 
    type(grid_type) :: img_time_grid 
    !> Number of occupied states (including core states) 
    integer, intent(in) :: n_occupied 
    !> G0 in bands representation (diagonal in KS basis) 
    complex(wp), intent(out) :: G0_bands(:, :, :) 
 
    !> Imaginary grid value 
    real(wp) :: tau 
    !> Number of k-points 
    integer :: n_k 
    !> Indices 
    integer :: i_tau, ik, in 
    !> Time grid limits 
    integer :: start, end 
 
    n_k = size(eigenvalues, 2) 
 
    ! Negative times and occupied states 
    call img_time_grid%negative_point_limits(start, end) 
 
    do i_tau = start, end 
      tau = img_time_grid%points(i_tau) 
      do ik = 1, n_k 
        do in = 1, n_occupied 
          G0_bands(in, ik, i_tau) = - exp(-(eigenvalues(in, ik) - mu) * tau) 
        enddo
      enddo 
    enddo 
  
    ! Zero-time for all states only one of terms (1st term on the RHS) contribute to G0 
    !> \[
    !>    G_{n n^{\prime}}^{0\bf k}(\tau = 0) = - \delta_{n n^{\prime}} 
    !> \]
    if (img_time_grid%has_zero_point()) then 
      i_tau = img_time_grid%zero_point_index 
      G0_bands(in, :, i_tau) = - 1._wp 
    endif 
     
    ! Positive times and unoccupied states 
    call img_time_grid%positive_point_limits(start, end) 
 
    do i_tau = start, end 
      tau = img_time_grid%points(i_tau) 
      do ik = 1, n_k 
        do in = n_occupied + 1, size(G0_bands, 1) 
          G0_bands(in, ik, i_tau) = - exp(-(eigenvalues(in, ik) - mu) * tau) 
        enddo
      enddo 
    enddo 
 
  end subroutine g0_bands_img_time 
 
end module greens_bands 
