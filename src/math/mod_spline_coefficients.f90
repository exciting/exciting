
module mod_spline_coefficients

    use precision, only: dp

    implicit none

    private
    public :: get_property_from_spline_coefficients

contains

!> Calculates a property at point $r_sub$ from spline coefficients at $r /leq r_sub$. 
pure function get_property_from_spline_coefficients(A_r, cf, r, r_sub) result(A_rsub)
    !> property at grid point r
    real(dp), intent(in) :: A_r
    !> cubic spline coefficients (1,2,3) and work space (4) at grid point r
    real (dp), intent(in)  :: cf(4)
    !> grid point 
    real(dp), intent (in) :: r
    !> subgrid point
    real(dp), intent(in) :: r_sub
    !> property at subgridpoint
    real(dp) :: A_rsub 

    ! local variables
    real(dp) :: delta_r

    delta_r = r_sub - r
    A_rsub = A_r + delta_r * (cf(1) + delta_r * (cf(2) + delta_r * cf(3)))

end function get_property_from_spline_coefficients

end module mod_spline_coefficients
