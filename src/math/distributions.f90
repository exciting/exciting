!> Distributions.
module distributions
    use precision, only: dp
    use constants, only: pi

  
    implicit none
    
    private
    public :: lorentzian
contains

  !> Constructs the Lorentzian distribution, defined as 
  !>
  !> \[
  !>    f(x; \mu,\gamma) = 
  !>        \frac{1}{\pi \gamma}\frac{\gamma^2}{(x - \mu)^2 + \gamma^2}
  !>                                                                          \]
  !>
  !> for an array of values \( x \), the location parameter \( \mu  \)
  !> and the broadening \(  \gamma \).
  !> See e.g. https://en.wikipedia.org/wiki/Cauchy_distribution.
  pure function lorentzian(x, gamma, mu) result(lorentzian_func)

    !> Grid on which the function is calculated
    real(dp), intent(in) :: x(:)
    !> Broadening \(\gamma\)
    real(dp), intent(in) :: gamma
    !> Location of the center of the Lorentzian
    real(dp), intent(in) :: mu

    real(dp), allocatable :: lorentzian_func(:)

    lorentzian_func = gamma/((x - mu)**2 + gamma**2) /pi
  end function

  !> Constructs the Gauss distribution, defined as 
  !>
  !> \[
  !>    f(x; \mu,\sigma) = 
  !>        \frac {1} {\sqrt{2 \pi \sigma^2}} \exp{ \frac {(x - \mu)^2} {\sqrt{2 \sigma^2}} }
  !>                                                                          \]
  !>
  !> for an array of values \( x \), the location parameter \( \mu  \)
  !> and the broadening \(  \sigma \).
  !> See e.g. https://en.wikipedia.org/wiki/Normal_distribution.
  pure function gaussian(x, sigma, mu) result(gaussian_func)

    !> Grid on which the function is calculated
    real(dp), intent(in) :: x(:)
    !> Broadening \(\sigma\)
    real(dp), intent(in) :: sigma
    !> Location of the center of the Gaussian
    real(dp), intent(in) :: mu

    real(dp), allocatable :: gaussian_func(:)
    
    gaussian_func = 1.0_dp / sqrt(2 * pi * sigma**2) * exp((x - mu)**2 / sqrt( 2* sigma**2))
  end function gaussian

end module distributions