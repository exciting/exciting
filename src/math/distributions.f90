!> Tools that are useful for defining grids or working with grids.
module distributions
    use precision, only: sp, dp

  
    implicit none
    
    private
    public :: lorentzian
contains
        !> Constructs the Lorentzian distribution, defined as 
    !>
    !> \[
    !>    f(\omega;\omega_0,\gamma) = 
    !>        \frac{1}{\pi\gamma}\frac{\gamma^2}{(\omega - \omega_0)^2 + \gamma^2}
    !>                                                                          \]
    !>
    !> for given frequencies \( \omega \), location parameter \( \omega_0  \)
    !> and broadening \(  \gamma \).
    !> See e.g. https://en.wikipedia.org/wiki/Cauchy_distribution.
    function lorentzian(broad, frequencies, freq_0) result(lorentzian_func)

        use constants, only: pi

        !> Broadening in eV
        real(dp) :: broad
        !> Location parameter, i.e. center of the Lorentzian
        real(dp) :: freq_0
        !> Frequencies of the frequency grid
        real(dp) :: frequencies(:)
        real(dp), allocatable :: lorentzian_func(:)
        !> Number of frequencies
        integer :: n_freqs
        !> Runnind index frequencies
        integer :: ifreq

        n_freqs = size(frequencies)

        allocate (lorentzian_func(size(frequencies)))

        lorentzian_func = broad/((frequencies - freq_0)**2 + broad**2) /pi



    end function
end 