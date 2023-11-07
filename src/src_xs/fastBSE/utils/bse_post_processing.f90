module bse_post_processing
  use precision, only: dp
  use asserts, only: assert
  use modmpi, only: terminate_if_false
  use distributions, only: lorentzian

  private
  public :: absorption_spectrum, symmetrize_quantity
  
  contains

  !> Calculate the absorption spectra from eigen values and oscillator strengths as
  !> \[
  !>     \Im \epsilon(\omega) = \sum_i=1 t_i (f_\broadening(\omega - \varepsilon_i) - f_\broadening(\omega + \varepsilon_i))
  !> \]
  !> where \( \varepsilon_i \) are the eigen energies and \( t_i \) the corresponding oscillator strengths, 'f' is a 
  !> distribution from [[distributions]].
  subroutine absorption_spectrum(energies, weights, omega, broadening, func, spectrum)
    !> Exciton eigen energies
    real(dp), intent(in) :: energies(:)
    !> Exciton oscillator strenghts
    real(dp), intent(in) :: weights(:)
    !> Frequency grid on which the spectrum is calculated
    real(dp), intent(in) ::  omega(:)
    !> Broadening of the peaks 
    real(dp), intent(in) :: broadening
    !> Peak function. Can be a function from [[distributions]] with the same interface.
    interface
      function func(x, gamma, mu) result(lorentzian_func)
        use precision, only: dp
        !> Grid on which the function is calculated
      real(dp), intent(in) :: x(:)
      !> Broadening \(\gamma\)
      real(dp), intent(in) :: gamma
      !> Location of the center of the Lorentzian
      real(dp), intent(in) :: mu

      real(dp), allocatable :: lorentzian_func(:)
      end function
    end interface  
    !> Spectrum.
    real(dp), intent(out) :: spectrum(:)

    integer :: i_energy

    call assert(size(energies) == size(weights), 'size(energies) /= size(weights).')
    call assert(size(omega) == size(spectrum), 'size(omega) /= size(spectrum).')

    spectrum = 0._dp

    do i_energy=1, size(energies)
      spectrum = spectrum + weights(i_energy) &
              * (func(omega, broadening, energies(i_energy)) - func(omega, broadening, -1._dp * energies(i_energy)))
    end do
    
  end subroutine absorption_spectrum

  subroutine symmetrize_quantity(spectrum)
    use modxs, only: symt2
    real(dp), intent(inout) :: spectrum(:, :)

    real(dp) :: symmetric_matrix(3, 3)

    integer :: n_spec, i, j


    call assert(size(spectrum, 2) == 3, 'size(spectrum, 2) /= 3')

    call init0

    do j=1, 3
      do i=1, 3
        symmetric_matrix(i, j) = symt2(i, i, j, j)
      end do 
    end do

    spectrum = transpose(matmul(symmetric_matrix, transpose(spectrum)))

  end subroutine 

end module bse_post_processing 