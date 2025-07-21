!> Collection of different occupation functions.
module occupation_functions
  use precision, only: dp
  use physical_constants, only: kboltz
  
  implicit none
  private

  public :: fermi_dirac, bose_einstein

contains

  !> Fermi-Dirac occupation
  !> \[ f(e,T) = \frac{1}{\exp\left(\frac{e}{k_{\rm B}\, T}\right) + 1} \]
  pure function fermi_dirac( e, T ) result( occ )
    !> energy in Hartree (relative to Fermi energy / chemical potential)
    real(dp), intent(in) :: e
    !> temperature in Kelvin (negative values will be treated as zero)
    real(dp), intent(in) :: T
    !> occupation
    real(dp) :: occ
  
    real(dp), parameter :: max_exp = log( huge( 1.0_dp ) )

    real(dp) :: et, x

    et = T * kboltz

    ! zero or negative temperature
    if (et < tiny(1.0_dp)) then
      if (e == 0.0_dp) then
        occ = 0.5_dp
      else if (e > 0.0_dp) then
        occ = 0.0_dp
      else
        occ = 1.0_dp
      end if
    ! finite temperature
    else
      x = e / et
      if (x > max_exp) then
        occ = 0.0_dp
      else if (x < -max_exp) then
        occ = 1.0_dp
      else
        occ = 1.0_dp / (exp( x ) + 1.0_dp)
      end if
    end if
  end function fermi_dirac

  !> Bose-Einstein occupation
  !> \[ n(e,T) = \frac{1}{\exp\left(\frac{e}{k_{\rm B}\, T}\right) - 1} \]
  pure function bose_einstein( e, T ) result( occ )
    !> energy in Hartree (relative to Fermi energy / chemical potential)
    real(dp), intent(in) :: e
    !> temperature in Kelvin (negative values will be treated as zero)
    real(dp), intent(in) :: T
    !> occupation
    real(dp) :: occ
  
    real(dp), parameter :: max_exp = log( huge( 1.0_dp ) )

    real(dp) :: et, x

    et = T * kboltz

    ! zero or negative temperature
    if (et < tiny(1.0_dp)) then
      occ = 0.0_dp
    ! finite temperature
    else
      x = e / et
      if (x > max_exp .or. x < 0.0_dp) then
        occ = 0.0_dp
      else if (x < tiny(1.0_dp)) then
        occ = huge(1.0_dp)
      else
        occ = 1.0_dp / (exp( x ) - 1.0_dp)
      end if
    end if
  end function bose_einstein

end module occupation_functions
