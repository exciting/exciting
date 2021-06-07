! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! For all new units, please:
! a) Define using CODATA 2018 (see reference below)
! b) Define using double precision (dp): elec_mass = 9.10938370e-31_dp;
! c) Use descriptive naming convention

! See the 20th May 2019 redefinition of the SI units and CODATA 2018.
! Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor (2019)
! "The 2018 CODATA Recommended Values of the Fundamental Physical Constants"
! (Web Version 8.0). Database developed by J. Baker, M. Douma, and S. Kotochigova.
! Available at http://physics.nist.gov/constants,
! National Institute of Standards and Technology, Gaithersburg, MD 20899.

!> Physical constants, defined according to [CODATA 2018])(http://physics.nist.gov/constants) 
module physical_constants
  use precision, only: dp
  use constants, only: pi
  implicit none
  private

  !> Planck constant in J.s (CODATA 2018)   
  real(dp), public, parameter :: hplanck_si = 6.62607015e-34_dp

  !> Reduced Planck constant in J.s (CODATA 2018)
  real(dp), public, parameter :: hbar_si = hplanck_si / (2._dp * pi)

  !> Elementary charge of an electron in C (CODATA 2018) 
  real(dp), public, parameter :: elec_charge = 1.602176634e-19_dp
  
  ! TODO(Alex) Issue #20. Update to CODATA 2018 physical units if required
  !> Boltzmann constant in Hartree/kelvin (CODATA 2006)
  Real(8), Public, Parameter :: kboltz = 3.166815343d-6

  !> Electron mass in kg (CODATA 2018)
  Real(dp), Public, Parameter :: elec_mass = 9.10938370e-31_dp;

  !> Speed of light in m/s (CODATA 2018)
  real(dp), private, parameter:: c0 = 299792458._dp

  !> Bohr radius in m (CODATA 2018)
  real(dp), private, parameter:: a0 = 0.52917721067e-10_dp

  !> 1 a.u. of time in s (CODATA 2018)
  real(dp), private, parameter:: t0 = 2.418884326509e-17_dp

  !> Speed of light (atomic units)
  real(dp), public, parameter :: c = c0*t0/a0

  !> Fine-structure constant (CODATA 2018)
  real(dp), public, parameter :: alpha = 1_dp/137.035999084_dp

  !> Electron g factor (CODATA 2018)
  real(dp), public, parameter :: ge = 2.00231930436256_dp

  !> Bohr radius in metres (CODATA 2018) 
  real(dp), public, parameter :: bohr_radius_si = 5.29177210903e-11_dp

end module physical_constants
