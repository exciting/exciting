!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!

!> Converts a real function \(r_{lm}\), expanded in terms of real spherical harmonics
!> into a complex spherical harmonics expansion \(z_{lm}\) with
!> \[ z_{lm} = \begin{cases}
!>    \frac{1}{\sqrt{2}} \left( (-1)^m\, r_{l-m} - {\rm i}\, r_{lm} \right) & m < 0 \\
!>    r_{lm} & m = 0 \\
!>    \frac{1}{\sqrt{2}} \left( r_{lm} + {\rm i}\, (-1)^m\, r_{l-m} \right) & m > 0
!>    \end{cases} \;. \]
subroutine rtozflm(lmax, rflm, zflm)
  use precision, only: dp
  use modmpi, only: terminate_if_false
  implicit none
  !> maximum angular momentum
  integer, intent(in) :: lmax
  !> coefficients of real spherical harmonic expansion
  real(dp), intent(in) :: rflm(*)
  !> coefficients of complex spherical harmonic expansion
  complex(dp), intent(out) :: zflm(*)

  real(dp), parameter :: one_over_root2 = 0.7071067811865475244_dp

  integer :: l, m, lm, lm_pos, lm_neg

  call terminate_if_false( lmax >= 0, &
    '(rtozflm) `lmax` must not be negative.')

  do l = 0, lmax
    lm = l**2 + l + 1
    ! m = 0
    zflm(lm) = cmplx( rflm(lm), 0._dp, dp)
    ! m /= 0
    do m = 1, l
      lm_pos = lm + m
      lm_neg = lm - m
      zflm(lm_pos) = one_over_root2 * cmplx(  rflm(lm_pos), (-1)**m * rflm(lm_neg), dp)
      zflm(lm_neg) = one_over_root2 * cmplx( (-1)**m * rflm(lm_pos), -rflm(lm_neg), dp)
    end do
  end do
end subroutine rtozflm

!> Compose two real functions \(r^R_{lm}\) and \(r^I_{lm}\) expanded in terms of real
!> spherical harmonics and representing real and imaginary part of a complex function
!> into a complex spherical harmonics expansion \(z_{lm}\) with
!> \[ z_{lm} = \begin{cases}
!>    \frac{1}{\sqrt{2}} \left( (-1)^m\, r^R_{l-m} + r^I_{lm} + 
!>      {\rm i} \left( -r^R_{lm} + (-1)**m r^I_{l-m} \right) \right) & m < 0 \\
!>    r^R_{lm} + {\rm i}\, r^I_{lm} & m = 0 \\
!>    \frac{1}{\sqrt{2}} \left( r^R_{lm} - (-1)**m r^I_{l-m} + 
!>       {\rm i} \left( (-1)**m r^R_{l-m} + r^I_{lm} \right) \right) & m > 0
!>    \end{cases} \;. \]
subroutine compose_zflm(lmax, rflm_real, rflm_imag, zflm)
  use precision, only: dp
  use modmpi, only: terminate_if_false
  implicit none
  !> maximum angular momentum
  integer, intent(in) :: lmax
  !> coefficients of real part real spherical harmonic expansion
  real(dp), intent(in) :: rflm_real(*)
  !> coefficients of imaginary part real spherical harmonic expansion
  real(dp), intent(in) :: rflm_imag(*)
  !> coefficients of complex spherical harmonic expansion
  complex(dp), intent(out) :: zflm(*)

  real(dp), parameter :: one_over_root2 = 0.7071067811865475244_dp

  integer :: l, m, lm, lm_pos, lm_neg

  call terminate_if_false( lmax >= 0, &
    '(compose_zflm) `lmax` must not be negative.')

  do l = 0, lmax
    lm = l**2 + l + 1
    ! m = 0
    zflm(lm) = cmplx( rflm_real(lm), rflm_imag(lm), dp)
    ! m /= 0
    do m = 1, l
      lm_pos = lm + m
      lm_neg = lm - m
      zflm(lm_pos) = one_over_root2 * cmplx(  rflm_real(lm_pos) - (-1)**m * rflm_imag(lm_neg), &
                                              (-1)**m * rflm_real(lm_neg) + rflm_imag(lm_pos), dp)
      zflm(lm_neg) = one_over_root2 * cmplx(  (-1)**m * rflm_real(lm_pos) + rflm_imag(lm_neg), &
                                             -rflm_real(lm_neg) + (-1)**m * rflm_imag(lm_pos), dp)
    end do
  end do
end subroutine compose_zflm
