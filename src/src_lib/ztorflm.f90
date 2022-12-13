!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!

!> Converts a real function \(z_{lm}\) expanded in terms of complex spherical harmonics
!> into a real spherical harmonics expansion \(r_{lm}\) with
!> \[ r_{lm} = \begin{cases}
!>    \frac{1}{\sqrt{2}} \Im \left\lbrace -z_{lm} + (-1)^m\, z_{l-m} \right\rbrace & m < 0 \\
!>    \Re \lbrace z_{lm \rbrace & m = 0 \\
!>    \frac{1}{\sqrt{2}} \Re \left\lbrace  z_{lm} + (-1)^m\, z_{l-m} \right\rbrace & m > 0
!>    \end{cases} \;. \]
subroutine ztorflm(lmax, zflm, rflm)
  use precision, only: dp
  use modmpi, only: terminate_if_false
  implicit none
  !> maximum angular momentum
  integer, intent(in) :: lmax
  !> coefficients of complex spherical harmonic expansion
  complex(dp), intent(in) :: zflm(*)
  !> coefficients of real spherical harmonic expansion
  real(dp), intent(out) :: rflm(*)

  real(dp), parameter :: one_over_root2 = 0.7071067811865475244_dp

  integer :: l, m, lm, lm_pos, lm_neg

  call terminate_if_false( lmax >= 0, &
    '(rtozflm) `lmax` must not be negative.')

  do l = 0, lmax
    lm = l**2 + l + 1
    ! m = 0
    rflm(lm) = dble(zflm(lm))
    ! m /= 0
    do m = 1, l
      lm_pos = lm + m
      lm_neg = lm - m
      rflm(lm_pos) = one_over_root2 * ( dble( zflm(lm_pos)) + (-1)**m * dble( zflm(lm_neg)))
      rflm(lm_neg) = one_over_root2 * (-aimag(zflm(lm_neg)) + (-1)**m * aimag(zflm(lm_pos)))
    end do
  end do
end subroutine ztorflm

!> Decompose a complex function \(z_{lm}\) expanded in terms of complex spherical harmonics
!> into two real spherical harmonic expansions \(r^R_{lm}\) and \(r^I_{lm}\) representing
!> the real and imaginary part of the complex function with
!> \[ r^R_{lm} = \begin{cases}
!>    \frac{1}{\sqrt{2}} \Im \left\lbrace -z_{lm} + (-1)^m\, z_{l-m} \right\rbrace & m < 0 \\
!>    \Re \lbrace z_{lm} \rbrace & m = 0 \\
!>    \frac{1}{\sqrt{2}} \Re \left\lbrace z_{lm} + (-1)^m\, z_{l-m} \right\rbrace & m > 0
!>    \end{cases} \, \]
!> and
!> \[ r^I_{lm} = \begin{cases}
!>    \frac{1}{\sqrt{2}} \Re \left\lbrace z_{lm} - (-1)^m\, z_{l-m} \right\rbrace & m < 0 \\
!>    \Im \lbrace z_{lm} \rbrace & m = 0 \\
!>    \frac{1}{\sqrt{2}} \Im \left\lbrace z_{lm} + (-1)^m\, z_{l-m} \right\rbrace & m > 0
!>    \end{cases} \;. \]
subroutine decompose_zflm(lmax, zflm, rflm_real, rflm_imag)
  use precision, only: dp
  use modmpi, only: terminate_if_false
  implicit none
  !> maximum angular momentum
  integer, intent(in) :: lmax
  !> coefficients of complex spherical harmonic expansion
  complex(dp), intent(in) :: zflm(*)
  !> coefficients of real part real spherical harmonic expansion
  real(dp), intent(out) :: rflm_real(*)
  !> coefficients of imaginary part real spherical harmonic expansion
  real(dp), intent(out) :: rflm_imag(*)

  real(dp), parameter :: one_over_root2 = 0.7071067811865475244_dp

  integer :: l, m, lm, lm_pos, lm_neg

  call terminate_if_false( lmax >= 0, &
    '(compose_zflm) `lmax` must not be negative.')

  do l = 0, lmax
    lm = l**2 + l + 1
    ! m = 0
    rflm_real(lm) = dble( zflm(lm))
    rflm_imag(lm) = aimag(zflm(lm))
    ! m /= 0
    do m = 1, l
      lm_pos = lm + m
      lm_neg = lm - m
      rflm_real(lm_pos) = one_over_root2 * ( dble( zflm(lm_pos)) + (-1)**m * dble( zflm(lm_neg)))
      rflm_real(lm_neg) = one_over_root2 * (-aimag(zflm(lm_neg)) + (-1)**m * aimag(zflm(lm_pos)))
      rflm_imag(lm_pos) = one_over_root2 * ( aimag(zflm(lm_pos)) + (-1)**m * aimag(zflm(lm_neg)))
      rflm_imag(lm_neg) = one_over_root2 * ( dble( zflm(lm_neg)) - (-1)**m * dble( zflm(lm_pos)))
    end do
  end do
end subroutine decompose_zflm
