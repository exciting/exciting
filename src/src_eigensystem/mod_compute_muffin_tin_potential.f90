! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) exciting Code, SOL group. 2025

! Created May 2025 (Mara Voiculescu)

!> The matrix elements of the muffin-tin (MT) potential contribution to the Hamiltonian matrix depend on terms of the
!> following form:
!> \[
!>      \langle \varphi^{\alpha}_{\lambda} | v^{\phantom{I}}_{\texttt{KS}} | \varphi^{\alpha}_{\lambda'} \rangle =
!>      \sum_{l m} G^{l m}_{l_{\lambda}, m_{\lambda}, l_{\lambda'}, m_{\lambda'}} \int_{0}^{r^{\texttt{MT}}_{\alpha}}
!>      \! f^{\alpha *}_{\lambda}(r) \, v_{l m}(r) \, f^{\alpha}_{\lambda'}(r) r^2 \, dr \; .
!> \]
!> where the functions \( f^{\alpha}_{\lambda}(r) \) refer to the radial part of the MT basis function and the
!> Gaunt coefficients  \( G^{l m}_{l_{\lambda}, m_{\lambda}, l_{\lambda'}, m_{\lambda'}} \) are defined as
!> \[
!>      G^{l m}_{l_{\lambda}, m_{\lambda}, l_{\lambda'}, m_{\lambda'}} =
!>      \int Y_{l_{\lambda} m_{\lambda}\!}(\hat{\mathbf{r}}) \, \mathcal{S}_{lm}(\hat{\mathbf{r}})
!>      \,Y_{l_{\lambda'} m_{\lambda'}\!}(\hat{\mathbf{r}}) d\mathbf{r} \; .
!> \]
!> This module gathers the necessary code for computing the radial integrals and the matrix elements in the
!> subroutine [[mt_pot]].
module mod_compute_muffin_tin_potential
  use precision, only: dp, i32

  implicit none
  private
  public :: integrate_mt_potential, build_mt_potential_matrix_element

contains

  !> This routine computes the radial integrals needed for the potential matrix elements.
  subroutine integrate_mt_potential(num_sph_harm, radial_mt_basis_function_bra, &
       radial_mt_basis_function_ket, num_rad_grid_points, mt_potential, rad_mesh, radial_integral)
    !> Number of spherical harmonics in potential expansion
    integer(i32), intent(in) :: num_sph_harm
    !> Radial part of MT basis function for a given l-value and index of atoms and species, corresponding to the
    !> "bra" state of the potential matrix element
    real(dp), intent(in), contiguous :: radial_mt_basis_function_bra(:)
    !> Radial part of MT basis function for a given l-value and index of atoms and species, corresponding to the
    !> "ket" state of the potential matrix element
    real(dp), intent(in), contiguous :: radial_mt_basis_function_ket(:)
    !> Number of radial grid points
    integer, intent(in) :: num_rad_grid_points
    !> Muffin-tin potential for a given index of atoms and species
    real(dp), intent(in), contiguous :: mt_potential(:, :)
    !> Radial-mesh array for a given species index
    real(dp), intent(in), contiguous :: rad_mesh(:)
    !> Output integral
    real(dp), intent(out) :: radial_integral(num_sph_harm)

    ! Local variables
    integer(i32) :: i_sh
    real(dp) :: integrand(num_rad_grid_points), integral(num_rad_grid_points), spline_coeffs(3, num_rad_grid_points)

    do i_sh = 1, num_sph_harm
       integrand(:) = radial_mt_basis_function_bra(1:num_rad_grid_points) * &
            radial_mt_basis_function_ket(1:num_rad_grid_points) * (rad_mesh(1:num_rad_grid_points)) ** 2 * &
            mt_potential(i_sh,1:num_rad_grid_points)
       call fderiv(-1, num_rad_grid_points, rad_mesh, integrand, integral, spline_coeffs)
       radial_integral(i_sh) = integral(num_rad_grid_points)
    end do

  end subroutine integrate_mt_potential

  !> This routine computes the values of the MT potential matrix elements.
  pure subroutine build_mt_potential_matrix_element(num_sph_harm, gaunt_coefficient, radial_integral, matrix_element)
    !> Number of spherical harmonics in potential expansion
    integer(i32), intent(in) :: num_sph_harm
    !> Gaunt coefficient array for given lm-values of the complex spherical harmonics
    complex(dp), intent(in), contiguous :: gaunt_coefficient(:)
    !> Radial integral array
    real(dp), intent(in), contiguous :: radial_integral(:)
    !> Output potential matrix element contribution
    complex(dp), intent(out) :: matrix_element

    matrix_element = sum(gaunt_coefficient(1:num_sph_harm) * radial_integral(1:num_sph_harm))

  end subroutine build_mt_potential_matrix_element

end module mod_compute_muffin_tin_potential
