module rttddft_Polarization  
  use constants, only: pi, twopi, zone, real_zero
  use linear_algebra_3d, only: triple_product
  use precision, only: i32, dp
  use rttddft_CurrentDensity, only: Current_Density_Field
  use rttddft_VectorField, only: Uniform_Vector_Field

  implicit none

  private

  !> Polarization vector \(\mathbf{P}\)
  type, public, extends(Uniform_Vector_Field) :: Polarization

  contains 
    procedure :: get_with_mtp, &
      obtain_j
  end type

contains

  !> Evaluate the current as a simple time-derivative \( J(t) = \frac{P(t) - P(t -\Delta t)}{\Delta t} \)
  pure subroutine obtain_j( this, p_vec_prev, dt, current )
    class(Polarization), intent(inout) :: this
    !> Polarization \( P(t -\Delta t) \)
    class(Polarization), intent(in) :: p_vec_prev
    !> Time step \( \Delta t \)
    real(dp), intent(in) :: dt
    !> Resulting current \( J(t) \)
    type(Current_Density_Field), intent(inout) :: current

    current%components = (this%components - p_vec_prev%components) / dt
  end subroutine

  !> Evaluate the infinite system polarization with the "modern theory of polarization"
  !> [Phys. Rev. B 47, 1651(R) (1993), Rev. Mod. Phys. 66, 899 (1994), Phys. Rev. B 69, 085106 (2004)].
  !> Macroscopic polarization vector \( \mathbf{P}_{\alpha} \) along the lattice vector \( \mathbf{a}_{\alpha} \) 
  !> is eveluated as
  !> \[
  !> \mathbf{P}_{\alpha} = \frac{f \mathbf{a}_{\alpha}}{2 \pi V N_{\mathbf{k}_{\alpha}^{\perp}}} 
  !> \sum_{\mathbf{k}_{\alpha}^{\perp}} {\rm Im \, ln} \, \prod_{i}^{N_{\mathbf{k}_{\alpha}} - 1} 
  !> {\rm det} \, S^{\mathbf{k}_i, \mathbf{k}_{i \alpha}^+},
  !> \]
  !> where \( V \) is the unit cell volume, \( N_{\mathbf{k}_{\alpha}^{\perp}} \) is the number of \( \mathbf{k} \) points 
  !> in a plane orthogonal (in the lattice coordinates) to direction \( \alpha \), and \( f \) is the spin degeneracy factor.
  pure subroutine get_with_mtp( this, td_overlap_det, k_grid_dimensions, k_3d_to_1d_map, &
      lattice_vectors, prev_phases, match_phases_with_prev )
    class(Polarization), intent(inout) :: this
    !> Determinants of the overlap matrix \( S^{\mathbf{k}_i, \mathbf{k}_{i \alpha}^+} \) 
    !> between the states corresponding to the neighbouring k points (3, nkpt)
    complex(dp), intent(in) :: td_overlap_det(:, :)
    !> Dimensions of the k points grid (3)
    integer(i32), intent(in) :: k_grid_dimensions(:)
    !> Array which maps 3D integer k point index to the 1D one
    integer(i32), intent(in) :: k_3d_to_1d_map(:, :, :)
    !> Set of the unit cell vectors \( \mathbf{a}_{\alpha} \)
    real(dp), intent(in) :: lattice_vectors(:, :)
    !> String phases evaluated at the previous time step. Used to 
    !> make sure the corect logarithm branch is being used (max_n_ort_plane, 3)
    real(dp), intent(inout) :: prev_phases(:, :)
    !> If .True., prev_phases is used for the phase matching
    logical, intent(in) :: match_phases_with_prev

    integer(i32), parameter :: n_cartesian = 3
    integer(i32) :: direction, n_ort_plane, ort_one, ort_two, i, j, t, cntr, iv(n_cartesian)
    real(dp) :: unit_cell_volume, pvec_lattice(n_cartesian)
    real(dp), allocatable :: phases(:)
    real(dp), parameter :: eps_det = 1.e-14_dp, spin_degeneracy = 2._dp
    complex(dp) :: tmp_c

    this%components = real_zero
    pvec_lattice = real_zero
    unit_cell_volume = triple_product( lattice_vectors )

    allocate( phases(size( prev_phases, 1 )) )
    ! cycle over lattice directions
    do direction = 1, n_cartesian
      if( direction == 1 ) then
        ort_one = 2
        ort_two = 3
      else if( direction == 2 ) then
        ort_one = 1
        ort_two = 3
      else if( direction == 3 ) then
        ort_one = 1
        ort_two = 2
      end if

      n_ort_plane = k_grid_dimensions(ort_one) * k_grid_dimensions(ort_two)
      phases = real_zero
      cntr = 0
      do i = 1, k_grid_dimensions(ort_one)
        iv(ort_one) = i
        do j = 1, k_grid_dimensions(ort_two)
          cntr = cntr + 1
          iv(ort_two) = j

          tmp_c = zone
          do t = 1, k_grid_dimensions(direction)
            iv(direction) = t
            tmp_c = tmp_c * td_overlap_det(direction, k_3d_to_1d_map(iv(1), iv(2), iv(3)))
          end do

          phases(cntr) = aimag( log( tmp_c ) ) ! always lies in (- pi, pi]
          if ( abs( tmp_c ) < eps_det ) phases(cntr) = real_zero

          if ( match_phases_with_prev ) call &
            fixJump( phases(cntr), prev_phases(cntr, direction), pi )

        end do ! orthogonal 1
      end do ! orthogonal 2
      
      ! average over strings
      do cntr = 1, n_ort_plane
        pvec_lattice(direction) = pvec_lattice(direction) + phases(cntr)
      end do
      pvec_lattice(direction) = pvec_lattice(direction) / real( n_ort_plane, dp )

      prev_phases(1 : n_ort_plane, direction) = phases(1 : n_ort_plane)

    end do ! direction

    pvec_lattice = spin_degeneracy * pvec_lattice / ( twopi * unit_cell_volume )

    ! transform to cartesian
    do direction = 1, n_cartesian
      this%components(direction) = this%components(direction) + &
      dot_product( pvec_lattice, lattice_vectors(direction, :) )
    end do

  end subroutine

  !> Match the phase between two values defined up to 2 \pi : 
  !> if they differ by threshold, adds the required number of 2 \pi jumps to the lesser
  pure subroutine fixJump( current, previous, threshold )
    !> Number to match
    real(dp), intent(inout) :: current
    !> Value the number is matched with
    real(dp), intent(in) :: previous
    !> Threshold value of the difference to consider the values 'not matched'
    real(dp), intent(in) :: threshold

    real(dp) :: delta

    delta = abs( current - previous )
    if ( delta > threshold ) then
      if ( current > previous ) then
        current = current - twopi * nint( delta / twopi )
      else
        current = current + twopi * nint( delta / twopi )
      end if
    end if

  end subroutine fixJump

end module