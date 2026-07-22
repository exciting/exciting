!> Module that deals with the Current Density in RT-TDDFT calculations
module rttddft_CurrentDensity
  use constants, only: real_zero, zzero
  use exciting_mpi, only: mpiinfo, xmpi_allreduce
  use mod_lattice, only: Omega
  use physical_constants, only: c
  use precision, only: i32, dp
  use rttddft_pmat, only: pmat_set
  use rttddft_VectorField, only: Uniform_Vector_Field
  use rttddft_VectorPotential, only: Vector_Potential_Field
  use xlapack, only: dot_multiply, hermitian_matrix_multiply, matrix_multiply
  use rttddft_Wavefunction, only: wavefunction_set

  implicit none

  private

  !> Number of cartesian directions
  integer(i32), parameter :: n_cartesian = 3

  !> Type to store any component of current density as a vector with `x`, `y`, and `z` components
  type, public, extends(Uniform_Vector_Field) :: Current_Density_Field
  end type

  !> Current density stored as a vector with the components along `x`, `y`, and `z`
  type, public :: Current_Density
    !> Diagmagnetic part of the current density
    !> \(\mathbf{J}_{\rm dia} = -\frac{N}{c\Omega}\mathbf{A}\)
    type(Current_Density_Field) :: diamagnetic
    !> Paramagnetic part of the current density
    type(Current_Density_Field) :: paramagnetic
    !> Spurious paramagnetic current density obtained at \( t = 0 \).
    type(Current_Density_Field), private :: spurious
    !> Paramagnetic current density contribution from the frozen states, to be evaluated only once.
    type(Current_Density_Field), private :: frozen
    !> If .true., frozen contribution was already calculated.
    logical, private :: frozen_calculated = .false.
  contains
    procedure, public :: evaluate_paramagnetic => Current_Density_evaluate_paramagnetic
    procedure, public :: evaluate_diamagnetic => Current_Density_evaluate_diamagnetic
    procedure, public :: set_spurious => Current_Density_set_spurious
    procedure, public :: total => Current_Density_total
    procedure, public :: total_components => Current_Density_total_components
  end type

contains
  !> Total current density = paramagnetic + diamagnetic
  !> Result is an object of type `Current_Density_Field`
  pure function Current_Density_total( this ) result( r )
    class(Current_Density), intent(in) :: this
    type(Current_Density_Field) :: r

    r%components = this%paramagnetic%components + this%diamagnetic%components
  end function

  !> Total current density = paramagnetic + diamagnetic
  !> Result is an array with the components
  pure function Current_Density_total_components( this ) result( r )
    class(Current_Density), intent(in) :: this
    real(dp) :: r(n_cartesian)

    r = this%paramagnetic%components + this%diamagnetic%components
  end function

  !> Evaluate the diamagnetic current density as
  !> \[ \mathbf{J}_{ind}(t) = - \frac{N_{\rm active} \mathbf{A}_{tot}(t)}{\Omega c} \]
  !> \(N_{\rm active}\) is the number of active valence electrons, \(c\) is the light speed, and
  !> \(\Omega\) is the unit cell volume
  pure subroutine Current_Density_evaluate_diamagnetic( this, N_active_el_per_volume, a_tot )
    class(Current_Density), intent(inout) :: this
    !> Number of active valence electrons per volume
    real(dp), intent(in) :: N_active_el_per_volume
    !> Vector potential
    class(Vector_Potential_Field), intent(in) :: a_tot

    this%diamagnetic%components = ( -N_active_el_per_volume / c ) * a_tot%components 
  end subroutine

  !> Calculate the paramagnetic part of the current density at time \( t \) as:
  !> \[
  !>    \mathbf{J}(t) = \frac{\mathrm{i}}{\Omega} \sum_{j\mathbf{k}}
  !>      w_{\mathbf{k}}f_{j\mathbf{k}} \left\langle \psi_{j\mathbf{k}}(t) \big|
  !>      \nabla \big|\psi_{j\mathbf{k}}(t)\right\rangle - \mathbf{J}_{\rm spurious}
  !>  \]
  !> where index \( j \) runs over all (frozen and active) states, 
  !> \( \Omega \) is the unit cell volume , \( w_{\mathbf{k}} \) is the 
  !> \( \mathbf{k} \)-point weight, \( f_{j\mathbf{k}} \) is the occupation number of the
  !> corresponding KS state, and \mathbf{J}_{\rm spurious} is the spurious current 
  !> attributed to numerical inaccuracy.
  subroutine Current_Density_evaluate_paramagnetic( this, psi, p_mat, mpi_env )
    class(Current_Density), intent(inout) :: this
    !> Basis-expansion coefficients of the KS-wavefunctions at time \( t \)
    class(wavefunction_set), intent(in) :: psi
    !> Momentum matrix elements
    class(pmat_set), intent(in) :: p_mat
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    integer(i32) :: ik, ist, j
    real(dp) :: aux(n_cartesian)
    real(dp), allocatable :: acc(:)
    complex(dp), allocatable :: draft(:, :)
    real(dp), parameter :: tol_default = 1e-6_dp
    logical :: assume_hermitian

    assume_hermitian = p_mat%is_hermitian()
    call wrapper_over_cartesian_loops( psi%active, psi%first_active(), aux )
    this%paramagnetic%components = aux / Omega
    call xmpi_allreduce( this%paramagnetic%components, mpi_env )

    if ( .not. this%frozen_calculated ) then
      if ( psi%n_frozen() > 0 ) then
        call wrapper_over_cartesian_loops( psi%frozen, 1, aux )
        this%frozen%components = aux / Omega
        call xmpi_allreduce( this%frozen%components, mpi_env )
      else
        this%frozen%components = real_zero
      end if
      this%frozen_calculated = .true.
    end if
    call this%paramagnetic%add_vector( this%frozen%components )
    call this%paramagnetic%add_vector( -this%spurious%components )
    
    contains
      subroutine wrapper_over_cartesian_loops( wavefunction, occs_offset, aux_out )
        complex(dp), contiguous, intent(in) :: wavefunction(:, :, :)
        integer(i32), intent(in) :: occs_offset
        real(dp), intent(out) :: aux_out(:)
        
        integer(i32) :: j, n_states
        
        n_states = size( wavefunction, 2 )
        allocate( draft(psi%n_basis(), n_states), acc(n_states) )

        do j = 1, n_cartesian
          call loop_for_selected_cartesian_component( &
              p_mat%components(j)%array, wavefunction, occs_offset, aux_out(j) )
        end do

        deallocate( draft, acc )
      end subroutine wrapper_over_cartesian_loops

      subroutine loop_for_selected_cartesian_component( p_mat_component, &
        wavefunction, occs_offset, aux_j )
        complex(dp), contiguous, intent(in) :: p_mat_component(:, :, :)
        complex(dp), contiguous, intent(in) :: wavefunction(:, :, :)
        integer(i32), intent(in) :: occs_offset
        real(dp), intent(out) :: aux_j
        integer(i32) :: shift, n_kpts, n_states
        
        shift = psi%first_kpt() - 1
        n_kpts = psi%last_kpt() - psi%first_kpt() + 1
        aux_j = real_zero
        n_states = size( wavefunction, 2 )
        !$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(ik, ist, draft, acc), REDUCTION(+:aux_j), &
        !$OMP& SHARED(assume_hermitian, n_kpts, n_states, occs_offset, wavefunction, p_mat_component, psi, shift)
        do ik = 1, n_kpts
          draft = zzero; acc = real_zero
          if ( assume_hermitian ) then
            call hermitian_matrix_multiply( p_mat_component(:, :, ik), wavefunction(:, :, ik), draft, 'U', 'L', tol_default )
          else
            call matrix_multiply( p_mat_component(:, :, ik), wavefunction(:, :, ik), draft )
          end if
          do ist = 1, n_states
            acc(ist) = real( dot_multiply( wavefunction(:, ist, ik), draft(:, ist), conjg_a=.true. ), dp )
          end do
          aux_j = aux_j - dot_multiply( psi%occupations(occs_offset : occs_offset + n_states - 1, &
            ik + shift), acc ) * psi%kset%wkpt(ik + shift)
        end do
        !$OMP END PARALLEL DO
      end subroutine
  end subroutine

  !> Set the spurious component to a user-provided value.
  pure subroutine Current_Density_set_spurious( this, j_para_spurious )
    class(Current_Density), intent(inout) :: this
    !> Spurious paramagnetic current density \mathbf{J}_{\rm spurious} obtained at \( t = 0 \)
    type(Current_Density_Field), intent(in) :: j_para_spurious

    this%spurious = j_para_spurious
  end subroutine

end module rttddft_CurrentDensity
