!> Module that deals with the Current Density in RT-TDDFT calculations
module rttddft_CurrentDensity
  use exciting_mpi, only: mpiinfo, xmpi_allreduce
  use mod_lattice, only: Omega
  use physical_constants, only: c
  use precision, only: i32, dp
  use rttddft_VectorField, only: Uniform_Vector_Field
  use rttddft_VectorPotential, only: Vector_Potential_Field
  use xlapack, only: dot_multiply, hermitian_matrix_multiply
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
  !> \[ \mathbf{J}_{ind}(t) = - \frac{N_{val} \mathbf{A}_{tot}(t)}{\Omega c} \]
  !> \(N_{val}\) is the number of valence electrons, \(c\) is the light speed, and
  !> \(\Omega\) is the unit cell volume
  pure subroutine Current_Density_evaluate_diamagnetic( this, Nel_per_volume, a_tot )
    class(Current_Density), intent(inout) :: this
    !> Number of valence electrons per volume
    real(dp), intent(in) :: Nel_per_volume
    !> Vector potential
    class(Vector_Potential_Field), intent(in) :: a_tot

    this%diamagnetic%components = ( -Nel_per_volume / c )*a_tot%components 
  end subroutine

  !> Calculate the paramagnetic part of the current density at time \( t \) as:
  !> \[
  !>    \mathbf{J}(t) = \frac{\mathrm{i}}{\Omega} \sum_{j\mathbf{k}}
  !>      w_{\mathbf{k}}f_{j\mathbf{k}} \left\langle \psi_{j\mathbf{k}}(t) \big|
  !>      \nabla \big|\psi_{j\mathbf{k}}(t)\right\rangle - \mathbf{J}_{\rm spurious}
  !>  \]
  !> where \( \Omega \) is the unit cell volume , \( w_{\mathbf{k}} \) is the 
  !> \( \mathbf{k} \)-point weight, \( f_{j\mathbf{k}} \) is the occupation number of the
  !> corresponding KS state, and \mathbf{J}_{\rm spurious} is the spurious current 
  !> attributed to numerical inaccuracy.
  subroutine Current_Density_evaluate_paramagnetic( this, psi, p_mat, mpi_env )
    class(Current_Density), intent(inout) :: this
    !> Basis-expansion coefficients of the KS-wavefunctions at time \( t \)
    class(wavefunction_set), intent(in) :: psi
    !> Momentum matrix elements
    complex(dp), intent(in) :: p_mat(:, :, :, :)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    integer(i32) :: ik, ist, j, n_states, n_basis, first_active, shift
    real(dp) :: aux(n_cartesian)
    real(dp), allocatable :: acc(:)
    complex(dp), allocatable :: draft(:, :)
    real(dp), parameter :: tol_default = 1e-6_dp

    first_active = psi%first_active()
    n_states = psi%n_active()
    n_basis = psi%n_basis()
    allocate( draft(n_basis, n_states), acc(n_states) )
    aux = 0._dp

    shift = 1-psi%first_kpt()
    !$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(ik, j, ist, draft, acc), REDUCTION(+:aux), &
    !$OMP& SHARED(n_states, psi, p_mat, first_active, shift)
    do ik = psi%first_kpt(), psi%last_kpt()
      ! For the x, y, and z components ...
      do j = 1, n_cartesian
        call hermitian_matrix_multiply( p_mat(:, :, j, ik+shift), psi%active(:, :, ik), draft, 'U', 'L', tol_default )
        do ist = 1, n_states
          acc(ist) = real( dot_multiply( psi%active(:, ist, ik), draft(:, ist), conjg_a=.true. ), dp )
        end do
        aux(j) = aux(j) - dot_multiply( psi%occupations(first_active: first_active + n_states - 1, ik), acc )*psi%kset%wkpt(ik)
      end do
    end do
    !$OMP END PARALLEL DO
    
    this%paramagnetic%components = aux / Omega
    call xmpi_allreduce( this%paramagnetic%components, mpi_env )
    call this%paramagnetic%add_vector( -this%spurious%components )
  end subroutine

  !> Set the spurious component to a use-provided value.
  pure subroutine Current_Density_set_spurious( this, j_para_spurious )
    class(Current_Density), intent(inout) :: this
    !> Spurious paramagnetic current density \mathbf{J}_{\rm spurious} obtained at \( t = 0 \)
    type(Current_Density_Field), intent(in) :: j_para_spurious

    this%spurious = j_para_spurious
  end subroutine

end module rttddft_CurrentDensity
