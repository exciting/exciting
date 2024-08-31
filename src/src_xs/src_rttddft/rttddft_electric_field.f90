module rttddft_electric_field
  use physical_constants, only: c
  use precision, only: dp
  use rttddft_VectorPotential, only: Vector_Potential_Field
  use rttddft_VectorField, only: Uniform_Vector_Field

  implicit none

  private

  public :: obtain_electric_field

  !> Electric Field \(\mathbf{E}\)
  type, public, extends(Uniform_Vector_Field) :: Electric_Field
  contains
    procedure :: obtain_electric_field
  end type

contains

!> Obtain the electric field as the time derivative of the vector potential
pure subroutine obtain_electric_field( this, dt, A_tot, A_tot_previous )
  class(Electric_Field), intent(inout) :: this
  !> time step
  real(dp), intent(in)  :: dt
  !> vector potential at time `t`
  type(Vector_Potential_Field), intent(in) :: A_tot
  !> vector potential at time `t-dt`
  type(Vector_Potential_Field), intent(in) :: A_tot_previous

  this%components = ( -1._dp / c / dt ) * ( A_tot%components - A_tot_previous%components )
end subroutine

end module