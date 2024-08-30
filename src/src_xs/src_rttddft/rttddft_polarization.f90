module rttddft_Polarization
  use rttddft_VectorField, only: Uniform_Vector_Field

  implicit none

  private

  !> Polarization vector \(\mathbf{P}\)
  type, public, extends(Uniform_Vector_Field) :: Polarization
  end type

end module