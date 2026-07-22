!> Driver for the MBD@rsSCS dispersion-energy correction.
module mbd_energy_module
  implicit none

  private
  public :: MBD_energy

contains

!> Compute the MBD@rsSCS dispersion energy and write it to `MBD.OUT` for
!> property calculations.
subroutine MBD_energy
  use mod_energy, only: e_disp
  use modinput, only: input
  use modMBD, only: longMBD
  use vdw_general_routines, only: set_default_vdW_parameters

  implicit none

  integer :: mbd_unit

  call set_default_vdW_parameters
  call longMBD(input%groundstate%MBDparameters%beta, input%groundstate%MBDparameters%d, e_disp)

  if (associated(input%properties)) then
     if (associated(input%properties%MBD)) then
        open(newunit=mbd_unit, file="MBD.OUT")
        write(mbd_unit,'(F18.8)') e_disp
        close(mbd_unit)
     end if
  end if

end subroutine MBD_energy

end module mbd_energy_module
