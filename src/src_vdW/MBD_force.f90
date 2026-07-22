!> Driver for the MBD@rsSCS dispersion-force correction.
module mbd_force_module
  implicit none

  private
  public :: MBD_force

contains

!> Compute the MBD@rsSCS dispersion force and store it in the global
!> dispersion-force array.
subroutine MBD_force
  use mod_force, only: force_disp
  use modinput, only: input
  use modMBD, only: longfMBD
  use vdw_general_routines, only: set_default_vdW_parameters

  implicit none

  call set_default_vdW_parameters
  call longfMBD(input%groundstate%MBDparameters%beta, input%groundstate%MBDparameters%d, force_disp)

end subroutine MBD_force

end module mbd_force_module
