! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! REVISION HISTORY:
! Created February 2021 (Ronaldo Rodrigues Pela)

!> Module to perform consistency checks and throw an error or raise a warning
!> if something is not consistent
module errors_warnings
   use modmpi, only: terminate_mpi_env
   implicit none
   private
   public :: terminate_if_false
   public :: terminate_if_true

! Note this isn't overloaded as the MPI environment should get initialised and
! passed around even for serial calculations
contains

!> Terminate exciting if consistent is false
subroutine terminate_if_false( mpiglobal, consistent, error_message )
   use modmpi, only: mpiinfo
   !> mpi environment
   type(mpiinfo), intent(inout) :: mpiglobal
   !> logical condition to be tested
   logical, intent(in) :: consistent
   !> Error message to be printed out
   character(len=*), intent(in) :: error_message

   if ( .not. consistent ) then
      call terminate_mpi_env( mpiglobal, error_message )
   endif
end subroutine

!> Terminate exciting if consistent is true
subroutine terminate_if_true( mpiglobal, consistent, error_message )
   use modmpi, only: mpiinfo
   !> mpi environment
   type(mpiinfo), intent(inout) :: mpiglobal
   !> logical condition to be tested
   logical, intent(in) :: consistent
   !> Error message to be printed out
   character(len=*), intent(in) :: error_message

   if ( consistent ) then
      call terminate_mpi_env( mpiglobal, error_message )
   endif
end subroutine
end module errors_warnings
