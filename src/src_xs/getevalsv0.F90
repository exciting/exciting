! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getevalsv0
! !INTERFACE:
subroutine getevalsv0(vpl, evalsvp)
! !USES:
  use mod_eigenvalue_occupancy, only: nstsv
  use mod_kpoint, only: nkptnr, vkl
  use mod_misc, only: task, filext
  use modxs, only: vkl0, filext0, usefilext0
! !INPUT/OUTPUT PARAMETERS:
! IN:
! real(8) :: vpl(3) ! k-point vector in lattice coordinates
! OUT:
! real(8) :: evalsvp(nstsv) ! Eigenvalues
! 
! !DESCTIPTION:
!   This routine reads the eigenvalues for a given k-vector
!   form the eigenvalue file assoziated with zero momentum 
!   transfer.
!
! !REVISION HISTORY:
!   Added to documentation scheme. 2016 (Aurich)
!
!EOP
!BOC
  implicit none

  real(8), intent(in) :: vpl(3)
  real(8), intent(out) :: evalsvp(nstsv)

  real(8), allocatable :: vklt(:, :)
  character(256) :: filextt

  ! Copy varialbes of k-grid to default variables
  ! and save the k' quantities in temporary arrays.
  allocate(vklt(3, nkptnr))
  vklt(:, :) = vkl(:, :)
  vkl(:, :) = vkl0(:, :)

  ! Save file extension for the k' quantities
  filextt = filext

  ! Change file extension to k file
  if(usefilext0) then 
    filext = filext0
  else
    call genfilextread(task)
  end if

  ! Call to getevalsv with changed (G+)k-point sets / matrix size
  call getevalsv(vpl, evalsvp)

  ! Restore original k+q variables
  vkl(:, :) = vklt(:, :)
  filext = filextt

  deallocate(vklt)
end subroutine getevalsv0
!EOC
