
! Copyright (C) 2007 S. Sagmeister J. K. Dewhurst, S. Sharma and 
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevalsv0(vpl,evalsvp)
  use modmain
  use modtddft
  ! arguments
  real(8), intent(in) :: vpl(3)
  real(8), intent(out) :: evalsvp(nstsv)
  ! local variables
  real(8), allocatable :: vklt(:,:)

  ! copy varialbes of k+(q=0) to default variables
  allocate(vklt(3,nkpt))
  vklt(:,:)=vkl(:,:); vkl(:,:)=vkl0(:,:)

  ! call to getevalsv with changed (G+)k-point sets / matrix size
  call getevalsv(vpl,evalsvp)

  ! restore original variables
  vkl(:,:)=vklt(:,:)
  deallocate(vklt)

end subroutine getevalsv0
