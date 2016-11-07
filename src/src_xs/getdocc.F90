! Copyright(C) 2010 W. Olovsson, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getdocc
! !INTERFACE:
subroutine getdocc(iq, ik, ikq, l1, u1, l2, u2, docc)
! !USES:
  use mod_eigenvalue_occupancy, only: nstsv
  use mod_kpoint, only: vkl
  use modxs, only: vkl0
! !DESCRIPTION:
!   Calculates occupation number differences for 
!   $f_{o_l k_i} - f_{u_m k_j}$. Currently the {\tt iq}
!   input argument does nothing.
!   WARNING: {\tt xssave0} has to be called in advance.
!
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC

  implicit none

  ! Arguments
  integer, intent(in) :: iq, ik, ikq, l1, u1, l2, u2
  real(8), intent(out) :: docc(u1-l1+1, u2-l2+1)

  ! Local variables
  integer :: ist, jst, iqt
  real(8), allocatable :: o0(:), o(:)

  iqt = iq

  ! Allocate occupancy arrays for k and k+q
  allocate(o0(nstsv), o(nstsv))

  ! Eigenvalues and occupancies for k-point
  ! Reads from OCCSV_QMT000.OUT
  call getoccsv0(vkl0(1:3, ik), o0)

  ! Eigenvalues and occupancies for k+q-point
  ! Reads form the OCCSV* file with the current file extension.
  ! So for no momentum transfer it should also be set to 
  ! OCCSV_QMT000.OUT.
  call getoccsv(vkl(1:3, ikq), o)

  ! Loop over band ranges    
  do ist = l1, u1
    do jst = l2, u2
      docc(ist-l1+1, jst-l2+1) = o0(ist) - o(jst)
    end do
  end do

  deallocate(o0, o)
end subroutine getdocc
!EOC
