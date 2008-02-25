
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symsci0(flag,scrnh0,scrnih0,scrnisym)
  implicit none
  ! arguments
  integer, intent(in) :: flag
  complex(8), intent(in) :: scrnh0(3)
  complex(8), intent(in) :: scrnih0(3)
  complex(8), intent(out) :: scrnisym
  select case(flag)
  case(0)
     ! Peter's original choice in his BSE implementation:
     ! average the screening diagonal tensor components and take
     ! inverse of this value
     scrnisym=1.d0/((sum(scrnh0))/3.d0)
  case(1)
     ! Treatment found in the BSE implementation of R. Laskowski, also
     ! mentioned in [M. Hybertsen, PRB 34, 5390 (1986), p5411]: average
     ! the diagonal tensor components of the inverse of the screening matrix
     scrnisym=sum(scrnih0)/3.d0
  case default
     write(*,*)
     write(*,'("Error(symsci0): not a valid choice for symmetrizing")')
     write(*,'(" the screened Coulomb interaction for q=0:",i8)') flag
     write(*,*)
     call terminate
  end select
end subroutine symsci0
