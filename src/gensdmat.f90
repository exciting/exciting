
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gensdmat
! !INTERFACE:
subroutine gensdmat(evecsv,sdmat)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   sdmat  : spin density matrices (out,complex(nspinor,nspinor,nstsv))
! !DESCRIPTION:
!   Computes the spin density matrices for a set of second-variational states.
!
! !REVISION HISTORY:
!   Created September 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: sdmat(nspinor,nspinor,nstsv)
! local variables
integer ispn,jspn,ist,j
complex(8) zt1,zt2
sdmat(:,:,:)=0.d0
do j=1,nstsv
  do ispn=1,nspinor
    do jspn=1,nspinor
      do ist=1,nstfv
        zt1=evecsv(ist+nstfv*(ispn-1),j)
        zt2=evecsv(ist+nstfv*(jspn-1),j)
        sdmat(ispn,jspn,j)=sdmat(ispn,jspn,j)+zt1*conjg(zt2)
      end do
    end do
  end do
end do
return
end subroutine
!EOC

