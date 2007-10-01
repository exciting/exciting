
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: spinchar
! !INTERFACE:
subroutine spinchar(ik,evecsv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point number (in,integer)
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Computes the spin character of a second-variational eigenstate. These are
!   given by
!   $$ \xi^{q}(\sigma)=\sum_p\left|C_{\sigma;p}^{q}\right|^2, $$
!   where $C_{\sigma;p}^{q}$ are the coefficients of the $q$th
!   second-variational state in the basis of first-variational states labelled
!   with $p$. The results are stored in the global array {\tt spnchr}.
!
! !REVISION HISTORY:
!   Created December 2005 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
integer ispn,ist,i,j
real(8) sum
if ((ik.le.0).or.(ik.gt.nkpt)) then
  write(*,*)
  write(*,'("Error(spinchar): ik out of range : ",I8)') ik
  write(*,*)
  stop
end if
if (spinpol) then
  do j=1,nstsv
    i=0
    do ispn=1,nspinor
      sum=0.d0
      do ist=1,nstfv
        i=i+1
        sum=sum+dble(evecsv(i,j))**2+aimag(evecsv(i,j))**2
      end do
      spnchr(ispn,j,ik)=sum
    end do
  end do
else
  spnchr(1,:,ik)=1.d0
end if
return
end subroutine
!EOC

