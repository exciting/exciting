
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rdmdedc(ik,dedc)
! calculate the derivative of total energy w.r.t. evecsv
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(out) :: dedc(nstsv,nstsv)
! local variables
integer ist
! allocatable arrays
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: c(:,:)
! allocate local arrays
allocate(evecsv(nstsv,nstsv))
allocate(c(nstsv,nstsv))
! get the eigenvectors from file
call getevecsv(vkl(1,ik),evecsv)
! kinetic and Coulomb potential contribution
call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,vclmat(1,1,ik),nstsv, &
 zzero,c,nstsv)
do ist=1,nstsv
  dedc(:,ist)=occsv(ist,ik)*(dkdc(:,ist,ik)+c(:,ist))
end do
! exchange-correlation contribution
call rdmdexcdc(ik,evecsv,dedc)
deallocate(evecsv,c)
return
end subroutine

