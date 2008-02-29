
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rdmdedn(dedn)
! calculates derivative of total energy w.r.t. occupation numbers
use modmain
implicit none
! arguments
real(8), intent(out) :: dedn(nstsv,nkpt)
! allocatable
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: c(:,:)
integer ik,ist
allocate(evecsv(nstsv,nstsv),c(nstsv,nstsv))
dedn(:,:)=0.d0
do ik=1,nkpt
! get evecsv from a file
  call getevecsv(vkl(1,ik),evecsv)
! kinetic contribution
  call zgemm('C','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,dkdc(1,1,ik),nstsv, &
   zzero,c,nstsv)
  do ist=1,nstsv
! include Coulomb contribution
    dedn(ist,ik)=dedn(ist,ik)-(dble(c(ist,ist))+dble(vclmat(ist,ist,ik)))
  end do
end do
! add exchange correlation contribution
call rdmdexcdn(dedn)
deallocate(evecsv,c)
return
end subroutine

