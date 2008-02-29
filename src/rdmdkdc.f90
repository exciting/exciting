
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rdmdkdc
! calculates the derivative of kinetic energy w.r.t. evecsv
use modmain
implicit none
! allocatable arrays
complex(8), allocatable :: evecsv(:,:)
integer ik
allocate(evecsv(nstsv,nstsv))
do ik=1,nkpt
  call getevecsv(vkl(1,ik),evecsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,kinmatc(1,1,ik),nstsv,evecsv, &
   nstsv,zzero,dkdc(1,1,ik),nstsv)
end do
deallocate(evecsv)
return
end subroutine

