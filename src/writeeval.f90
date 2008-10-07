
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeeval
! !INTERFACE:
subroutine writeeval
! !USES:
use modmain
! !DESCRIPTION:
!   Outputs the second-variational eigenvalues and occupation numbers to the
!   file {\tt EIGVAL.OUT}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,is,ia,ias
! write out the valence eigenvalues
open(50,file='EIGVAL'//trim(filext),action='WRITE',form='FORMATTED')
write(50,'(I6," : nkpt")') nkpt
write(50,'(I6," : nstsv")') nstsv
do ik=1,nkpt
  write(50,*)
  write(50,'(I6,3G18.10," : k-point, vkl")') ik,vkl(:,ik)
  write(50,'(" (state, eigenvalue and occupancy below)")')
  do ist=1,nstsv
    write(50,'(I6,2G18.10)') ist,evalsv(ist,ik),occsv(ist,ik)
  end do
  write(50,*)
end do
close(50)
! write out the core eigenvalues
open(50,file='EVALCORE'//trim(filext),action='WRITE',form='FORMATTED')
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    do ist=1,spnst(is)
      if (spcore(ist,is)) then
        write(50,'(" n = ",I2,", l = ",I2,", k = ",I2," : ",G18.10)') &
         spn(ist,is),spl(ist,is),spk(ist,is),evalcr(ist,ias)
      end if
    end do
  end do
end do
close(50)
return
end subroutine
!EOC
