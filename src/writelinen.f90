
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writelinen
! !INTERFACE:
subroutine writelinen
! !USES:
use modmain
! !DESCRIPTION:
!   Writes the linearisation energies for all APW and local-orbital functions to
!   the file {\tt LINENGY.OUT}.
!
! !REVISION HISTORY:
!   Created February 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,l,io,ilo
open(50,file='LINENGY'//trim(filext),action='WRITE',form='FORMATTED')
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    write(50,'(" APW functions :")')
    do l=0,lmaxapw
      do io=1,apword(l,is)
        write(50,'("  l = ",I2,", order = ",I2," : ",G18.10)') l,io, &
         apwe(io,l,ias)
      end do
    end do
    write(50,'(" local-orbital functions :")')
    do ilo=1,nlorb(is)
      do io=1,lorbord(ilo,is)
        write(50,'("  l.o. = ",I2,", l = ",I2,", order = ",I2," : ",G18.10)') &
         ilo,lorbl(ilo,is),io,lorbe(io,ilo,ias)
      end do
    end do
  end do
end do
close(50)
return
end subroutine
!EOC
