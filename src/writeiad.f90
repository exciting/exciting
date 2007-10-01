
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeiad
! !INTERFACE:
subroutine writeiad
! !USES:
use modmain
! !DESCRIPTION:
!   Outputs the interatomic distances to the file {\tt IAD.OUT}.
!
! !REVISION HISTORY:
!   Created May 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is1,is2,ia1,ia2
integer i1,i2,i3
real(8) d,dmin,v(3)
! external functions
real(8) r3dist
external r3dist
open(50,file='IAD'//trim(filext),action='WRITE',form='FORMATTED')
do is1=1,nspecies
  do ia1=1,natoms(is1)
    write(50,*)
    write(50,'("Distance between is =",I4,", ia =",I4," and")') is1,ia1
    do is2=1,nspecies
      do ia2=1,natoms(is2)
        dmin=1.d8
        do i1=-1,1
          do i2=-1,1
            do i3=-1,1
              v(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3) &
               +atposc(:,ia2,is2)
              d=r3dist(atposc(1,ia1,is1),v)
              dmin=min(d,dmin)
            end do
          end do
        end do
        write(50,'(" is =",I4,", ia =",I4," : ",G18.10)') is2,ia2,dmin
      end do
    end do
  end do
end do
close(50)
return
end subroutine
!EOC

