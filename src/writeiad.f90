
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
!   Outputs the interatomic distances to the file {\tt IADIST.OUT}.
!
! !REVISION HISTORY:
!   Created May 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,js,ia,ja
integer i1,i2,i3
real(8) d,dmin,v(3)
! external functions
real(8) r3dist
external r3dist
open(50,file='IADIST'//trim(filext),action='WRITE',form='FORMATTED')
do is=1,nspecies
  do ia=1,natoms(is)
    write(50,*)
    write(50,'("Distance between is =",I4," (",A,"), ia =",I4," and")') is, &
     trim(spsymb(is)),ia
    do js=1,nspecies
      do ja=1,natoms(js)
        dmin=1.d8
        do i1=-1,1
          do i2=-1,1
            do i3=-1,1
              v(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3) &
               +atposc(:,ja,js)
              d=r3dist(atposc(1,ia,is),v)
              dmin=min(d,dmin)
            end do
          end do
        end do
        write(50,'(" is =",I4," (",A,"), ia =",I4," : ",G18.10)') js, &
         trim(spsymb(js)),ja,dmin
      end do
    end do
  end do
end do
close(50)
return
end subroutine
!EOC

