
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: checkmt
! !INTERFACE:
subroutine checkmt
! !USES:
use modmain
! !DESCRIPTION:
!   Checks for overlapping muffin-tins. If any muffin-tins are found to
!   intersect the program is terminated with error.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ia,is,ja,js
integer i1,i2,i3
real(8) v1(3),v2(3),v3(3)
real(8) t1,t2
do i1=-1,1
  do i2=-1,1
    do i3=-1,1
      v1(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3)
      do is=1,nspecies
        do ia=1,natoms(is)
          v2(:)=v1(:)+atposc(:,ia,is)
          do js=1,nspecies
            t1=(rmt(is)+rmt(js))**2
            do ja=1,natoms(js)
              if ((i1.ne.0).or.(i2.ne.0).or.(i3.ne.0).or.(is.ne.js).or. &
               (ia.ne.ja)) then
                v3(:)=v2(:)-atposc(:,ja,js)
                t2=v3(1)**2+v3(2)**2+v3(3)**2
                if (t1.gt.t2) then
                  write(*,*)
                  write(*,'("Error(checkmt): muffin-tin tins overlap for")')
                  write(*,'("     species ",I4," atom ",I4)') is,ia
                  write(*,'(" and species ",I4," atom ",I4)') js,ja
                  write(*,*)
                  write(*,'("Sum of muffin-tin radii : ",G18.10)') sqrt(t1)
                  write(*,'("Distance between atoms  : ",G18.10)') sqrt(t2)
                  write(*,*)
                  stop
                end if
              end if
            end do
          end do
        end do
      end do
    end do
  end do
end do
return
end subroutine
!EOC
