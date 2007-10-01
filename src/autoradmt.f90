
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: autoradmt
! !INTERFACE:
subroutine autoradmt
! !USES:
use modmain
! !DESCRIPTION:
!   Automatically determines the muffin-tin radii from the formula
!   $$ R_i\propto 1+\zeta|Z_i|^{1/3}, $$
!   where $Z_i$ is the atomic number of the $i$th species, $\zeta$ is a
!   user-supplied constant ($\sim 0.625$). The parameter $\zeta$ is stored in
!   {\tt rmtapm(1)} and the value which governs the distance between the
!   muffin-tins is stored in {\tt rmtapm(2)}. When {\tt rmtapm(2)} $=1$, the
!   closest muffin-tins will touch.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!   Changed the formula, September 2006 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is1,is2,ia1,ia2
integer i1,i2,i3
real(8) s,v1(3),v2(3),t1,t2
! external functions
real(8) r3dist
external r3dist
! initial muffin-tin radii
do is1=1,nspecies
  rmt(is1)=1.d0+rmtapm(1)*abs(spzn(is1))**(1.d0/3.d0)
end do
! determine scaling factor
s=1.d6
do i1=-1,1
  do i2=-1,1
    do i3=-1,1
      v1(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3)
      do is1=1,nspecies
        do ia1=1,natoms(is1)
          v2(:)=v1(:)+atposc(:,ia1,is1)
          do is2=1,nspecies
            t1=rmt(is1)+rmt(is2)
            do ia2=1,natoms(is2)
              if ((i1.ne.0).or.(i2.ne.0).or.(i3.ne.0).or.(is1.ne.is2).or. &
               (ia1.ne.ia2)) then
                t2=r3dist(v2,atposc(1,ia2,is2))
                s=min(s,t2/t1)
              end if
            end do
          end do
        end do
      end do
    end do
  end do
end do
s=s*rmtapm(2)
! scale all radii
do is1=1,nspecies
  t1=s*rmt(is1)*10000.d0
  t1=dble(int(t1))/10000.d0
  rmt(is1)=t1
end do
return
end subroutine
!EOC

