
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
integer is,js,ia,ja,i1,i2,i3
real(8) s,v1(3),v2(3),t1,t2,t3
! external functions
real(8) r3dist
external r3dist
! initial muffin-tin radii
do is=1,nspecies
  rmt(is)=1.d0+rmtapm(1)*abs(spzn(is))**(1.d0/3.d0)
end do
! determine scaling factor
s=1.d8
do i1=-1,1
  do i2=-1,1
    do i3=-1,1
      v1(:)=dble(i1)*avec(:,1)+dble(i2)*avec(:,2)+dble(i3)*avec(:,3)
      do is=1,nspecies
        do ia=1,natoms(is)
          v2(:)=v1(:)+atposc(:,ia,is)
          do js=1,nspecies
            t1=1.d0/(rmt(is)+rmt(js))
            do ja=1,natoms(js)
              if ((i1.ne.0).or.(i2.ne.0).or.(i3.ne.0).or.(is.ne.js).or. &
               (ia.ne.ja)) then
                t2=r3dist(v2,atposc(1,ja,js))
                t3=t1*t2
                if (t3.lt.s) s=t3
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
do is=1,nspecies
! limit number of decimal digits
  t1=s*rmt(is)*10000.d0
  t1=dble(int(t1))/10000.d0
  rmt(is)=t1
end do
return
end subroutine
!EOC

