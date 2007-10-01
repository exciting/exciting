
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: genppts
! !INTERFACE:
subroutine genppts(reducep,ngridp,vploff,nppt,ipmap,ivp,vpl,vpc,wppt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   reducep : .true. if p-point set is to be reduced (in,logical)
!   ngridp  : p-point grid size (in,integer(3))
!   vploff  : offset of p-point grid in lattice coordinates (in,real(3))
!   nppt    : total number of p-points (out,integer)
!   ipmap   : map from integer grid to p-point index
!             (out,integer(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1))
!   ivp     : integer coordinates of the p-points
!             (out,integer(3,ngridp(1)*ngridp(2)*ngridp(3)))
!   vpl     : lattice coordinates of each p-point
!             (out,real(3,ngridp(1)*ngridp(2)*ngridp(3)))
!   vpc     : Cartesian coordinates of each p-point
!             (out,real(3,ngridp(1)*ngridp(2)*ngridp(3)))
!   wppt    : weights of each p-point (out,real(ngridp(1)*ngridp(2)*ngridp(3)))
! !DESCRIPTION:
!   This routine is used for generating $k$-point or $q$-point sets. Since these
!   are stored in global arrays, the points passed to this and other routines
!   are referred to as $p$-points. If {\tt reducep} is {\tt .true.} the set is
!   reduced with the spatial part of the crystal symmetries. In lattice
!   coordinates the vectors are given by
!   $$ {\bf p}=(\frac{i_1}{n_1},\frac{i_2}{n_2},\frac{i_3}{n_3})+
!    {\bf v}_{\rm off}, $$
!   where $i_j$ runs from 0 to $n_j-1$ and $0\le{\bf v}_{{\rm off};j}<1$ for
!   $j=1,2,3$. The $p$-point weights are stored in {\tt wppt} and the array
!   {\tt ipmap} contains the map from the integer coordinates to the reduced
!   index.
!
! !REVISION HISTORY:
!   Created August 2002 (JKD)
!   Updated April 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: reducep
integer, intent(in) :: ngridp(3)
real(8), intent(in) :: vploff(3)
integer, intent(out) :: nppt
integer, intent(out) :: ipmap(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1)
integer, intent(out) :: ivp(3,ngridp(1)*ngridp(2)*ngridp(3))
real(8), intent(out) :: vpl(3,ngridp(1)*ngridp(2)*ngridp(3))
real(8), intent(out) :: vpc(3,ngridp(1)*ngridp(2)*ngridp(3))
real(8), intent(out) :: wppt(ngridp(1)*ngridp(2)*ngridp(3))
! local variables
integer i1,i2,i3,ip,jp
integer isym,lspl,iv(3)
real(8) v1(3),v2(3)
real(8) s(3,3),t1,t2
! external functions
real(8) r3taxi
external r3taxi
if ((ngridp(1).le.0).or.(ngridp(2).le.0).or.(ngridp(3).le.0)) then
  write(*,*)
  write(*,'("Error(genppts): invalid ngridp : ",3I8)') ngridp
  write(*,*)
  stop
end if
if ((vploff(1).lt.0.d0).or.(vploff(1).gt.1.d0-epslat).or. &
    (vploff(2).lt.0.d0).or.(vploff(2).gt.1.d0-epslat).or. &
    (vploff(3).lt.0.d0).or.(vploff(3).gt.1.d0-epslat)) then
  write(*,*)
  write(*,'("Error(genppts): vploff not in [0,1) interval : ",3G18.10)') &
   vploff
  write(*,*)
  stop
end if
t1=1.d0/dble(ngridp(1)*ngridp(2)*ngridp(3))
ip=0
do i3=0,ngridp(3)-1
  v1(3)=(dble(i3)+vploff(3))/dble(ngridp(3))
  do i2=0,ngridp(2)-1
    v1(2)=(dble(i2)+vploff(2))/dble(ngridp(2))
    do i1=0,ngridp(1)-1
      v1(1)=(dble(i1)+vploff(1))/dble(ngridp(1))
      if (reducep) then
! determine if this point is equivalent to one already in the set
        do isym=1,nsymcrys
          lspl=lsplsymc(isym)
          s(:,:)=dble(symlat(:,:,lspl))
          call r3mtv(s,v1,v2)
          call r3frac(epslat,v2,iv)
          do jp=1,ip
            t2=r3taxi(vpl(1,jp),v2)
            if (t2.lt.epslat) then
! equivalent k-point found so add to current weight
              ipmap(i1,i2,i3)=jp
              wppt(jp)=wppt(jp)+t1
              goto 10
            end if
          end do
        end do
      end if
! add new point to set
      ip=ip+1
      ipmap(i1,i2,i3)=ip
      ivp(1,ip)=i1; ivp(2,ip)=i2; ivp(3,ip)=i3
      vpl(:,ip)=v1(:)
      wppt(ip)=t1
10 continue
    end do
  end do
end do
nppt=ip
! determine the Cartesian coordinates of the p-points
do ip=1,nppt
  call r3mv(bvec,vpl(1,ip),vpc(1,ip))
end do
return
end subroutine
!EOC

