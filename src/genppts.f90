
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: genppts
! !INTERFACE:
subroutine genppts(reducep,tfbz,ngridp,boxl,nppt,ipmap,ivp,vpl,vpc,wppt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   reducep : .true. if p-point set is to be reduced (in,logical)
!   tfbz    : .true. if vpl and vpc should be mapped to the first Brillouin
!             zone (in,logical)
!   ngridp  : p-point grid size (in,integer(3))
!   boxl    : corners of box containing p-points in lattice coordinates, the
!             first vector is the origin (in,real(3,4))
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
!   coordinates, the ${\bf p}$ vectors are given by
!   $$ {\bf p}=\left(\begin{matrix} & & \\
!     {\bf B}_2-{\bf B}_1 & {\bf B}_3-{\bf B}_1 & {\bf B}_4-{\bf B}_1 \\
!       & & \end{matrix}\right)
!     \left(\begin{matrix}i_1/n_1 \\ i_2/n_2 \\ i_3/n_3 \end{matrix}\right)
!     +{\bf B}_1 $$
!   where $i_j$ runs from 0 to $n_j-1$, and the ${\bf B}$ vectors define the
!   corners of a box with ${\bf B}_1$ as the origin. If {\tt tfbz} is
!   {\tt .true.} then the vectors {\tt vpl} (and {\tt vpc}) are mapped to the
!   first Brillouin zone. If {\tt tfbz} is {\tt .false.} and {\tt reducep} is
!   {\tt .true.} then the coordinates of {\tt vpl} are mapped to the $[0,1)$
!   interval. The $p$-point weights are stored in {\tt wppt} and the array
!   {\tt ipmap} contains the map from the integer coordinates to the reduced
!   index.
!
! !REVISION HISTORY:
!   Created August 2002 (JKD)
!   Updated April 2007 (JKD)
!   Added mapping to the first Brillouin zone, September 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: reducep
logical, intent(in) :: tfbz
integer, intent(in) :: ngridp(3)
real(8), intent(in) :: boxl(3,4)
integer, intent(out) :: nppt
integer, intent(out) :: ipmap(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1)
integer, intent(out) :: ivp(3,ngridp(1)*ngridp(2)*ngridp(3))
real(8), intent(out) :: vpl(3,ngridp(1)*ngridp(2)*ngridp(3))
real(8), intent(out) :: vpc(3,ngridp(1)*ngridp(2)*ngridp(3))
real(8), intent(out) :: wppt(ngridp(1)*ngridp(2)*ngridp(3))
! local variables
integer i1,i2,i3,ip,jp
integer isym,lspl,iv(3)
real(8) v1(3),v2(3),v3(3)
real(8) b(3,3),s(3,3),t1,t2
if ((ngridp(1).le.0).or.(ngridp(2).le.0).or.(ngridp(3).le.0)) then
  write(*,*)
  write(*,'("Error(genppts): invalid ngridp : ",3I8)') ngridp
  write(*,*)
  stop
end if
! box vector matrix
b(:,1)=boxl(:,2)-boxl(:,1)
b(:,2)=boxl(:,3)-boxl(:,1)
b(:,3)=boxl(:,4)-boxl(:,1)
t1=1.d0/dble(ngridp(1)*ngridp(2)*ngridp(3))
ip=0
do i3=0,ngridp(3)-1
  v1(3)=dble(i3)/dble(ngridp(3))
  do i2=0,ngridp(2)-1
    v1(2)=dble(i2)/dble(ngridp(2))
    do i1=0,ngridp(1)-1
      v1(1)=dble(i1)/dble(ngridp(1))
      call r3mv(b,v1,v2)
      v2(:)=v2(:)+boxl(:,1)
      if (reducep) then
        call r3frac(epslat,v2,iv)
! determine if this point is equivalent to one already in the set
        do isym=1,nsymcrys
          lspl=lsplsymc(isym)
          s(:,:)=dble(symlat(:,:,lspl))
          call r3mtv(s,v2,v3)
          call r3frac(epslat,v3,iv)
          do jp=1,ip
            t2=abs(vpl(1,jp)-v3(1))+abs(vpl(2,jp)-v3(2))+abs(vpl(3,jp)-v3(3))
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
      vpl(:,ip)=v2(:)
      wppt(ip)=t1
10 continue
    end do
  end do
end do
nppt=ip
do ip=1,nppt
! map vpl to the first Brillouin zone if required
  if (tfbz) call vecfbz(epslat,bvec,vpl(:,ip),iv)
! determine the Cartesian coordinates of the p-points
  call r3mv(bvec,vpl(:,ip),vpc(:,ip))
end do
return
end subroutine
!EOC

