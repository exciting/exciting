
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: connect
! !INTERFACE:
subroutine connect(cvec,nv,np,vvl,vpl,dv,dp)
! !INPUT/OUTPUT PARAMETERS:
!   cvec : matrix of (reciprocal) lattice vectors stored column-wise
!         (in,real(3,3))
!   nv   : number of vertices (in,integer)
!   np   : number of connecting points (in,integer)
!   vvl  : vertex vectors in lattice coordinates (in,real(3,nv))
!   vpl  : connecting point vectors in lattice coordinates (out,real(3,np))
!   dv   : cummulative distance to each vertex (out,real(nv))
!   dp   : cummulative distance to each connecting point (out,real(np))
! !DESCRIPTION:
!   Generates a set of points which interpolate between a given set of vertices.
!   Vertex points are supplied in lattice coordinates in the array {\tt vvl} and
!   converted to Cartesian coordinates with the matrix {\tt cvec}. Interpolating
!   points are stored in the array {\tt vpl}. The cummulative distances to the
!   vertices and points along the path are stored in arrays {\tt dv} and
!   {\tt dp}, respectively.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: cvec(3,3)
integer, intent(in) :: nv
integer, intent(in) :: np
real(8), intent(in) :: vvl(3,nv)
real(8), intent(out) :: vpl(3,np)
real(8), intent(out) :: dv(nv)
real(8), intent(out) :: dp(np)
! local variables
integer iv,ip
real(8) st,sv,vl(3),vc(3),f
! alloctable arrays
real(8), allocatable :: seg(:)
if (nv.lt.1) then
  write(*,*)
  write(*,'("Error(connect): nv < 1 : ",I8)') nv
  write(*,*)
  stop
end if
if (np.lt.nv) then
  write(*,*)
  write(*,'("Error(connect): np < nv : ",2I8)') np,nv
  write(*,*)
  stop
end if
if (np.eq.1) then
  vpl(:,1)=vvl(:,1)
  dv(1)=0.d0
  dp(1)=0.d0
  return
end if
allocate(seg(nv))
! find the total distance and the length of each segment
st=0.d0
do iv=1,nv-1
  dv(iv)=st
  vl(:)=vvl(:,iv+1)-vvl(:,iv)
  vc(:)=vl(1)*cvec(:,1)+vl(2)*cvec(:,2)+vl(3)*cvec(:,3)
  seg(iv)=sqrt(vc(1)**2+vc(2)**2+vc(3)**2)
  st=st+seg(iv)
end do
dv(nv)=st
! compute the interpolating points
do ip=1,np
  dp(ip)=st*dble(ip-1)/dble(np-1)
  sv=st
! determine the segment of the current point
  do iv=nv-1,1,-1
    sv=sv-seg(iv)
    if (dp(ip).gt.sv-1.d-7) then
      if (seg(iv).gt.1.d-7) then
        f=(dp(ip)-sv)/seg(iv)
      else
        f=0.d0
      end if
      vpl(:,ip)=vvl(:,iv)*(1.d0-f)+vvl(:,iv+1)*f
      goto 10
    end if
  end do
10 continue
end do
deallocate(seg)
return
end subroutine
!EOC
