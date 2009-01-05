
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
!   Improved September 2007 (JKD)
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
integer iv,ip,ip0,ip1,n
real(8) vl(3),vc(3)
real(8) dt,f,t1
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
dt=0.d0
do iv=1,nv-1
  dv(iv)=dt
  vl(:)=vvl(:,iv+1)-vvl(:,iv)
  call r3mv(cvec,vl,vc)
  seg(iv)=sqrt(vc(1)**2+vc(2)**2+vc(3)**2)
  dt=dt+seg(iv)
end do
dv(nv)=dt
if (dt.lt.1.d-8) then
  do ip=1,np
    vpl(:,ip)=vvl(:,1)
    dp(ip)=0.d0
  end do
else
  do iv=1,nv-1
    t1=dble(np)*dv(iv)/dt
    ip0=nint(t1)+1
    if (ip0.lt.1) ip0=1
    t1=dble(np)*dv(iv+1)/dt
    ip1=nint(t1)
    if (ip1.gt.np) ip1=np
    n=ip1-ip0
    if (n.le.0) n=1
    do ip=ip0,ip1
      f=dble(ip-ip0)/dble(n)
      dp(ip)=f*seg(iv)+dv(iv)
      vpl(:,ip)=vvl(:,iv)*(1.d0-f)+vvl(:,iv+1)*f
    end do
  end do
end if
deallocate(seg)
return
end subroutine
!EOC
