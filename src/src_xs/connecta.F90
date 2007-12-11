
! Copyright (C) 2002-2007 S. Sagmeister, J. K. Dewhurst, S. Sharma and 
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: connecta
! !INTERFACE:
subroutine connecta(cvec,nv,np,vvl,vpl,dv,dp)
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
!   {\tt dp}, respectively. Based upon the routine {\tt connect}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!   Modifications 2007 (Sagmeister)
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
integer iv,ip,j,c,spts,npi
real(8) st,sv,vl(3),vc(3),v1(3),v2(3)
! alloctable arrays
real(8), allocatable :: seg(:)
integer, allocatable :: idx(:), segpts(:)
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
allocate(seg(nv-1),idx(nv-1),segpts(nv-1)) !!!
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
! sort segments according to their length in descending order
call sortidx(nv-1,seg,idx)
npi=np-nv
segpts(:) = 0
do ip=1,npi
   iv=mod(ip,nv-1)
   if (iv.eq.0) iv=nv-1
   segpts(idx(nv-iv)) = segpts(idx(nv-iv)) + 1
end do
! loop over vertices
c=1
do iv=1,nv-1
   v1(:) = vvl(:,iv)
   v2(:) = vvl(:,iv+1)   
   vpl(:,c) = v1(:)
   dp(c) = dv(iv)
   c=c+1
   spts=segpts(iv)
   ! linear interplation along segment
   call linterplin(spts,v1,v2,vpl(1,c))
   dp(c:c+spts-1)=(/(dp(c-1)+seg(iv)*dble(j)/dble(spts+1),j=1,spts)/)
   c=c+spts
   if (iv.eq.nv-1) then
      vpl(:,c) = v2(:)
      dp(c) = dv(iv+1)
      c=c+1
   end if
end do
deallocate(seg,idx,segpts)
end subroutine
!EOC

subroutine linterplin(np,v1,v2,vintp)
  implicit none
  ! arguments
  integer, intent(in) :: np
  real(8), intent(in) :: v1(3), v2(3)
  real(8), intent(out) :: vintp(3,np)
  ! local variables
  real(8) :: np2,lam
  integer :: j
  np2=np+2
  do j=1,np
     lam=dble(j)/(np2-1)
     vintp(:,j)=v1(:)*(1-lam)+v2(:)*lam
  end do
end subroutine linterplin

