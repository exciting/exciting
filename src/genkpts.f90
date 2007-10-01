
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: genkpts
! !INTERFACE:
subroutine genkpts(eps,nsym,sym,ngridk,bvec,vkloff,nkpt,ikmap,ivk,vkl,vkc,wkpt)
! !INPUT/OUTPUT PARAMETERS:
!   eps    : zero vector tolerance (in,real)
!   nsym   : number of symmetry matrices (in,integer)
!   sym    : symmetry matrices (in,integer(3,3,nsym))
!   ngridk : k-point grid size (in,integer(3))
!   bvec   : reciprocal lattice vectors (in,real(3,3))
!   vkloff : offset of k-point grid in lattice coordinates (in,real(3))
!   nkpt   : total number of k-points (out,integer)
!   ikmap  : map from grid to k-point set
!            (out,integer(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
!   ivk    : integer coordinates of the k-points
!            (out,integer(3,ngridk(1)*ngridk(2)*ngridk(3)))
!   vkl    : lattice coordinates of each k-point
!            (out,real(3,ngridk(1)*ngridk(2)*ngridk(3)))
!   vkc    : Cartesian coordinates of each k-point
!            (out,real(3,ngridk(1)*ngridk(2)*ngridk(3)))
!   wkpt   : weights of each k-point (out,real(ngridk(1)*ngridk(2)*ngridk(3)))
! !DESCRIPTION:
!   This routine is used for generating $k$-point or $q$-point sets. The set is
!   reduced with the symmetries in the array {\tt sym}. If no reduction is
!   required then {\tt nsym} should be set to 1. In lattice coordinates the
!   vectors are given by
!   $$ {\bf k}=(\frac{i_1}{n_1},\frac{i_2}{n_2},\frac{i_3}{n_3})+
!    {\bf v}_{\rm off}, $$
!   where $i_j$ runs from 0 to $n_j-1$ and $0\le{\bf v}_{{\rm off};j}<1$ for
!   $j=1,2,3$. The array {\tt ikmap} contains the map from the $k$-point integer
!   coordinates to the reduced index.
!
! !REVISION HISTORY:
!   Created August 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
integer, intent(in) :: nsym
integer, intent(in) :: sym(3,3,nsym)
integer, intent(in) :: ngridk(3)
real(8), intent(in) :: bvec(3,3)
real(8), intent(in) :: vkloff(3)
integer, intent(out) :: nkpt
integer, intent(out) :: ikmap(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1)
integer, intent(out) :: ivk(3,ngridk(1)*ngridk(2)*ngridk(3))
real(8), intent(out) :: vkl(3,ngridk(1)*ngridk(2)*ngridk(3))
real(8), intent(out) :: vkc(3,ngridk(1)*ngridk(2)*ngridk(3))
real(8), intent(out) :: wkpt(ngridk(1)*ngridk(2)*ngridk(3))
! local variables
integer i1,i2,i3,ik1,ik2
integer isym,id(3),nkptnr
real(8) t1,v1(3),s(3,3)
! allocatable arrays
logical, allocatable :: used(:)
integer, allocatable :: ivknr(:,:)
real(8), allocatable :: vklnr(:,:)
! external functions
real(8) r3taxi
external r3taxi
allocate(used(ngridk(1)*ngridk(2)*ngridk(3)))
allocate(ivknr(3,ngridk(1)*ngridk(2)*ngridk(3)))
allocate(vklnr(3,ngridk(1)*ngridk(2)*ngridk(3)))
if ((nsym.lt.1).or.(nsym.gt.48)) then
  write(*,*)
  write(*,'("Error(genkpts): invalid nsym : ",I8)') nsym
  write(*,*)
  stop
end if
if ((ngridk(1).le.0).or.(ngridk(2).le.0).or.(ngridk(3).le.0)) then
  write(*,*)
  write(*,'("Error(genkpts): invalid ngridk : ",3I8)') ngridk
  write(*,*)
  stop
end if
if ((vkloff(1).lt.0.d0).or.(vkloff(1).gt.1.d0-eps).or. &
    (vkloff(2).lt.0.d0).or.(vkloff(2).gt.1.d0-eps).or. &
    (vkloff(3).lt.0.d0).or.(vkloff(3).gt.1.d0-eps)) then
  write(*,*)
  write(*,'("Error(genkpts): vkloff not in [0,1) interval : ",3G18.10)') &
   vkloff
  write(*,*)
  stop
end if
ik1=0
do i3=0,ngridk(3)-1
  do i2=0,ngridk(2)-1
    do i1=0,ngridk(1)-1
      ik1=ik1+1
      ivknr(1,ik1)=i1
      ivknr(2,ik1)=i2
      ivknr(3,ik1)=i3
      vklnr(:,ik1)=(dble(ivknr(:,ik1))+vkloff(:))/dble(ngridk(:))
      ikmap(i1,i2,i3)=ik1
    end do
  end do
end do
nkptnr=ik1
t1=1.d0/dble(ngridk(1)*ngridk(2)*ngridk(3))
if (nsym.gt.1) then
  used(:)=.false.
  wkpt(:)=0.d0
  nkpt=0
  do ik1=1,nkptnr
    if (.not.used(ik1)) then
      nkpt=nkpt+1
      ivk(:,nkpt)=ivknr(:,ik1)
      vkl(:,nkpt)=vklnr(:,ik1)
      do isym=1,nsym
        s(:,:)=dble(sym(:,:,isym))
        call r3mtv(s,vkl(1,nkpt),v1)
        call r3frac(eps,v1,id)
        do ik2=ik1,nkptnr
          if (.not.used(ik2)) then
            if (r3taxi(v1,vklnr(1,ik2)).lt.eps) then
              ikmap(ivknr(1,ik2),ivknr(2,ik2),ivknr(3,ik2))=nkpt
              wkpt(nkpt)=wkpt(nkpt)+t1
              used(ik2)=.true.
              goto 10
            end if
          end if
        end do
10 continue
      end do
    end if
  end do
else
  nkpt=nkptnr
  ivk(:,:)=ivknr(:,:)
  vkl(:,:)=vklnr(:,:)
  wkpt(:)=t1
end if
! generate k-points in Cartesian coordinates
do ik1=1,nkpt
  vkc(:,ik1)=vkl(1,ik1)*bvec(:,1)+vkl(2,ik1)*bvec(:,2)+vkl(3,ik1)*bvec(:,3)
end do
deallocate(used,ivknr,vklnr)
return
end subroutine
!EOC
