
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gengvec
! !INTERFACE:
subroutine gengvec
! !USES:
use modmain
! !DESCRIPTION:
!   Generates a set of ${\bf G}$-vectors used for the Fourier transform of the
!   charge density and potential and sorts them according to length. Integers
!   corresponding to the vectors in lattice coordinates are stored, as well as
!   the map from these integer coordinates to the ${\bf G}$-vector index. A map
!   from the ${\bf G}$-vector set to the standard FFT array structure is also
!   generated. Finally, the number of ${\bf G}$-vectors with magnitude less than
!   {\tt gmaxvr} is determined.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!   Increased number of G-vectors to ngrtot, July 2007 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ig,i1,i2,i3,j1,j2,j3,k
real(8) v(3),t1
! allocatable arrays
integer, allocatable :: idx(:)
integer, allocatable :: iar(:)
real(8), allocatable :: rar(:)
! allocate local arrays
allocate(idx(ngrtot))
allocate(iar(ngrtot))
allocate(rar(ngrtot))
! allocate global G-vector arrays
if (allocated(ivg)) deallocate(ivg)
allocate(ivg(3,ngrtot))
if (allocated(ivgig)) deallocate(ivgig)
allocate(ivgig(intgv(1,1):intgv(1,2),intgv(2,1):intgv(2,2), &
 intgv(3,1):intgv(3,2)))
if (allocated(igfft)) deallocate(igfft)
allocate(igfft(ngrtot))
if (allocated(vgc)) deallocate(vgc)
allocate(vgc(3,ngrtot))
if (allocated(gc)) deallocate(gc)
allocate(gc(ngrtot))
ig=0
do i1=intgv(1,1),intgv(1,2)
  do i2=intgv(2,1),intgv(2,2)
    do i3=intgv(3,1),intgv(3,2)
      v(:)=dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
      t1=v(1)**2+v(2)**2+v(3)**2
      ig=ig+1
! map from G-vector to (i1,i2,i3) index
      ivg(1,ig)=i1
      ivg(2,ig)=i2
      ivg(3,ig)=i3
! length of each G-vector
      gc(ig)=sqrt(t1)
    end do
  end do
end do
! sort by vector length
call sortidx(ngrtot,gc,idx)
! re-order arrays
do ig=1,ngrtot
  rar(ig)=gc(ig)
end do
do ig=1,ngrtot
  gc(ig)=rar(idx(ig))
end do
do k=1,3
  do ig=1,ngrtot
    iar(ig)=ivg(k,ig)
  end do
  do ig=1,ngrtot
    ivg(k,ig)=iar(idx(ig))
  end do
end do
ivgig(:,:,:)=0
do ig=1,ngrtot
  i1=ivg(1,ig)
  i2=ivg(2,ig)
  i3=ivg(3,ig)
! map from (i1,i2,i3) index to G-vector
  ivgig(i1,i2,i3)=ig
! assign G-vectors to global array
  vgc(:,ig)=dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
! Fourier transform index
  if (i1.ge.0) then
    j1=i1+1
  else
    j1=ngrid(1)+i1+1
  end if
  if (i2.ge.0) then
    j2=i2+1
  else
    j2=ngrid(2)+i2+1
  end if
  if (i3.ge.0) then
    j3=i3+1
  else
    j3=ngrid(3)+i3+1
  end if
  igfft(ig)=(j3-1)*ngrid(2)*ngrid(1)+(j2-1)*ngrid(1)+j1
end do
! find the number of vectors with G < gmaxvr
ngvec=1
do ig=1,ngrtot
  if (gc(ig).ge.gmaxvr) then
    ngvec=ig
    goto 10
  end if
end do
10 continue
deallocate(idx,iar,rar)
return
end subroutine
!EOC

