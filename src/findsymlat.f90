
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: findsymlat
subroutine findsymlat
! !USES:
use modmain
! !DESCRIPTION:
!   Finds the point group symmetries which leave the Bravais lattice invariant.
!   Let $A$ be the matrix consisting of the lattice vectors in columns, then
!   $$ g=A^{\rm T}A $$
!   is the metric tensor. Any $3\times 3$ matrix $S$ with elements $-1$, 0 or 1
!   is a point group symmetry of the lattice if $\det(S)$ is $-1$ or 1, and
!   $$ S^{\rm T}gS=g. $$
!   The first matrix in the set returned is the identity.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!   Removed arguments and simplified, April 2007 (JKD)
!EOP
!BOC
implicit none
! local variables
integer i,j,md,sym(3,3)
integer i11,i12,i13,i21,i22,i23,i31,i32,i33
real(8) s(3,3),g(3,3),sgs(3,3),v(3),t1
! external functions
integer i3mdet
real(8) r3taxi
external i3mdet,r3taxi
! determine metric tensor
call r3mtm(avec,avec,g)
! loop over all possible symmetry matrices
nsymlat=0
do i11=-1,1; do i12=-1,1; do i13=-1,1
do i21=-1,1; do i22=-1,1; do i23=-1,1
do i31=-1,1; do i32=-1,1; do i33=-1,1
  sym(1,1)=i11; sym(1,2)=i12; sym(1,3)=i13;
  sym(2,1)=i21; sym(2,2)=i22; sym(2,3)=i23;
  sym(3,1)=i31; sym(3,2)=i32; sym(3,3)=i33;
! determinant of matrix
  md=i3mdet(sym)
! matrix should be orthogonal
  if (abs(md).ne.1) goto 10
! check invariance of metric tensor
  s(:,:)=dble(sym(:,:))
  call r3mtm(s,g,sgs)
  call r3mm(sgs,s,sgs)
  do i=1,3
    do j=1,3
      if (abs(sgs(i,j)-g(i,j)).gt.epslat) goto 10
    end do
  end do
! check invariance of spin-spiral q-vector if required
  if (spinsprl) then
    call r3mv(s,vqlss,v)
    t1=r3taxi(vqlss,v)
    if (t1.gt.epslat) goto 10
  end if
  nsymlat=nsymlat+1
  if (nsymlat.gt.48) then
    write(*,*)
    write(*,'("Error(findsymlat): more than 48 symmetries found")')
    write(*,'(" (lattice vectors may be linearly dependent)")')
    write(*,*)
    stop
  end if
  symlat(:,:,nsymlat)=sym(:,:)
10 continue
end do; end do; end do
end do; end do; end do
end do; end do; end do
if (nsymlat.eq.0) then
  write(*,*)
  write(*,'("Error(findsymlat): no symmetries found")')
  write(*,*)
  stop
end if
! make the first symmetry the identity
do i=1,nsymlat
  if ((symlat(1,1,i).eq.1).and.(symlat(1,2,i).eq.0).and.(symlat(1,3,i).eq.0) &
 .and.(symlat(2,1,i).eq.0).and.(symlat(2,2,i).eq.1).and.(symlat(2,3,i).eq.0) &
 .and.(symlat(3,1,i).eq.0).and.(symlat(3,2,i).eq.0).and.(symlat(3,3,i).eq.1)) &
  then
    sym(:,:)=symlat(:,:,1)
    symlat(:,:,1)=symlat(:,:,i)
    symlat(:,:,i)=sym(:,:)
    goto 20
  end if
end do
20 continue
return
end subroutine
!EOC

