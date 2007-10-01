
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: findprim
! !INTERFACE:
subroutine findprim(eps,avec,nspecies,natoms,ld,atposl,plrvc)
! !INPUT/OUTPUT PARAMETERS:
!   eps      : zero vector tolerance (in,real)
!   avec     : lattice vectors stored column-wise (inout,real(3,3))
!   nspecies : number of species (in,integer)
!   natoms   : number atoms for each species (inout,integer(nspecies))
!   ld       : leading dimension (in,integer)
!   atposl   : atomic positions in lattice coordinates
!              (inout,real(3,ld,nspecies))
!   plrvc    : polarisation vector at each atom site in lattice coordinates
!              (inout,real(3,ld,nspecies))
! !DESCRIPTION:
!   Given an input unit cell, this routine finds the smallest primitive cell
!   which produces the same crystal structure. This is done by searching
!   through all the all the vectors which connect atomic positions and finding
!   those which leave the crystal invariant. Of these, the three shortest which
!   produce a non-zero unit cell volume are returned.
!
! !REVISION HISTORY:
!   Created July 2005 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(inout) :: avec(3,3)
integer, intent(in) :: nspecies
integer, intent(inout) :: natoms(nspecies)
integer, intent(in) :: ld
real(8), intent(inout) :: atposl(3,ld,nspecies)
real(8), intent(inout) :: plrvc(3,ld,nspecies)
! local variables
integer is1,is2,ia1,ia2,ia3,ia4
integer i1,i2,i3,j1,j2,j3,k,n
integer np,npmax,id(3)
real(8) eps2,smin,omega,t1
real(8) v1(3),v2(3),v3(3),v4(3)
real(8) a(3,3),ai(3,3),b(3),c(3,3)
! allocatable arrays
real(8), allocatable :: apl(:,:)
real(8), allocatable :: avl(:,:)
real(8), allocatable :: avc(:,:)
real(8), allocatable :: ac(:)
real(8), allocatable :: pvc(:,:)
eps2=eps**2
! find species with minimum number of atoms
! also find maximum atoms over all species
n=0
is1=1
do is2=1,nspecies
  n=max(n,natoms(is2))
  if (natoms(is2).lt.natoms(is1)) is1=is2
  do ia1=1,natoms(is2)
! map atomic positions to [0,1)
    call r3frac(eps,atposl(1,ia1,is2),id)
  end do
end do
allocate(apl(3,n))
n=27*(natoms(is1)**2)
allocate(avl(3,n))
allocate(avc(3,n))
allocate(ac(n))
allocate(pvc(3,n))
! generate set of possible lattice vectors
n=0
do ia1=1,natoms(is1)
  do ia2=1,natoms(is1)
    v1(:)=atposl(:,ia1,is1)-atposl(:,ia2,is1)
    do i1=-1,1
      v2(1)=v1(1)+dble(i1)
      do i2=-1,1
        v2(2)=v1(2)+dble(i2)
        do i3=-1,1
          v2(3)=v1(3)+dble(i3)
! check if v2 is a lattice vector
          if ((v2(1)**2+v2(2)**2+v2(3)**2).gt.eps2) then
            do is2=1,nspecies
              do ia3=1,natoms(is2)
                v3(:)=atposl(:,ia3,is2)+v2(:)
                call r3frac(eps,v3,id)
                do ia4=1,natoms(is2)
! check positions are the same modulo a primitive translation
                  v4(:)=v3(:)-atposl(:,ia4,is2)
                  if ((v4(1)**2+v4(2)**2+v4(3)**2).lt.eps2) then
! check polarisation vectors are the same
                    v4(:)=plrvc(:,ia3,is2)-plrvc(:,ia4,is2)
                    if ((v4(1)**2+v4(2)**2+v4(3)**2).lt.eps2) goto 10
                  end if
                end do
                goto 20
10 continue
              end do
            end do
! check if lattice vector is already in array
            do j1=1,n
              v3(:)=v2(:)-avl(:,j1)
              if ((v3(1)**2+v3(2)**2+v3(3)**2).lt.eps2) goto 20
            end do
! add lattice vector to array
            n=n+1
            avl(:,n)=v2(:)
20 continue
          end if
        end do
      end do
    end do
  end do
end do
if (n.eq.0) then
  write(*,*)
  write(*,'("Error(findprim): no lattice vectors found")')
  write(*,*)
  stop
end if
! find the Cartesian coordinates and lengths of all lattice vectors
do i1=1,n
  call r3mv(avec,avl(1,i1),avc(1,i1))
  ac(i1)=sqrt(avc(1,i1)**2+avc(2,i1)**2+avc(3,i1)**2)
end do
! find the unit cell with smallest sum of lattice vector lengths and maximum
! number of positive coordinates
smin=1.d10
npmax=-1
j1=0; j2=0; j3=0
do i1=1,n
  do i2=1,n
    if (i1.ne.i2) then
      call r3cross(avc(1,i1),avc(1,i2),b)
      do i3=1,n
        if ((i1.ne.i3).and.(i2.ne.i3)) then
! unit cell volume
          omega=abs(avc(1,i3)*b(1)+avc(2,i3)*b(2)+avc(3,i3)*b(3))
          if (omega.gt.eps) then
! sum of lengths
            t1=ac(i1)+ac(i2)+ac(i3)
            if (t1.lt.smin+eps) then
! count the number of positive coordinates
              np=0
              do k=1,3
                if (avc(k,i1).gt.-eps) np=np+1
                if (avc(k,i2).gt.-eps) np=np+1
                if (avc(k,i3).gt.-eps) np=np+1
              end do
! reject if cell size is the same and number of positive coordinates is smaller
              if ((abs(t1-smin).lt.eps).and.(np.lt.npmax)) goto 30
              j1=i1
              j2=i2
              j3=i3
              smin=t1
              npmax=np
            end if
          end if
        end if
30 continue
      end do
    end if
  end do
end do
if (npmax.eq.-1) then
  write(*,*)
  write(*,'("Error(findprim): could not determine primitive cell")')
  write(*,*)
  stop
end if
a(:,1)=avc(:,j1)
a(:,2)=avc(:,j2)
a(:,3)=avc(:,j3)
! transformation matrix from old to new lattice coordinates
call r3minv(a,ai)
call r3mm(ai,avec,c)
! remove redundant atoms
do is1=1,nspecies
  n=0
  do ia1=1,natoms(is1)
    call r3mv(c,atposl(1,ia1,is1),v1)
    call r3frac(eps,v1,id)
    do ia2=1,n
      v2(:)=v1(:)-apl(:,ia2)
      if ((v2(1)**2+v2(2)**2+v2(3)**2).lt.eps2) goto 40
    end do
    n=n+1
    apl(:,n)=v1(:)
    pvc(:,n)=plrvc(:,ia1,is1)
40 continue
  end do
  atposl(:,1:n,is1)=apl(:,1:n)
  plrvc(:,1:n,is1)=pvc(:,1:n)
  natoms(is1)=n
end do
avec(:,:)=a(:,:)
deallocate(apl,avl,avc,ac,pvc)
return
end subroutine
!EOC

