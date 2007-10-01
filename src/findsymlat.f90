
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: findsymlat
! !INTERFACE:
subroutine findsymlat(eps,avec,tpolar,plrvl,nsymlat,symlat)
! !INPUT/OUTPUT PARAMETERS:
!   eps     : zero vector tolerance (in,real)
!   avec    : lattice vectors stored column-wise (in,real(3,3))
!   tpolar  : .true. if plrvl is a polar vector, .false. if it is an axial
!             vector (in,logical)
!   plrvl   : polarisation vector (in,real(3))
!   nsymlat : number of lattice point group symmetries found (out,integer)
!   symlat  : lattice point group symmetries found (out,integer(3,3,48))
! !DESCRIPTION:
!   Finds the point group symmetries which leave the Bravais lattice and a
!   polarisation vector, ${\bf p}$, invariant. If {\tt tpolar} is {\tt .true.}
!   then ${\bf p}$ is assumed to be a polar vector, otherwise it is axial and
!   invariant under inversion. Let $A$ be the matrix consisting of lattice
!   vectors in columns, then
!   $$ g=A^{\rm T}A $$
!   is the metric tensor. Any $3\times 3$ matrix $S$ with elements $-1$, 0 or 1
!   is a symmetry of the lattice if the following conditions hold:
!   \begin{align*}
!    S^{\rm T}gS&=g, \\
!    \kappa S{\bf p}&={\bf p},
!   \end{align*}
!   where $\kappa=1$ if ${\bf p}$ is polar, and $\kappa=|S|$ if ${\bf p}$ is
!   axial. Vectors which are within a distance $\epsilon$ of each other are
!   considered to be equal. The first matrix in the set returned is the
!   identity.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(in) :: avec(3,3)
logical, intent(in) :: tpolar
real(8), intent(in) :: plrvl(3)
integer, intent(out) :: nsymlat
integer, intent(out) :: symlat(3,3,48)
! local variables
integer i,j,md
integer i11,i12,i13,i21,i22,i23,i31,i32,i33,sym(3,3)
equivalence (i11,sym(1,1)),(i12,sym(1,2)),(i13,sym(1,3)), &
            (i21,sym(2,1)),(i22,sym(2,2)),(i23,sym(2,3)), &
            (i31,sym(3,1)),(i32,sym(3,2)),(i33,sym(3,3))
real(8) s(3,3),g(3,3),stg(3,3),stgs(3,3),v(3)
! external functions
integer i3mdet
real(8) r3taxi
external i3mdet,r3taxi
! define metric tensor
call r3mtm(avec,avec,g)
! find symmetries which leave metric invariant
nsymlat=0
do i11=-1,1; do i12=-1,1; do i13=-1,1
do i21=-1,1; do i22=-1,1; do i23=-1,1
do i31=-1,1; do i32=-1,1; do i33=-1,1
! determinant of matrix
  md=i3mdet(sym)
! matrix should be unitary
  if (abs(md).ne.1) goto 10
! check plrvl invariance
  s(:,:)=dble(sym(:,:))
  call r3mv(s,plrvl,v)
  if (.not.tpolar) v(:)=dble(md)*v(:)
  if (r3taxi(plrvl,v).gt.eps) goto 10
! check metric tensor invariance
  call r3mtm(s,g,stg)
  call r3mm(stg,s,stgs)
  do i=1,3
    do j=1,3
      if (abs(stgs(i,j)-g(i,j)).gt.eps) goto 10
    end do
  end do
  if (nsymlat.ge.48) then
    write(*,*)
    write(*,'("Error(findsymlat): more than 48 symmetries!")')
    write(*,'(" (lattice vectors may be linearly dependent)")')
    write(*,*)
    stop
  end if
  nsymlat=nsymlat+1
  symlat(:,:,nsymlat)=sym(:,:)
10 continue
end do; end do; end do
end do; end do; end do
end do; end do; end do
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
