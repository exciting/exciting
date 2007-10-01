
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: findsymcrysctr
! !INTERFACE:
subroutine findsymcrysctr(eps,nspecies,natoms,nsymlat,symlat,ld,atposl,tpolar, &
 plrvl,ctrvl,nsymcrys,symcrys)
! !INPUT/OUTPUT PARAMETERS:
!   eps      : zero vector tolerance (in,real)
!   nspecies : number of species (in,integer)
!   natoms   : number atoms for each species (in,integer(nspecies))
!   nsymlat  : number of lattice point group symmetries (in,integer)
!   symlat   : lattice point group symmetries (in,integer(3,3,48))
!   ld       : leading dimension (in,integer)
!   atposl   : atomic positions in lattice coordinates (in,real(3,ld,nspecies))
!   tpolar   : .true. if plrvl are polar vectors, .false. if they are axial
!              vectors (in,logical)
!   plrvl    : polarisation vector at each atom site in lattice coordinates
!              (in,real(3,ld,nspecies))
!   ctrvl    : symmetry center in lattice coordinates (in,real(3))
!   nsymcrys : number of crystal point group symmetries found (out,integer)
!   symcrys  : crystal point group symmetries found (out,integer(3,3,48))
! !DESCRIPTION:
!   For a given symmetry center, this routine finds those point group symmetries
!   which leave the crystal structure, including atomic positions and
!   polarisation vectors at each site, invariant. If {\tt tpolar} is
!   {\tt .true.} then the polarisation vectors are taken to be polar, otherwise
!   they are axial and invariant under inversion. The routine requires the
!   symmetries which leave the Bravais lattice invariant. Vectors which are
!   within a distance $\epsilon$ of each other are considered to be equal.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
integer, intent(in) :: nspecies
integer, intent(in) :: natoms(nspecies)
integer, intent(in) :: nsymlat
integer, intent(in) :: symlat(3,3,48)
integer, intent(in) :: ld
real(8), intent(in) :: atposl(3,ld,nspecies)
logical, intent(in) :: tpolar
real(8), intent(in) :: plrvl(3,ld,nspecies)
real(8), intent(in) :: ctrvl(3)
integer, intent(out) :: nsymcrys
integer, intent(out) :: symcrys(3,3,48)
! local variables
integer isym,is,ia1,ia2,id(3),md
real(8) s(3,3)
! allocatable arrays
real(8), allocatable :: av1(:,:),av2(:,:)
real(8), allocatable :: pv(:,:)
! external functions
integer i3mdet
real(8) r3taxi
external i3mdet,r3taxi
if (nspecies.lt.0) then
  write(*,*)
  write(*,'("Error(findsymcrysctr): invalid nspecies : ",I8)') nspecies
  write(*,*)
  stop
end if
if ((nsymlat.le.0).or.(nsymlat.gt.48)) then
  write(*,*)
  write(*,'("Error(findsymcrysctr): invalid nsymlat : ",I8)') nsymlat
  write(*,*)
  stop
end if
allocate(av1(3,ld),av2(3,ld))
allocate(pv(3,ld))
nsymcrys=0
do isym=1,nsymlat
! determinant of symmetry matrix
  md=i3mdet(symlat(1,1,isym))
  s(:,:)=dble(symlat(:,:,isym))
  do is=1,nspecies
    do ia1=1,natoms(is)
      av1(:,ia1)=atposl(:,ia1,is)-ctrvl(:)
      call r3frac(eps,av1(1,ia1),id)
      call r3mv(s,av1(1,ia1),av2(1,ia1))
      call r3frac(eps,av2(1,ia1),id)
      call r3mv(s,plrvl(1,ia1,is),pv(1,ia1))
      if (.not.tpolar) pv(:,ia1)=dble(md)*pv(:,ia1)
    end do
    do ia1=1,natoms(is)
      do ia2=1,natoms(is)
        if ((r3taxi(av1(1,ia1),av2(1,ia2)).lt.eps).and. &
         (r3taxi(plrvl(1,ia1,is),pv(1,ia2)).lt.eps)) goto 10
      end do
      goto 20
10 continue
    end do
  end do
  nsymcrys=nsymcrys+1
  symcrys(:,:,nsymcrys)=symlat(:,:,isym)
20 continue
end do
deallocate(av1,av2,pv)
return
end subroutine
!EOC
