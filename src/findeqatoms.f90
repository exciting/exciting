
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: findeqatoms
! !INTERFACE:
subroutine findeqatoms(eps,nspecies,natoms,nsymlat,symlat,ld1,atposl,tpolar, &
 plrvl,ld2,nsymeqat,symeqat,tvleqat)
! !INPUT/OUTPUT PARAMETERS:
!   eps      : zero vector tolerance (in,real)
!   nspecies : number of species (in,integer)
!   natoms   : number atoms for each species (in,integer(nspecies))
!   nsymlat  : number of lattice point group symmetries (in,integer)
!   symlat   : lattice point group symmetries (in,integer(3,3,48))
!   ld1      : leading dimension 1 (in,integer)
!   atposl   : atomic positions in lattice coordinates (in,real(3,ld1,nspecies))
!   tpolar   : .true. if plrvl are polar vectors, .false. if they are axial
!              vectors (in,logical)
!   plrvl    : polarisation vector at each atom site in lattice coordinates
!              (in,real(3,ld1,nspecies))
!   ld2      : leading dimension 2 (in,integer)
!   nsymeqat : number of symmetries which map the equivalent atoms
!              (out,integer(ld2,ld2,nspecies))
!   symeqat  : symmetries in symlat which map the equivalent atoms
!              (out,integer(48,ld2,ld2,nspecies))
!   tvleqat  : translation vectors in lattice coordinates corresponding to each
!              symeqat applied about origin (out,real(3,48,ld2,ld2,nspecies))
! !DESCRIPTION:
!   This routine finds equivalent atoms and related symmetries. Let
!   $\{({\bf r}_i,{\bf p}_i):i=1,2,\ldots\}$ be the set of all atomic positions
!   and polarisation vectors in the crystal. If {\tt tpolar} is {\tt .true.}
!   then the polarisation vectors are taken to be polar, otherwise they are
!   axial and invariant under inversion. Atom $\beta$ is equivalent to atom
!   $\alpha$ if they are of the same species,
!   $$ ({\bf r}_{\alpha},{\bf p}_{\alpha})=(S{\bf r}_{\beta}+{\bf t},
!    \kappa S{\bf p}_{\beta}) $$
!   and
!   $$ \{({\bf r}_i,{\bf p}_i):i=1,2,\ldots\}=\{(S{\bf r}_i+{\bf t},
!    \kappa S{\bf p}_i):i=1,2,\ldots\}, $$
!   for some lattice point symmetry $S$ and translation ${\bf t}$, where
!   $\kappa=1$ if ${\bf p}$ is polar, and $\kappa=|S|$ if ${\bf p}$ is axial.
!   Vectors which are within a distance $\epsilon$ of each other are considered
!   to be equal.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
integer, intent(in) :: nspecies
integer, intent(in) :: natoms(nspecies)
integer, intent(in) :: nsymlat
integer, intent(in) :: symlat(3,3,48)
integer, intent(in) :: ld1
real(8), intent(in) :: atposl(3,ld1,nspecies)
logical, intent(in) :: tpolar
real(8), intent(in) :: plrvl(3,ld1,nspecies)
integer, intent(in) :: ld2
integer, intent(out) :: nsymeqat(ld2,ld2,nspecies)
integer, intent(out) :: symeqat(48,ld2,ld2,nspecies)
real(8), intent(out) :: tvleqat(3,48,ld2,ld2,nspecies)
! local variables
integer is1,is2,ia1,ia2,ia3,ia4
integer isym,i,md,id(3)
real(8) s(3,3),v1(3),v2(3)
! allocatable arrays
logical, allocatable :: used(:)
real(8), allocatable :: av1(:,:,:)
real(8), allocatable :: av2(:,:,:)
! external functions
integer i3mdet
real(8) r3taxi
external i3mdet,r3taxi
! identity only
if (nsymlat.eq.1) then
  do is1=1,nspecies
    do ia1=1,natoms(is1)
      do ia2=1,natoms(is1)
        nsymeqat(ia1,ia2,is1)=0
      end do
      nsymeqat(ia1,ia1,is1)=1
      symeqat(1,ia1,ia1,is1)=1
      tvleqat(:,1,ia1,ia1,is1)=0.d0
    end do
  end do
  return
end if
allocate(used(ld2))
allocate(av1(3,ld2,nspecies))
allocate(av2(3,ld2,nspecies))
do is1=1,nspecies
  do ia1=1,natoms(is1)
! generate a set of positions with ia1 at the origin
    do is2=1,nspecies
      do ia2=1,natoms(is2)
        av1(:,ia2,is2)=atposl(:,ia2,is2)-atposl(:,ia1,is1)
        call r3frac(eps,av1(1,ia2,is2),id)
      end do
    end do
    do ia2=1,natoms(is1)
! generate a set of positions with ia2 at the origin
      do is2=1,nspecies
        do ia3=1,natoms(is2)
          av2(:,ia3,is2)=atposl(:,ia3,is2)-atposl(:,ia2,is1)
          call r3frac(eps,av2(1,ia3,is2),id)
        end do
      end do
      i=0
! determine if ia1 and ia2 are equivalent for some symmetry
      do isym=1,nsymlat
        s(:,:)=dble(symlat(:,:,isym))
! determinant of symmetry matrix
        md=i3mdet(symlat(1,1,isym))
        do is2=1,nspecies
          used(:)=.false.
          do ia3=1,natoms(is2)
            call r3mv(s,av2(1,ia3,is2),v1)
            call r3frac(eps,v1,id)
            call r3mv(s,plrvl(1,ia3,is2),v2)
            if (.not.tpolar) v2(:)=dble(md)*v2(:)
            do ia4=1,natoms(is2)
              if (.not.used(ia4)) then
                if ((r3taxi(v1,av1(1,ia4,is2)).lt.eps).and. &
                 (r3taxi(v2,plrvl(1,ia4,is2)).lt.eps)) then
                  used(ia4)=.true.
                  goto 10
                end if
              end if
            end do
            goto 20
10 continue
          end do
        end do
        i=i+1
        symeqat(i,ia2,ia1,is1)=isym
! find translation vector for symmetry applied about the origin
        call r3mv(s,atposl(1,ia2,is1),v1)
        tvleqat(:,i,ia2,ia1,is1)=atposl(:,ia1,is1)-v1(:)
20 continue
! end loop over lattice symmetries
      end do
      nsymeqat(ia2,ia1,is1)=i
! end loop over ia2
    end do
! end loop over ia1
  end do
! end loop over species
end do
deallocate(used,av1,av2)
return
end subroutine
!EOC
