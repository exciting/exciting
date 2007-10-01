
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: findsymcrys
! !INTERFACE:
subroutine findsymcrys(tsymctr,eps,nspecies,natoms,nsymlat,symlat,ld,atposl, &
 tpolar,plrvl,nsymcrys,symcrys)
! !INPUT/OUTPUT PARAMETERS:
!   tsymctr  : .true. if the crystal should be shifted to the optimal symmetry
!              center (in,logical)
!   eps      : zero vector tolerance (in,real)
!   nspecies : number of species (in,integer)
!   natoms   : number atoms for each species (in,integer(nspecies))
!   nsymlat  : number of lattice point group symmetries (in,integer)
!   symlat   : lattice point group symmetries (in,integer(3,3,48))
!   ld       : leading dimension (in,integer)
!   atposl   : atomic positions in lattice coordinates
!              (inout,real(3,ld,nspecies))
!   tpolar   : .true. if plrvl are polar vectors, .false. if they are axial
!              vectors (in,logical)
!   plrvl    : polarisation vector at each atom site in lattice coordinates
!              (in,real(3,ld,nspecies))
!   nsymcrys : number of crystal point group symmetries found (out,integer)
!   symcrys  : crystal point group symmetries found (out,integer(3,3,48))
! !DESCRIPTION:
!   Finds the largest number of point group symmetries which leave the crystal
!   structure, including atomic positions and polarisation vectors at each site,
!   invariant. If {\tt tpolar} is {\tt .true.} then the polarisation vectors are
!   taken to be polar, otherwise they are axial and invariant under inversion.
!   All atomic positions as well as the mid-points between two and three atoms
!   are checked as possible symmetry centers. The atomic positions are then
!   shifted so that the new symmetry center lies at $(0,0,0)$. The routine
!   requires the symmetries which leave the Bravais lattice invariant. See
!   routine {\tt findsymcrysctr}.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsymctr
real(8), intent(in) :: eps
integer, intent(in) :: nspecies
integer, intent(in) :: natoms(nspecies)
integer, intent(in) :: nsymlat
integer, intent(in) :: symlat(3,3,48)
integer, intent(in) :: ld
real(8), intent(inout) :: atposl(3,ld,nspecies)
logical, intent(in) :: tpolar
real(8), intent(in) :: plrvl(3,ld,nspecies)
integer, intent(out) :: nsymcrys
integer, intent(out) :: symcrys(3,3,48)
! local variables
integer is,ia1,ia2,ia3,i1,i2,i3
integer nsym,maxsym
integer sym(3,3,48)
real(8) v(3),ctrvl1(3),ctrvl2(3)
if (nspecies.lt.0) then
  write(*,*)
  write(*,'("Error(findsymcrys): invalid nspecies : ",I8)') nspecies
  write(*,*)
  stop
end if
if ((nsymlat.le.0).or.(nsymlat.gt.48)) then
  write(*,*)
  write(*,'("Error(findsymcrys): invalid nsymlat : ",I8)') nsymlat
  write(*,*)
  stop
end if
! identity only
if (nsymlat.eq.1) then
  nsymcrys=1
  symcrys(:,:,1)=symlat(:,:,1)
  return
end if
if (.not.tsymctr) goto 10
! find the optimal symmetry center
ctrvl1(:)=0.d0
call findsymcrysctr(eps,nspecies,natoms,nsymlat,symlat,ld,atposl,tpolar,plrvl, &
 ctrvl1,maxsym,sym)
do is=1,nspecies
! primitive translations
  do i1=-1,1
    v(1)=dble(i1)
    do i2=-1,1
      v(2)=dble(i2)
      do i3=-1,1
        v(3)=dble(i3)
        do ia1=1,natoms(is)
! check atom sites
          ctrvl2(:)=atposl(:,ia1,is)+v(:)
          call findsymcrysctr(eps,nspecies,natoms,nsymlat,symlat,ld,atposl, &
           tpolar,plrvl,ctrvl2,nsym,sym)
          if (nsym.gt.maxsym) then
            maxsym=nsym
            ctrvl1(:)=ctrvl2(:)
          end if
! check mid-point between two atoms
          do ia2=1,ia1-1
            ctrvl2(:)=(atposl(:,ia1,is)+atposl(:,ia2,is)+v(:))/2.d0
            call findsymcrysctr(eps,nspecies,natoms,nsymlat,symlat,ld,atposl, &
             tpolar,plrvl,ctrvl2,nsym,sym)
            if (nsym.gt.maxsym) then
              maxsym=nsym
              ctrvl1(:)=ctrvl2(:)
            end if
! check mid-point between three atoms
            do ia3=1,ia2-1
              ctrvl2(:)=(atposl(:,ia1,is)+atposl(:,ia2,is)+atposl(:,ia3,is) &
               +v(:))/3.d0
              call findsymcrysctr(eps,nspecies,natoms,nsymlat,symlat,ld, &
               atposl,tpolar,plrvl,ctrvl2,nsym,sym)
              if (nsym.gt.maxsym) then
                maxsym=nsym
                ctrvl1(:)=ctrvl2(:)
              end if
            end do
          end do
        end do
! end loops over i1, i2 and i3
      end do
    end do
  end do
! end loop over species
end do
! re-center the atomic coordinates
do is=1,nspecies
  do ia1=1,natoms(is)
    atposl(:,ia1,is)=atposl(:,ia1,is)-ctrvl1(:)
  end do
end do
10 continue
ctrvl2(:)=0.d0
call findsymcrysctr(eps,nspecies,natoms,nsymlat,symlat,ld,atposl,tpolar,plrvl, &
 ctrvl2,nsymcrys,symcrys)
return
end subroutine
!EOC

