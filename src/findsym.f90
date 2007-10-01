
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findsym
! !INTERFACE:
subroutine findsym(apl1,apl2,nsym,lspl,lspn,iea)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   apl1 : first set of atomic positions in lattice coordinates
!          (in,real(3,maxatoms,maxspecies))
!   apl2 : second set of atomic positions in lattice coordinates
!          (in,real(3,maxatoms,maxspecies))
!   nsym : number of symmetries (out,integer)
!   lspl : spatial rotation element in lattice point group for each symmetry
!          (out,integer(48))
!   lspn : spin rotation element in lattice point group for each symmetry
!          (out,integer(48))
!   iea  : equivalent atom index for each symmetry
!          (out,integer(iea(natmmax,nspecies,48))
! !DESCRIPTION:
!   Finds the symmetries which rotate one set of atomic positions into another.
!   Both sets of positions differ only by a translation vector and have the same
!   muffin-tin magnetic fields (stored in the global array {\tt bfcmt}). Any
!   symmetry element consists of a spatial rotation of the atomic position
!   vectors followed by a global magnetic rotation: $\{\alpha_S|\alpha_R\}$. In
!   the case of spin-orbit coupling $\alpha_S=\alpha_R$. The symmetries are
!   returned as indices of elements in the Bravais lattice point group. An
!   index to equivalent atoms is stored in the array {\tt iea}.
!
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: apl1(3,maxatoms,maxspecies)
real(8), intent(in) :: apl2(3,maxatoms,maxspecies)
integer, intent(out) :: nsym
integer, intent(out) :: lspl(48)
integer, intent(out) :: lspn(48)
integer, intent(out) :: iea(natmmax,nspecies,48)
! local variables
integer isym,jsym,jsym0,jsym1
integer is,ia,ja,iv(3),md
real(8) v(3),t1
real(8) sl(3,3,48),sc(3,3,48)
! automatic arrays
integer jea(natmmax,nspecies)
real(8) apl3(3,natmmax)
! external functions
integer i3mdet
real(8) r3taxi
external i3mdet,r3taxi
do isym=1,nsymlat
! make real copy of lattice rotation symmetries
  sl(:,:,isym)=dble(symlat(:,:,isym))
! determine lattice symmetry matrix in Cartesian coordinates
  call r3mm(sl(1,1,isym),ainv,sc(1,1,isym))
  call r3mm(avec,sc(1,1,isym),sc(1,1,isym))
end do
nsym=0
! loop over lattice symmetries (spatial rotations)
do isym=1,nsymlat
  do is=1,nspecies
! map apl1 coordinates to [0,1) and store in apl3
    do ia=1,natoms(is)
      apl3(:,ia)=apl1(:,ia,is)
      call r3frac(epslat,apl3(1,ia),iv)
    end do
    do ja=1,natoms(is)
! apply lattice symmetry to atomic positions
      call r3mv(sl(1,1,isym),apl2(1,ja,is),v)
! map coordinates to [0,1)
      call r3frac(epslat,v,iv)
! check if atomic positions are invariant
      do ia=1,natoms(is)
        t1=r3taxi(apl3(1,ia),v)
        if (t1.lt.epslat) then
! equivalent atom index
          jea(ia,is)=ja
          goto 10
        end if
      end do
! not invariant so try new spatial rotation
      goto 40
10 continue
    end do
  end do
! all atomic positions invariant at this point
  jsym=1
! spin polarised case
  if (spinpol) then
! check invariance of magnetic fields under global spin rotation
    if (spinorb) then
! with spin-orbit coupling spin rotation equals spatial rotation
      jsym0=isym
      jsym1=isym
    else
! without spin-orbit coupling spin rotation independent of spatial rotation
      jsym0=1
      jsym1=nsymlat
    end if
    do jsym=jsym0,jsym1
! only use proper rotations for spin
      md=i3mdet(symlat(1,1,jsym))
      if (md.lt.0) goto 20
! rotate global field and check invariance
      call r3mv(sc(1,1,jsym),bfieldc,v)
      t1=r3taxi(bfieldc,v)
! if not invariant try a different global spin rotation
      if (t1.gt.epslat) goto 20
! rotate muffin-tin magnetic fields and check invariance
      do is=1,nspecies
        do ia=1,natoms(is)
! equivalent atom
          ja=jea(ia,is)
          call r3mv(sc(1,1,jsym),bfcmt(1,ja,is),v)
          t1=r3taxi(bfcmt(1,ia,is),v)
! if not invariant try a different global spin rotation
          if (t1.gt.epslat) goto 20
        end do
      end do
! all fields invariant
      goto 30
20 continue
! end loop over global spin rotations
    end do
! magnetic fields not invariant so try different spatial rotation
    goto 40
  end if
30 continue
! everything invariant so add symmetry to set
  nsym=nsym+1
  lspl(nsym)=isym
  lspn(nsym)=jsym
  do is=1,nspecies
    do ia=1,natoms(is)
      iea(ia,is,nsym)=jea(ia,is)
    end do
  end do
40 continue
! end loop over spatial rotations
end do
return
end subroutine
!EOC

