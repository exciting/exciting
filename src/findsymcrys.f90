
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findsymcrys
! !INTERFACE:
subroutine findsymcrys
! !USES:
use modmain
! !DESCRIPTION:
!   Finds the complete set of symmetries which leave the crystal structure
!   (including the magnetic fields) invariant. A crystal symmetry is of the
!   form $\{\alpha_S|\alpha_R|{\bf t}\}$, where ${\bf t}$ is a translation
!   vector, $\alpha_R$ is a spatial rotation operation and $\alpha_S$ is a
!   global spin rotation. Note that the order of operations is important and
!   defined to be from right to left, i.e. translation followed by spatial
!   rotation followed by spin rotation. In the case of spin-orbit coupling
!   $\alpha_S=\alpha_R$. The translation vectors are determined by checking
!   all displacements of the form $(\frac{i}{12},\frac{j}{12},\frac{k}{12})$
!   for $i,j,k=0\ldots 11$. See L. M. Sandratskii and P. G. Guletskii, {\it J.
!   Phys. F: Met. Phys.} {\bf 16}, L43 (1986) and the routine {\tt findsym}.
!
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!EOP
!BOC
implicit none
! local variables
integer, parameter :: ndiv=12
integer ia,ja,is,js
integer i1,i2,i3
integer isym,nsym
integer lspl(48),lspn(48)
real(8) v(3),vtl(3),t1,t2
real(8) apl(3,maxatoms,maxspecies)
! allocatable arrays
integer, allocatable :: iea(:,:,:)
! allocate local arrays
allocate(iea(natmmax,nspecies,48))
! allocate equivalent atom arrays
if (allocated(ieqatom)) deallocate(ieqatom)
allocate(ieqatom(natmmax,nspecies,maxsymcrys))
if (allocated(eqatoms)) deallocate(eqatoms)
allocate(eqatoms(natmmax,natmmax,nspecies))
if (tshift) then
! shift basis so that the atom closest to the origin is exactly at the origin
  js=1
  ja=1
  t1=1.d8
  do is=1,nspecies
    do ia=1,natoms(is)
      t2=sqrt(atposc(1,ia,is)**2+atposc(2,ia,is)**2+atposc(3,ia,is)**2)
      if (t2.lt.t1+epslat) then
        js=is
        ja=ia
        t1=t2
      end if
    end do
  end do
  v(:)=atposl(:,ja,js)
  do is=1,nspecies
    do ia=1,natoms(is)
      atposl(:,ia,is)=atposl(:,ia,is)-v(:)
      call r3mv(avec,atposl(1,ia,is),atposc(1,ia,is))
    end do
  end do
end if
eqatoms(:,:,:)=.false.
nsymcrys=0
! loop over all possible translations
do i1=0,ndiv-1
  vtl(1)=dble(i1)/dble(ndiv)
  do i2=0,ndiv-1
    vtl(2)=dble(i2)/dble(ndiv)
    do i3=0,ndiv-1
      vtl(3)=dble(i3)/dble(ndiv)
! construct new array with translated positions
      do is=1,nspecies
        do ia=1,natoms(is)
          apl(:,ia,is)=atposl(:,ia,is)+vtl(:)
        end do
      end do
! find the symmetries for current translation
      call findsym(atposl,apl,nsym,lspl,lspn,iea)
      do isym=1,nsym
        nsymcrys=nsymcrys+1
        if (nsymcrys.gt.maxsymcrys) then
          write(*,*)
          write(*,'("Error(findsymcrys): too many symmetries")')
          write(*,'(" Adjust maxsymcrys in modmain and recompile code")')
          write(*,*)
          stop
        end if
        vtlsymc(:,nsymcrys)=vtl(:)
        lsplsymc(nsymcrys)=lspl(isym)
        lspnsymc(nsymcrys)=lspn(isym)
        do is=1,nspecies
          do ia=1,natoms(is)
            ja=iea(ia,is,isym)
            ieqatom(ia,is,nsymcrys)=ja
            eqatoms(ia,ja,is)=.true.
            eqatoms(ja,ia,is)=.true.
          end do
        end do
      end do
    end do
  end do
end do
deallocate(iea)
return
end subroutine
!EOC

