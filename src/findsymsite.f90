
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findsymsite
use modmain
implicit none
! local variables
integer is,js,ia,ja,ias
real(8) apl(3,maxatoms,maxspecies)
! automatic arrays
real(8) iea(natmmax,nspecies,48)
! allocate the site symmetry arrays
if (allocated(nsymsite)) deallocate(nsymsite)
allocate(nsymsite(natmtot))
if (allocated(lsplsyms)) deallocate(lsplsyms)
allocate(lsplsyms(48,natmtot))
if (allocated(lspnsyms)) deallocate(lspnsyms)
allocate(lspnsyms(48,natmtot))
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do js=1,nspecies
      do ja=1,natoms(js)
        apl(:,ja,js)=atposl(:,ja,js)-atposl(:,ia,is)
      end do
    end do
    call findsym(apl,apl,nsymsite(ias),lsplsyms(:,ias),lspnsyms(:,ias),iea)
  end do
end do
return
end subroutine

