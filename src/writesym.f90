
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writesym
! !INTERFACE:
subroutine writesym
! !USES:
use modmain
! !DESCRIPTION:
!   Outputs the Bravais, crystal and site symmetry matrices to files
!   {\tt SYMLAT.OUT}, {\tt SYMCRYS.OUT} and {\tt SYMSITE.OUT}, respectively.
!   Also writes out equivalent atoms and related crystal symmetries to
!   {\tt EQATOMS.OUT}.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ja,ias,i
integer isym,lspl,lspn
! output the Bravais lattice symmetries
open(50,file='SYMLAT'//trim(filext),action='WRITE',form='FORMATTED')
write(50,'(I4," : nsymlat")') nsymlat
do isym=1,nsymlat
  write(50,*)
  write(50,'(I4)') isym
  do i=1,3
    write(50,'(3I4)') symlat(i,:,isym)
  end do
end do
close(50)
! output the crystal symmetries
open(50,file='SYMCRYS'//trim(filext),action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("(translation vectors and rotation matrices are in lattice &
 &coordinates)")')
write(50,*)
write(50,'(I4," : nsymcrys")') nsymcrys
do isym=1,nsymcrys
  write(50,*)
  write(50,'("Crystal symmetry : ",I4)') isym
  write(50,'(" spatial translation :")')
  write(50,'(3G18.10)') vtlsymc(:,isym)
  write(50,'(" spatial rotation :")')
  lspl=lsplsymc(isym)
  do i=1,3
    write(50,'(3I4)') symlat(i,:,lspl)
  end do
  write(50,'(" global spin rotation :")')
  lspn=lspnsymc(isym)
  do i=1,3
    write(50,'(3I4)') symlat(i,:,lspn)
  end do
end do
close(50)
! output the site symmetries
open(50,file='SYMSITE'//trim(filext),action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("(rotation matrices are in lattice coordinates)")')
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    write(50,'(I4," : nsymsite")') nsymsite(ias)
    do isym=1,nsymsite(ias)
      write(50,*)
      write(50,'(" Site symmetry : ",I4)') isym
      write(50,'("  spatial rotation :")')
      lspl=lsplsyms(isym,ias)
      do i=1,3
        write(50,'(3I4)') symlat(i,:,lspl)
      end do
      write(50,'("  global spin rotation :")')
      lspn=lspnsyms(isym,ias)
      do i=1,3
        write(50,'(3I4)') symlat(i,:,lspn)
      end do
    end do
  end do
end do
close(50)
! output the equivalent atoms and related symmetries
open(50,file='EQATOMS'//trim(filext),action='WRITE',form='FORMATTED')
do is=1,nspecies
  write(50,*)
  write(50,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
  do ia=1,natoms(is)
    write(50,'(" atom ",I4," is equivalent to atom(s)")') ia
    i=0
    do ja=1,natoms(is)
      if (eqatoms(ia,ja,is)) then
        if ((i.gt.0).and.(mod(i,20).eq.0)) write(50,*)
        write(50,'(I4)',advance='NO') ja
        i=i+1
      end if
    end do
    write(50,*)
  end do
end do
close(50)
return
end subroutine
!EOC

