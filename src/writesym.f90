
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
integer isym,i,is,ia1,ia2,ias
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
write(50,'(I4," : nsymcrys")') nsymcrys
do isym=1,nsymcrys
  write(50,*)
  write(50,'(I4)') isym
  do i=1,3
    write(50,'(3I4)') symcrys(i,:,isym)
  end do
end do
close(50)
! output the site symmetries
open(50,file='SYMSITE'//trim(filext),action='WRITE',form='FORMATTED')
do is=1,nspecies
  do ia1=1,natoms(is)
    ias=idxas(ia1,is)
    write(50,*)
    write(50,'(2I4," : atom, species")') ia1,is
    write(50,'(I4," : nsymsite")') nsymsite(ias)
    do isym=1,nsymsite(ias)
      write(50,*)
      write(50,'(I4)') isym
      do i=1,3
        write(50,'(3I4)') symsite(i,:,isym,ias)
      end do
    end do
  end do
end do
close(50)
! output the equivalent atoms and related symmetries
open(50,file='EQATOMS'//trim(filext),action='WRITE',form='FORMATTED')
write(50,'("Operations are represented by a lattice point group symmetry &
 &applied about the")')
write(50,'("origin, followed by a translation vector in lattice &
 &coordinates")')
do is=1,nspecies
  write(50,*)
  write(50,'("Species : ",I4,", ",A)') is,trim(spsymb(is))
  write(50,*)
  do ia1=1,natoms(is)
    do ia2=1,natoms(is)
      write(50,'(" atom ",I4," can be mapped to atom ",I4, " with ",I4,&
       &" operations")') ia2,ia1,nsymeqat(ia2,ia1,is)
      do i=1,nsymeqat(ia2,ia1,is)
        write(50,'(I4,3F12.8)') symeqat(i,ia2,ia1,is),tvleqat(:,i,ia2,ia1,is)
      end do
      write(50,*)
    end do
  end do
end do
close(50)
return
end subroutine
!EOC
