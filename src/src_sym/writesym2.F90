
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writesym2
  use modmain
  use modsym
  implicit none
  ! local varaibles
  integer :: isym,jsym,igenr
  character(32) :: str
  ! write out multiplication table
  open(50,file='SYMMULT'//trim(filext),action='WRITE',form='FORMATTED')
  write(50,*)
  if (abelsg) write(50,'("The symmetry group is Abelian (commutative)")')
  write(50,'(" (first and second group element and product below)")')
  do isym=1,nsymcrys
     do jsym=1,nsymcrys
        write(50,'(2i6,5x,i6)') isym,jsym,sgmut(isym,jsym)
     end do
  end do
  close(50)
  open(50,file='SYMMULT_TABLE'//trim(filext),action='WRITE',form='FORMATTED')
  write(str,*) maxsymcrys
  write(50,*)
  write(50,'(" (symmety group multiplication table below)")')
  do isym=1,nsymcrys
     write(50,'('//trim(str)//'i3.2)') sgmut(isym,:)
  end do
  close(50)
  ! write out generators
  open(50,file='SYMGENR'//trim(filext),action='WRITE',form='FORMATTED')
  write(50,*)
  write(50,'("Number of elements in symmetry group    : ",i6)') nsymcrys
  write(50,'("Number of generators for symmetry group : ",i6)') ngenr
  write(50,*)
  do igenr=1,ngenr
     write(50,'("generating element:",i4," ,number of elemnts in orbit:",i4)')&
       genr(igenr),negenr(igenr)
     write(50,'(" (orbit of generator below)")')
     write(50,'('//trim(str)//'i4)') orbgenr(igenr,:negenr(igenr))
     write(50,*)
  end do
  close(50)
end subroutine writesym2
