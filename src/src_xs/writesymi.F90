
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writesymi
! !INTERFACE:
subroutine writesymi
! !USES:
  use modmain
  use modxs
! !DESCRIPTION:
!   Outputs the crystal and symmetries including their inverse
!   elements to file {\tt SYMINV.OUT}
!
! !REVISION HISTORY:
!   Created December 2007 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  integer :: i,isym,isymi,lspl,lspli,lspn,lspni
  real(8) :: vtlc(3),vtlic(3)
  ! output the crystal symmetries and inverses
  open(50,file='SYMINV'//trim(filext),action='WRITE',form='FORMATTED')
  write(50,'(I4," : nsymcrys")') nsymcrys
  do isym=1,nsymcrys
     ! inverse crystal symmetry
     isymi=scimap(isym)
     ! spatial rotation
     lspl=lsplsymc(isym)
     lspn=lspnsymc(isym)
     ! global spin rotation
     lspli=lsplsymc(isymi)
     lspni=lspnsymc(isymi)
     ! spatial translation
     call r3mv(avec,vtlsymc(1,isym),vtlc)
     call r3mv(avec,vtlsymc(1,isymi),vtlic)
     write(50,*)
     write(50,'("Crystal symmetry, Bravais symmetry : ",2I4)') isym,lspl
     write(50,'(" inverse operations                : ",2I4)') isymi,lspli
     write(50,'(" spatial rotation and inverse (lattice coordinates):")')
     do i=1,3
        write(50,'(3I4,5x,3I4)') symlat(i,:,lspl),symlat(i,:,lspli)
     end do
     write(50,'(" spatial rotation and inverse (Cartesian coordinates):")')
     do i=1,3
        write(50,'(3G18.10,5x,3G18.10)') symlatc(i,:,lspl),symlatc(i,:,lspli)
     end do
     write(50,'(" spatial translation and inverse (lattice coordinates) :")')
     write(50,'(3G18.10,5x,3G18.10)') vtlsymc(:,isym),vtlsymc(:,isymi)
     write(50,'(" spatial translation and inverse (Cartesian coordinates) :")')
     write(50,'(3G18.10,5x,3G18.10)') vtlc,vtlic
     write(50,'(" global spin rotation and inverse (lattice coordinates) :")')
     do i=1,3
        write(50,'(3I4,5x,3I4)') symlat(i,:,lspn),symlat(i,:,lspni)
     end do
     write(50,'(" global spin rotation and inverse (Cartesian coordinates) :")')
     do i=1,3
        write(50,'(3G18.10,5x,3G18.10)') symlatc(i,:,lspn),symlatc(i,:,lspni)
     end do
  end do
  close(50)
end subroutine writesymi
!EOC

