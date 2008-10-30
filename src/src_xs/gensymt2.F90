
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gensymt2(maxsymcrys,nsymcrys,symlatc,lsplsymc,symt2)
  implicit none
  ! arguments
  integer, intent(in) :: maxsymcrys
  integer, intent(in) :: nsymcrys
  real(8) :: symlatc(3,3,48)  
  integer, intent(in) :: lsplsymc(maxsymcrys)
  real(8), intent(out) :: symt2(3,3,3,3)
  ! local variables
  integer :: isym,i,j,iop1,iop2
  real(8) :: s(3,3),sc(3,3)
  do iop1=1,3
     do iop2=1,3
        s(:,:)=0.d0
        do isym=1,nsymcrys
           sc(:,:)=dble(symlatc(:,:,lsplsymc(isym)))
           do i=1,3
              do j=1,3
                 s(i,j)=s(i,j)+sc(i,iop1)*sc(j,iop2)
              end do
           end do
        end do
        symt2(iop1,iop2,:,:)=s(:,:)/dble(nsymcrys)
     end do
  end do
end subroutine gensymt2
