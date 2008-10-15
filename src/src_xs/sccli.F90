
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_sccli
  implicit none
contains
  subroutine sccli(s,lb,ub,emat13,emat24,sci,scimat)
    use indices
    implicit none
    ! arguments
    real(8), intent(in) :: s
    integer, intent(in) :: lb(4),ub(4)
    complex(8), intent(in) :: emat13(:,:),emat24(:,:),sci(:,:)
    complex(8), intent(out) :: scimat(:,:)
    ! local variables    
    integer :: siz(4),j12,j34,j13,j24,n12,n34,n13,n24
    complex(8), allocatable :: scm(:,:)
    siz(:)=ub(:)-lb(:)+1
    n12=siz(1)*siz(2); n34=siz(3)*siz(4); n13=siz(1)*siz(3); n24=siz(2)*siz(4)
    allocate(scm(n13,n24))
    ! carry out the double summation
    scm(:,:)=s*matmul(transpose(conjg(emat13)),matmul(sci,emat24))
    do j13=1,n13
       do j24=1,n24
          call idxt1324(j12,j34,j13,j24,lb,ub,inv=.true.)
          scimat(j12,j34)=scm(j13,j24)
       end do
    end do
    deallocate(scm)
  end subroutine sccli
end module m_sccli
