
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module indices
  implicit none
contains

  integer function lidx(idx,lb,ub,rev)
    ! Linear index map. First index of index array a runs fastest
    ! If trev is .true., last index runs fastest
    ! lidx = idx1-lb1+1 + n1[idx2-lb2+1] + ...
    ! where nj=ubj-lbj+1
    implicit none
    ! arguments
    integer, intent(in) :: idx(:),lb(:),ub(:)
    logical, optional :: rev
    ! local variables
    logical :: trev
    integer :: n,i,j
    integer, allocatable :: siz(:),prd(:,:)
    trev=.false.
    if (present(rev)) trev=.true.
    n=size(idx)
    allocate(siz(n),prd(2,n))
    siz(:)=ub(:)-lb(:)+1
    prd(1,1)=1
    prd(2,1)=1
    do i=2,n
       ! n(1)*...*n(i-1)
       prd(1,i)=prd(1,i-1)*siz(i-1)
       ! n(n)*...*n(n-i+2)
       prd(2,i)=prd(2,i-1)*siz(n-i+2)
    end do
    lidx=1
    if (trev) then
       do i=1,n
          lidx=lidx+(idx(n-i+1)-lb(n-i+1))*prd(2,i)
       end do
    else
       do i=1,n
          lidx=lidx+(idx(i)-lb(i))*prd(1,i)
       end do
    end if
    deallocate(siz,prd)
  end function lidx

  subroutine ilidx(idx,lidx,lb,ub,rev)
    ! Linear index map. First index of index array a runs fastest
    ! If trev is .true., last index runs fastest
    ! idx_j=mod(lidx-1,n_1*...*n_j)/(n_1*...*n_(j-1)), j=1,n-1
    ! idx_n=(lidx-1)/(n_1*...*n_(n-1))
    ! where nj=ubj-lbj+1
    implicit none
    ! arguments
    integer, intent(out) :: idx(:)
    integer, intent(in) :: lidx,lb(:),ub(:)
    logical, optional :: rev
    ! local variables
    logical :: trev
    integer :: n,i,j,lt
    integer, allocatable :: siz(:),prd(:,:)
    trev=.false.
    if (present(rev)) trev=.true.
    n=size(idx)
    allocate(siz(n),prd(2,n))
    siz(:)=ub(:)-lb(:)+1
    lt=lidx-1
    prd(1,1)=1
    prd(2,1)=1
    do i=2,n
       ! n(1)*...*n(i-1)
       prd(1,i)=prd(1,i-1)*siz(i-1)
       ! n(n)*...*n(n-i+2)
       prd(2,i)=prd(2,i-1)*siz(n-i+2)
    end do
    if (trev) then
       do i=1,n-1
          idx(n-i+1)=mod(lt,prd(2,i+1))/prd(2,i)+lb(n-i+1)
       end do
       idx(1)=lt/prd(2,n)+lb(1)
    else
       do i=1,n-1
          idx(i)=mod(lt,prd(1,i+1))/prd(1,i)+lb(i)
       end do
       idx(n)=lt/prd(1,n)+lb(n)
    end if
    deallocate(siz,prd)
  end subroutine ilidx

!!$  subroutine test_indices
!!$    implicit none
!!$    integer :: c,j1,j2,jj(2),i1,i2,l1,l2,u1,u2,li,ii(2)
!!$    integer :: ll(2),uu(2),aa(2),bb(2)
!!$    logical :: pass
!!$    l1=-1; u1=7
!!$    l2=-2; u2=5
!!$    ll=(/l1,l2/)
!!$    uu=(/u1,u2/)
!!$    c=0
!!$    do i2=l2,u2
!!$       do i1=l1,u1
!!$          c=c+1
!!$          ii=(/i1,i2/)
!!$          li=lidx(ii,ll,uu)
!!$          call ilidx(jj,c,ll,uu)
!!$          j1=jj(1); j2=jj(2)
!!$          pass=(i1.eq.j1).and.(i2.eq.j2).and.(c.eq.li)
!!$          write(*,'(a,3i4,3x,3i4,3x,l)') 'loop:i1,i2,c|ifun:j1,j2,li:', &
!!$               i1,i2,c,j1,j2,li,pass
!!$       end do
!!$    end do       
!!$  end subroutine test_indices

end module indices
