

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module indices
  implicit none
contains

  ! generic routines -----------------------------------------------------------

!BOP
! !ROUTINE: lidx
! !INTERFACE:  
  integer function lidx(idx, lb, ub, rev)
! !DESCRIPTION:
!   Linear index map. The first (last) index of index array a runs fastest.
!   The map is given by
!   $$ L(a_1,\ldots,a_n) = 1+\sum_{i=1}^n (a_i-l_i)\prod_{j=1}^{i-1}m_j, $$
!   where $n$ is the dimension of the index-array $l_i$ and $u_i$ are the lower
!   and upper bounds $l_i\le a_i \le u_i$, respectively and $m_i=u_i-l_i+1$.
!
! !REVISION HISTORY:
!   Created October 2008 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in) :: idx(:), lb(:), ub(:)
    logical, optional, intent(in) :: rev
    ! local variables
    logical :: trev
    integer :: n, i
    integer, allocatable :: siz(:), prd(:, :)
    trev=.false.
    if (present(rev)) trev=rev
    n=size(idx)
    allocate(siz(n), prd(2, n))
    siz(:)=ub(:)-lb(:)+1
    prd(:, 1)=1
    do i=2, n
       prd(1, i)=prd(1, i-1)*siz(i-1)
       prd(2, i)=prd(2, i-1)*siz(n-i+2)
    end do
    lidx=1
    if (trev) then
       do i=1, n
	  lidx=lidx+(idx(n-i+1)-lb(n-i+1))*prd(2, i)
       end do
    else
       do i=1, n
	  lidx=lidx+(idx(i)-lb(i))*prd(1, i)
       end do
    end if
    deallocate(siz, prd)
  end function lidx
!EOC

!BOP
! !ROUTINE: ilidx
! !INTERFACE:  


subroutine ilidx(idx, lidx, lb, ub, rev)
! !DESCRIPTION:
!   Inverse linear index map. The first (last) index of index array a runs
!   fastest. The map is given by
!   \begin{align}
!    a_i &=\frac{\mod(L-1,\prod_{j=1}^i m_j)}{\prod_{j=1}^{i-1} m_j}, 
!      i=1,\ldots,n-1 \&
!    a_n &=\frac{L-1}{\prod_{j=1}^{n-1} m_j}, 
!   \end{align}
!   where $n$ is the dimension of the index-array $l_i$ and $u_i$ are the lower
!   and upper bounds $l_i\le a_i \le u_i$, respectively and $m_i=u_i-l_i+1$.
!
! !REVISION HISTORY:
!   Created October 2008 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(out) :: idx(:)
    integer, intent(in) :: lidx, lb(:), ub(:)
    logical, optional, intent(in) :: rev
    ! local variables
    logical :: trev
    integer :: n, i, lt
    integer, allocatable :: siz(:), prd(:, :)
    trev=.false.
    if (present(rev)) trev=rev
    n=size(idx)
    allocate(siz(n), prd(2, n))
    siz(:)=ub(:)-lb(:)+1
    lt=lidx-1
    prd(:, 1)=1
    do i=2, n
       prd(1, i)=prd(1, i-1)*siz(i-1)
       prd(2, i)=prd(2, i-1)*siz(n-i+2)
    end do
    if (trev) then
       do i=1, n-1
	  idx(n-i+1)=mod(lt, prd(2, i+1))/prd(2, i)+lb(n-i+1)
       end do
       idx(1)=lt/prd(2, n)+lb(1)
    else
       do i=1, n-1
	  idx(i)=mod(lt, prd(1, i+1))/prd(1, i)+lb(i)
       end do
       idx(n)=lt/prd(1, n)+lb(n)
    end if
    deallocate(siz, prd)
  end subroutine ilidx
!EOC

  ! wrapper routines -----------------------------------------------------------

  integer function lidx2(i1, i2, l1, l2, u1, u2, rev)
    implicit none
    ! arguments
    integer, intent(in) :: i1, i2, l1, l2, u1, u2
    logical, optional :: rev
    if (present(rev)) then
       lidx2=lidx((/i1, i2/), (/l1, l2/), (/u1, u2/), rev)
    else
       lidx2=lidx((/i1, i2/), (/l1, l2/), (/u1, u2/))
    end if
  end function lidx2


subroutine ilidx2(i1, i2, il, l1, l2, u1, u2, rev)
    implicit none
    ! arguments
    integer, intent(out) :: i1, i2
    integer, intent(in) :: il, l1, l2, u1, u2
    logical, optional, intent(in) :: rev
    ! local variables
    integer :: idx(2)
    call ilidx(idx, il, (/l1, l2/), (/u1, u2/), rev)
    i1=idx(1)
    i2=idx(2)
  end subroutine ilidx2


subroutine idxt1324(j13, j24, i12, i34, lb, ub, inv)
    implicit none
    ! arguments
    integer, intent(out) :: j13, j24
    integer, intent(in) :: i12, i34, lb(4), ub(4)
    logical, optional, intent(in) :: inv
    ! local variables
    integer :: s1, s2, s3, s4, lbt(4), ubt(4)
    logical :: tinv
    tinv=.false.
    ! inverse transformation
    if (present(inv)) tinv=inv
    lbt(:)=lb(:); ubt(:)=ub(:)
    if (tinv) then
       lbt(:)=lb((/1, 3, 2, 4/))
       ubt(:)=ub((/1, 3, 2, 4/))
    end if
    call ilidx2(s1, s2, i12, lbt(1), lbt(2), ubt(1), ubt(2))
    call ilidx2(s3, s4, i34, lbt(3), lbt(4), ubt(3), ubt(4))
    j13=lidx2(s1, s3, lbt(1), lbt(3), ubt(1), ubt(3))
    j24=lidx2(s2, s4, lbt(2), lbt(4), ubt(2), ubt(4))
  end subroutine idxt1324

  ! tests ----------------------------------------------------------------------


subroutine test_indices
    implicit none
    integer :: c, j1, j2, jj(2), i1, i2, l1, l2, u1, u2, li, ii(2)
    integer :: ll(2), uu(2)
    logical :: pass
    l1=-1; u1=7; l2=-2; u2=5
    ll=(/l1, l2/); uu=(/u1, u2/)
    c=0
    do i2=l2, u2
       do i1=l1, u1
	  c=c+1
	  ii=(/i1, i2/)
	  li=lidx(ii, ll, uu)
	  call ilidx(jj, c, ll, uu)
	  j1=jj(1); j2=jj(2)
	  pass=(i1.eq.j1).and.(i2.eq.j2).and.(c.eq.li)
	  write( * , '(a, 3i4, 3x, 3i4, 3x, l)') 'loop:i1, i2, c|ifun:j1, j2, li:', &
	       i1, i2, c, j1, j2, li, pass
	  call ilidx2(j1, j2, c, l1, l2, u1, u2)
	  li=lidx2(i1, i2, l1, l2, u1, u2)
	  pass=(i1.eq.j1).and.(i2.eq.j2).and.(c.eq.li)
	  write(*, '(3i5, 3x, l)') j1, j2, li, pass
       end do
    end do
  end subroutine test_indices


subroutine test_indices2
    implicit none
    integer :: i1, i2, i3, i4, l1, l2, l3, l4, u1, u2, u3, u4, c12, c34, j12, j34, j13, j24
    integer :: lb(4), ub(4)
    logical :: pass
    l1=-1; u1=7; l2=-2; u2=5
!    l3=-1; u3=7; l4=-2; u4=5    ! same limits as for 1,2 indices
    l3=-2; u3=4; l4=3; u4=5    ! different limits as for 1, 2 indices
    lb=(/l1, l2, l3, l4/)
    ub=(/u1, u2, u3, u4/)
    c34=0
    do i4=l4, u4
       do i3=l3, u3
	  c34=c34+1
	  c12=0
	  do i2=l2, u2
	     do i1=l1, u1
		c12=c12+1
		call idxt1324(j13, j24, c12, c34, lb, ub)
		call idxt1324(j12, j34, j13, j24, lb, ub, inv=.true.)
!                call idxt1324(j12,j34,j13,j24,lb((/1,3,2,4/)),ub((/1,3,2,4/)))
		pass=(c12.eq.j12).and.(c34.eq.j34)
		write(*, '(2i4, 3x, 2i4, 3x, 2i4, 3x, l)') c12, c34, j12, j34, j13, j24, pass
	     end do
	  end do
       end do
    end do
  end subroutine test_indices2

end module indices
