
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module indices
  implicit none
contains

  integer function idx_lex(idx,le,lo,up,rev)
    implicit none
    ! arguments
    integer, intent(in) :: idx(:)
    integer, optional, intent(in) :: lo(:),up(:),le(:)
    logical, optional, intent(in) :: rev
    ! local variables
    integer :: s_idx(1),s_lo(1),s_up(1),s_le(1)
    integer :: ilo(size(idx)),iup(size(idx)),ile(size(idx)),i,j,p,n,s
    logical :: pass,tlou,tle,re
    ! check for either lower and upper bounds or lengths for index range
    pass=.false.
    if (present(lo).and.present(up)) tlou=.true.
    if (present(le)) tle=.true.
    if ((tlou.and.(.not.tle)).or.((.not.tlou).and.tle)) pass=.true.
    if (.not.pass) then
       write(*,*)
       write(*,'("Error(idx_lex): specify either lower and upper bounds or &
            & lengts of index ranges")')
       write(*,*)
       stop
    end if
    pass=.true.
    s_idx=shape(idx)   
    if (tlou) then
       s_lo=shape(lo)
       s_up=shape(up)
       if (any(s_lo.ne.s_up)) pass=.false.
       if (any(s_lo.ne.s_idx)) pass=.false.
       if (.not.pass) then
          write(*,*)
          write(*,'("Error(idx_lex): sizes of indices and bounds not &
               &consistent")')
          write(*,*)
          stop
       end if
    end if
    if (tle) then
       s_le=shape(le)
       if (any(s_le.ne.s_idx)) pass=.false.
       if (.not.pass) then
          write(*,*)
          write(*,'("Error(idx_lex): sizes of indices and lengths not &
               &consistent")')
          write(*,*)
          stop
       end if
    end if
    if (tlou) then
       ilo=lo
       iup=up
       ile=up-lo+1
    end if
    if (tle) then
       ilo=1
       iup=le
       ile=le
    end if
    re=.false.
    if (present(rev)) re=rev
    ! dimension
    n=s_idx(1)
    if (re) then
       ! reverse lexikographical order
       s=0
       do i=1,n
          p=1
          do j=1,i-1
             p=p*ile(n-j+1)
          end do
          s=s+p*(idx(n-i+1)-lo(n-i+1)+1)
       end do       
    else
       ! lexikographical order
       s=0
       do i=1,n
          p=1
          do j=1,i-1
             p=p*ile(j)
          end do
          s=s+p*(idx(i)-lo(i)+1)
       end do       
    end if
    idx_lex=s
  end function idx_lex

  !/////////////////////////////////////////////////////////////////////////////

  subroutine idxi_lex(l,idx,le,lo,up,rev)
    implicit none
    ! arguments
    integer, intent(in) :: l
    integer, intent(out) :: idx(:)
    integer, optional, intent(in) :: lo(:),up(:),le(:)
    logical, optional, intent(in) :: rev
    ! local variables
  end subroutine idxi_lex


end module indices

!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////

integer function idxkkp2(ik,ikp,n)
  use indices
  implicit none
  integer, intent(in) :: ik,ikp,n
  integer :: idx(2),lo(2),up(2),m(2)
  idx=(/ik,ikp/)
  m(:)=n
  idxkkp2=idx_lex(idx,m)
end function idxkkp2

