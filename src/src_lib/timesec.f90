
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine timesec(ts)
implicit none
! arguments
real(8), intent(out) :: ts
! local variables
integer count,count_rate
call system_clock(count=count,count_rate=count_rate)
ts=dble(count)/dble(count_rate)
return
end subroutine

