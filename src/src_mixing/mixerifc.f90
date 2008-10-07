
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mixerifc(mtype,n,v,dv,nwork,work)
use modmain
use modmixadapt
implicit none
! arguments
integer, intent(in) :: mtype
integer, intent(in) :: n
real(8), intent(inout) :: v(n)
real(8), intent(out) :: dv
integer, intent(inout) :: nwork
real(8), intent(inout) :: work(*)
select case(mtype)
case(1)
! adaptive linear mixing
  if (nwork.le.0) then
	call    init_mixadapt_arrays(n)
    return
  end if
  call mixadapt(iscl,beta0,betainc,betadec,n,v,dv)
case(2)
  call mixmsec(iscl,v,dv,n)
case default
  write(*,*)
  write(*,'("Error(mixerifc): mtype not defined : ",I8)') mtype
  write(*,*)
  stop
end select
return
end subroutine

