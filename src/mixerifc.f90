
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mixerifc(mtype,n,v,dv,nwork,work)
use modmain
implicit none
! arguments
integer, intent(in) :: mtype
integer, intent(in) :: n
real(8), intent(inout) :: v(n)
real(8), intent(out) :: dv
integer, intent(inout) :: nwork
real(8), intent(inout) :: work(*)
! local variables
! maximum subspace dimension for the Pulay mixer
integer, parameter :: maxsd=3
select case(mtype)
case(1)
! adaptive linear mixing
  if (nwork.le.0) then
    nwork=3*n
    return
  end if
  call mixadapt(iscl,beta0,betainc,betadec,n,v,work,work(n+1),work(2*n+1),dv)
case(2)
! Pulay mixing
  if (spinpol) then
    write(*,*)
    write(*,'("Warning(mixerifc): Pulay mixing problematic with spin-polarised&
     & calculations")')
  end if
  if (nwork.le.0) then
    nwork=2*n*maxsd
    return
  end if
  call mixpulay(iscl,n,maxsd,v,work,work(n*maxsd+1),dv)
case default
  write(*,*)
  write(*,'("Error(mixerifc): mtype not defined : ",I8)') mtype
  write(*,*)
  stop
end select
return
end subroutine

