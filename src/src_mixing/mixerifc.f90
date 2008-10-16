
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mixerifc(mtype,n,v,dv,nwork)
use modmain
use modmixermsec
implicit none
! arguments
integer, intent(in) :: mtype
integer, intent(in) :: n
real(8), intent(inout) :: v(n)
real(8), intent(out) :: dv
integer, intent(inout) :: nwork

select case(mtype)
case(1)
! adaptive linear mixing
! calculate memmory requirement if nwork negative
  if (nwork .eq. -1) then
    nwork=3*n
    if(allocated(work))deallocate(work)
    allocate(work(nwork))
    return
  end if
    if (nwork .eq. -2) then
    deallocate(work)
    return
  end if
!--
  call mixadapt(iscl,beta0,betainc,betadec,n,v,work,work(n+1),work(2*n+1),dv)

case(2)
 ! multicecant broyden
  if (nwork .eq. -1) then
    call initmixermsec(n)
    nwork=0
    return
   end if
    if (nwork .eq. -2) then
    call freearraysmixermsec()
     return
  end if
 call  mixmsec(iscl,v,dv,n)
case default
  write(*,*)
  write(*,'("Error(mixerifc): mtype not defined : ",I8)') mtype
  write(*,*)
  stop
end select
return
end subroutine

