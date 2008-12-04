
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mixerifc(mtype,n,v,dv,mode)
use modmain
use modmixermsec
use modmixadapt
implicit none
! arguments
integer, intent(in) :: mtype
integer, intent(in) :: n
real(8), intent(inout) :: v(n)
real(8), intent(out) :: dv
integer, intent(inout) :: mode
!mode: 	-1 call initialisation routines,
!		-2:call destructor
!		else ignore
integer, parameter :: maxsd=3
select case(mtype)
case(1)
! adaptive linear mixing
! calculate memmory requirement if mode negative
  if (mode .eq. -1) then
    mode=0
    write(60,*)"Using Adaptive step size linear potential mixing (1)"
    if(allocated(work))deallocate(work)
    allocate(work(3*n))
    return
  end if
    if (mode .eq. -2) then
    deallocate(work)
    return
  end if
!--
  call mixadapt(iscl,beta0,betainc,betadec,n,v,work,work(n+1),work(2*n+1),dv)

case(2)
 ! multicecant broyden
  if (mode .eq. -1) then
    call initmixermsec(n)
	mode =0
	 write(60,*)"Using Muttisecant Broyden potential mixing (2)"
    return
   end if
    if (mode .eq. -2) then
    call freearraysmixermsec()
     return
  end if
 call  mixmsec(iscl,v,dv,n)
 case(3)
 ! Pulay mixing
   if (spinpol) then
     write(*,*)
     write(*,'("Warning(mixerifc): Pulay mixing problematic with spin-polarised&
      & calculations")')
   end if
   if (mode.eq.-1) then
 write(60,*)"Using Pulay potential mixing (3)"
     allocate(work(2*n*maxsd))
        mode=0
     return
   end if
      if (mode.eq.-2) then
     deallocate(work)

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

