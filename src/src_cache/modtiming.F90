module timing
implicit none

type timer
private
logical :: running
real(8) :: mark,elapsed
end type timer

private :: get_system_time
private :: get_cpu_time

contains

subroutine new(t)
implicit none
type(timer) :: t
t%running=.false.
t%elapsed=0.d0
end subroutine

subroutine delete(t)
implicit none
type(timer) :: t
call new(t)
end subroutine

subroutine tic(t)
implicit none
type(timer) :: t
real(8) :: s
s=get_system_time()
if (t%running) then
  t%elapsed=t%elapsed+s-t%mark
  t%mark=s
else
  t%mark=s
end if
end subroutine

subroutine toc(t)
implicit none
type(timer) :: t
real(8) :: s
if (.not.t%running) then
  write(*,'("Error(modtiming:toc): timer already stopped")')
  stop
end if
s=get_system_time()
t%mark=s
t%elapsed=t%elapsed+s-t%mark
t%running=.false.
end subroutine

!%%%%%%%%%%%%%%%%
!private routines

real(8) function get_system_time()
implicit none
integer :: cnt,cntr
call system_clock(count=cnt,count_rate=cntr)
get_system_time=dble(cnt)/dble(cntr)
end function

real(8) function get_cpu_time()
implicit none
call cpu_time(get_cpu_time)
end function

end module timing



