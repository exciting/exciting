
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modtimer
  implicit none
  type timer
     private
     ! name of timer
     character(256) :: name
     ! type of timer (system clock or CPU clock)
     logical :: systemclock
     ! true if the timer is running
     logical :: running
     ! mark the current time
     real(8) :: mark
     ! the elapsed time
     real(8) :: elapsed
     ! number of measured intervals
     integer :: intervals
     ! the smallest measured interval
     real(8) :: int_small
     ! the largest measured interval
     real(8) :: int_large
     ! the mean value of the measured intervals
     real(8) :: int_mean
  end type timer
  private :: get_system_time
  private :: get_cpu_time
  private :: get_time
  private :: new_interval
contains

  subroutine new_timer(t,typ)
    implicit none
    type(timer) :: t
    character(*), intent(in) :: typ
    t%name=trim(adjustl(typ))
    t%running=.false.
    t%elapsed=0.d0
    t%intervals=0
    t%int_small=1.d10
    t%int_large=-1.d10
    t%int_mean=0.d0
    select case (trim(adjustl(typ)))
       case('system','wall')
          t%systemclock=.true.
       case('cpu')
          t%systemclock=.false.
       case default
          write(*,*)
          write(*,'("Error(timing:new): unknown clocktype: ",a)') trim(typ)
          write(*,*)
          stop
    end select
  end subroutine new_timer

  subroutine tic(t)
    implicit none
    type(timer) :: t
    if (t%running) then
       call new_interval(t,get_time(t))
    else
       t%running=.true.
       t%mark=get_time(t)
    end if
  end subroutine tic

  subroutine toc(t)
    implicit none
    type(timer) :: t
    if (.not.t%running) then
       write(*,'("Error(modtiming:toc): timer already stopped")')
       stop
    end if
    call new_interval(t,get_time(t))
    t%running=.false.
  end subroutine toc

  subroutine report(t)
    implicit none
    type(timer) :: t
    write(*,'("Timer: ",a)') trim(t%name)
    if (t%running) then
       write(*,'(" timer is running")')
    else
       write(*,'(" timer is stopped")')
    end if
    if (t%systemclock) then
       write(*,'(" timer is using system-clock")')
    else
       write(*,'(" timer is using CPU-clock")')
    end if
    write(*,'(" elapsed time      : ",g18.10)') t%elapsed
    write(*,'(" intervals         : ",i8)') t%intervals
    write(*,'(" mean interval     : ",g18.10)') t%int_mean
    write(*,'(" smallest interval : ",g18.10)') t%int_small
    write(*,'(" largest  interval : ",g18.10)') t%int_large
    write(*,'(" current mark      : ",g18.10)') t%mark
  end subroutine report

  !private routines
  subroutine new_interval(t,s)
    implicit none
    type(timer) :: t
    real(8), intent(in) :: s
    real(8) :: sint
    sint=s-t%mark
    t%elapsed=t%elapsed+sint
    t%mark=s
    t%intervals=t%intervals+1
    t%int_small=min(t%int_small,sint)
    t%int_large=max(t%int_large,sint)
    t%int_mean=(dble(t%intervals-1)*t%int_mean+sint)/dble(t%intervals)
  end subroutine new_interval

  real(8) function get_time(t)
    implicit none
    type(timer) :: t
    if (t%systemclock) then
       get_time=get_system_time()
    else
       get_time=get_cpu_time()
    end if
  end function get_time

  real(8) function get_system_time()
    implicit none
    integer :: cnt,cntr
    call system_clock(count=cnt,count_rate=cntr)
    get_system_time=dble(cnt)/dble(cntr)
  end function get_system_time

  real(8) function get_cpu_time()
    implicit none
    call cpu_time(get_cpu_time)
  end function get_cpu_time
end module modtimer


!///////////////////////////////////////////////////////////////////////////


module modtimer2
  use modtimer
  type timer2
     private
     type(timer) :: sys
     type(timer) :: cpu
  end type timer2

contains

  subroutine new_timer2(t)
    implicit none
    type(timer2) :: t
    call new_timer(t%sys)
    call new_timer(t%cpu)
  end subroutine new_timer2

end module modtimer2


!///////////////////////////////////////////////////////////////////////////


program test_timing
  use modtimer
  use modtimer2
  implicit none
  complex(8) :: a(100,100)

  type(timer) :: mysystimer
  type(timer) :: mycputimer

stop

  call new(mysystimer,'system')
  call new(mycputimer,'cpu')
  call report(mysystimer)
  call report(mycputimer)

  call tic(mysystimer)
  call tic(mycputimer)
  call report(mysystimer)
  call report(mycputimer)
  call sleep(1)

a=matmul(a,a)

  call tic(mysystimer)
  call tic(mycputimer)
  call report(mysystimer)
  call report(mycputimer)
  call sleep(2)

  call tic(mysystimer)
  call tic(mycputimer)
  call report(mysystimer)
  call report(mycputimer)
  call sleep(4)

  call toc(mysystimer)
  call toc(mycputimer)
  call report(mysystimer)
  call report(mycputimer)
  
end program test_timing


