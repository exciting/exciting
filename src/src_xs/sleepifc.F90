
subroutine sleepifc(sec)
  implicit none
  ! arguments
  integer, intent(in) :: sec
  ! local variables
  integer, parameter :: sleepfac=200000
  integer :: i,j,k

  ! interface to sleep routine
  ! use the ifort intrinsic for the moment
  !call sleep(sec)

  ! mathematical sleep implementation
  ! depending on the the performance of the cpu
  ! and on the load
  call sleepm(sec*sleepfac)

contains

  subroutine sleepm(c)
    implicit none
    ! arguments
    integer, intent(in) :: c
    ! local variables
    character(10) dat, tim
    integer :: i,sys0,sys1,cnt
    complex(8) :: m(10,10), n(10,10)

!!$    call system_clock(COUNT_RATE=cnt)
!!$    call system_clock(COUNT=sys0)

    do i=1,c
       n=matmul(m,m)
    end do

!!$    call system_clock(COUNT=sys1)
!!$    write(*,*) 'Time elapsed   : ', dble(sys1-sys0)/dble(cnt)

  end subroutine sleepm

end subroutine sleepifc

