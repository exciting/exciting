
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine sleepifc(sec)
  implicit none
  ! arguments
  integer, intent(in) :: sec
  ! local variables
  integer, parameter :: sleepfac=200000

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
    integer :: i
!!$    integer :: sys0,sys1,cnt
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

