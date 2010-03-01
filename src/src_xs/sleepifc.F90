!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine sleepifc (sec)
      Implicit None
  ! arguments
      Integer, Intent (In) :: sec
  ! local variables
      Integer, Parameter :: sleepfac = 200000
!
  ! interface to sleep routine
  ! use the ifort intrinsic for the moment
  !call sleep(sec)
!
  ! mathematical sleep implementation
  ! depending on the the performance of the cpu
  ! and on the load
      Call sleepm (sec*sleepfac)
!
Contains
!
!
      Subroutine sleepm (c)
         Implicit None
    ! arguments
         Integer, Intent (In) :: c
    ! local variables
         Integer :: i
         Real (8) :: m (10, 10), n (10, 10)
         call random_number(m)
         Do i = 1, c
            n = matmul (m, m)
         End Do
      End Subroutine sleepm
!
End Subroutine sleepifc
