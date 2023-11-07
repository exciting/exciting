
! Copyright (C) 2009-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!> Select exciting features to run
subroutine tasklauncher()
    use modinput
    use xs_task_launcher
    
    implicit none

    ! Note that the order of the calls below may be important!
    if (associated(input%groundstate)) &
        call groundstatetasklauncher()

    if (associated(input%properties)) &
        call propertylauncher()

    if (associated(input%phonons)) &
        call phononstasklauncher()

    if (associated(input%gw)) &
        call gwtasklauncher()

    if (associated(input%xs)) then
        call xstasklauncher()
    endif 

end subroutine
