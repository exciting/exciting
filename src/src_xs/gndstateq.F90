
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Module m_gndstateq
      Implicit None
Contains
      Subroutine gndstateq (voff, filxt)
         Use modmain
         Use modinput
         Implicit None
    ! arguments
         Real (8), Intent (In) :: voff (:)
         Character (*), Intent (In) :: filxt
    ! local varialbes
         Character (256) :: filext_save
         Real (8) :: vkloff_save (3)
         Integer :: task_save, maxscl_save
    ! save original values
         filext_save = trim (filext)
         vkloff_save = input%groundstate%vkloff
         task_save = task
         maxscl_save = input%groundstate%maxscl
    ! one iteration, new offset, special file extension
         filext = trim (filxt)
         input%groundstate%vkloff = voff
         task = 1
         input%groundstate%maxscl = 1
    ! call with the above parameters changed
         Call gndstate
    ! restore original parameters
         filext = trim (filext_save)
         input%groundstate%vkloff = vkloff_save
         task = task_save
         input%groundstate%maxscl = maxscl_save
      End Subroutine gndstateq
End Module m_gndstateq
