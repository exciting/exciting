
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Module m_gndstateq
      Implicit None
Contains
      Subroutine gndstateq (voff, filxt)
         Use modmain
         Use modinput
          Use modxs, only: skipgnd
         Implicit None
    ! arguments
         Real (8), Intent (In) :: voff (:)
         Character (*), Intent (In) :: filxt
    ! local varialbes
         Character (256) :: filext_save
         Real (8) :: vkloff_save (3)
!         Real (8) ::  bfieldc_save(3)
         Integer :: task_save
!         Integer :: maxscl_save
    ! save original values
         filext_save = trim (filext)
         vkloff_save = input%groundstate%vkloff
!         If (associated(input%groundstate%spin)) Then
!            bfieldc_save = input%groundstate%spin%bfieldc
!         End If
         task_save = task
 !        maxscl_save = input%groundstate%maxscl
    ! one iteration, new offset, special file extension
         filext = trim (filxt)
         input%groundstate%vkloff = voff
         task = 1

         !ran gs from scratch if spin polarized calculation is needed.
         If (input%xs%dogroundstate .Eq. "fromscratch") Then
           task=0
        End If
    ! call with the above parameters changed
        If (.not.(skipgnd)) Then 
               Call gndstate
        End If  
         Call rewritesorted
    ! restore original parameters
         If (input%xs%dogroundstate .Ne. "fromscratch") Then
            filext = trim (filext_save)
         End If
         input%groundstate%vkloff = vkloff_save
  !       If (associated(input%groundstate%spin)) Then
  !          input%groundstate%spin%bfieldc = bfieldc_save
  !       End If
         task = task_save
!         input%groundstate%maxscl = maxscl_save
      End Subroutine gndstateq
End Module m_gndstateq
