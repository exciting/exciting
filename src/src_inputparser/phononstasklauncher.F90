
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine phononstasklauncher
      Use modinput
      Use modmain, Only: task
      Use inputdom
!
      Implicit None
      character(512) :: part
      integer :: i

      If (associated(input%phonons%parts)) then
       Do i = 1, size (input%phonons%parts%dopartarray)
         part = trim(adjustl(input%phonons%parts%dopartarray(i)%dopart%id))
         write(*,'("Info(phononstasklauncher): executing part Nr. ",i6," : ",a)') i, trim(part)
         Select Case (trim(part))
         Case ('reformatdynamicalmatrices')
            task=200
            call reformatdynamicalmatrices()
         case ('writephn')
            task=230
            call writephn()
         case ('debug')
            task=200
            call phonondebug()
         case default
            Write (*,*)
            Write (*,*) 'Error(phononstasklauncher): id not defined: ', trim(part)
            Write (*,*)
            stop
         end select
       end do
      else
        task=200
        ! task 201 is only a dry-run and will not be considered here
        If (input%phonons%do .Ne. "skip") then
            Call phonon
            ! restore input parameters manipulated by supercell setup
            call rereadinput
        end if
        if (associated(input%phonons%reformatdynmat)) then
            call reformatdynamicalmatrices
        end if
        if (associated(input%phonons%qpointset)) then
            task=230
            call writephn
        end if
        if (associated(input%phonons%phonondos)) then
            task=210
            call phdos
        end if
        if (associated(input%phonons%phonondispplot)) then
            task=220
            call phdisp
        end if
        if (associated(input%phonons%interpolate)) then
            call phononinterpolate
        end if

      end if

End Subroutine
