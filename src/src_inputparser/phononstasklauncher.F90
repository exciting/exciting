
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine phononstasklauncher
      Use modinput
      Use modmain, Only: task
      use dfpt, only: dfpt_launcher
      Use inputdom
!
      Implicit None
      character(512) :: part
      integer :: i

      if( input%phonons%method == 'sc' .and. associated( input%phonons%parts)) then

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
        if( input%phonons%method == 'dfpt') then
          if( input%phonons%do /= 'skip') then
            call dfpt_launcher
          end if
        else
          if( input%phonons%do /= 'skip') then
            Call phonon
            ! restore input parameters manipulated by supercell setup
            call rereadinput
            call writeborncharges
          end if
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
