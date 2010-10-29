
! Copyright (C) 2009-2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine groundstatetasklauncher
      Use modinput
      Use modmain, Only: task,xctype
      Use inputdom
      Implicit None
      xctype(1)=  input%groundstate%xctypenumber
	!    if we use libxc then
       if(associated(input%groundstate%libxc)) then
	       input%groundstate%xctypenumber=100
	       xctype(1)=100
	       input%groundstate%xctype="LibXC"
	       xctype(2)=input%groundstate%libxc%exchangenumber
	       xctype(3)=input%groundstate%libxc%correlationnumber
           if (input%groundstate%libxc%xcnumber .ne. 0)then
				xctype(2)=input%groundstate%libxc%xcnumber
				xctype(3)=0
		   endif
       endif
      If ( .Not. (associated(input%groundstate%solver))) Then
        ! set the default values if solver element not present
         input%groundstate%solver => getstructsolver (emptynode)
      End If
      If (input%groundstate%do .Eq. "fromscratch") Then
         If (associated(input%structureoptimization)) Then
            task = 2
         Else
            task = 0
         End If
      Else
         If (associated(input%structureoptimization)) Then
            task = 3
         Else
            task = 1
         End If
      End If
      If (input%groundstate%do .Ne. "skip") then
        Call gndstate
        ! do conversion to XML format if requested
        if (associated(input%groundstate%output)) then
          if (input%groundstate%output%state .eq. "XML") call portstate(1)
        end if
      end if

end subroutine
