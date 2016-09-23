
! Copyright (C) 2010 C. Meisenbichler, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine initlibxc
  Use modinput
  Use modmain, Only: xctype
  Implicit None
    xctype(1)=  input%groundstate%xctypenumber
!    if we use libxc then
   if (associated(input%groundstate%libxc)) then
       input%groundstate%xctypenumber=100
       xctype(1)=100
       input%groundstate%xctype="LibXC"
       xctype(2)=input%groundstate%libxc%exchangenumber
       xctype(3)=input%groundstate%libxc%correlationnumber
       if (input%groundstate%libxc%xcnumber .ne. 0)then
            xctype(2)=input%groundstate%libxc%xcnumber
            xctype(3)=0
       end if
   end if
end subroutine
