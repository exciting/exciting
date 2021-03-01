
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

program eos
use modmain
implicit none
call readinput
call getedata(etype,nparam,ename,pname)
call fitdata
call output
stop
end program

