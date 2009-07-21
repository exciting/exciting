! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! main routine for the EXCITING code

subroutine tasklauncher
use modinput
use modmain,only:task
use inputdom

implicit none

  if(associated(input%groundstate).and. .not. input%groundstate%do .eq. "skipp") then
  	if(.not.(associated(input%groundstate%solver)))then
  		! set the default values if solver element not present
  		input%groundstate%solver=>getstructsolver(emptynode)
  	endif
   	 if(input%groundstate%do.eq."fromscratch") then
	    if(associated(input%structureoptimization)) then
	        task=2
	 	else
	        task=0
	   	endif
  	else
	    if(associated(input%structureoptimization)) then
	 	    task=3
	    else
	   		task=1
	   	endif
    endif
  	call gndstate
  endif


if(associated(input%properties)) then
   	call propertylauncher
endif

if(associated(input%phonons)) then
	call phonon
endif

if(associated(input%xs)) then
	call xstasklauncher()
endif
end subroutine
