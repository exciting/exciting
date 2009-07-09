subroutine xstasklauncher
use modinput
use modmain,only:task
use inputdom

if(.not.(associated(input%xs%tddft)))then
  ! set the default values if solver element not present
  input%xs%tddft=>getstructtddft(emptynode)
endif
if(.not.(associated(input%xs%screening)))then
  ! set the default values if solver element not present
  input%xs%screening=>getstructscreening(emptynode)
endif
if(.not.(associated(input%xs%BSE)))then
  		! set the default values if solver element not present
  		input%xs%BSE=>getstructBSE(emptynode)
  	endif

if(associated(input%xs%plan)) then
	call xsmain(input%xs%plan)
else if(trim(input%xs%xstype).eq."TDDFT") then
 	if(input%xs%tddft%fxctypenumber.eq.7 .or. input%xs%tddft%fxctypenumber.eq.8) then
   	    task=401
 		call scrgeneigvec
 		task=420
 		call scrwritepmat
 			if(associated(input%xs%tetra)) then
		task=410
		call scrtetcalccw
		endif
		task=430
		call screen
		task=440
		call scrcoulint
		task=450
		call kernxc_bse
    endif
    task=301
    call xsgeneigvec
    task=320
	call writepmatxs
	task=330
	call writeemat
	if(associated(input%xs%tetra)) then
	task=310
		call tetcalccw
	endif
	task=340
	call df
else if(trim(input%xs%xstype).eq."BSE")then

    task=401
 	call scrgeneigvec
 	task=420
 	call scrwritepmat

 	task=441
 	call exccoulint
 	if(associated(input%xs%tetra)) then
		task=410
		call scrtetcalccw
	endif

	task=430
    call screen
    task=440
	call scrcoulint
	task=445
    call bse
else
!error
endif


end subroutine
