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

if(.not.(associated(input%xs%tetra)))then
  ! set the default values if solver element not present
	input%xs%tetra=>getstructtetra(emptynode)
endif

call backup0
call backup1
call backup2

if(associated(input%xs%plan)) then
	call xsmain(input%xs%plan)
else if(trim(input%xs%xstype).eq."TDDFT") then

    task=301
    call xsinit
    call xsgeneigvec
    call xsfinit
    
    if((input%xs%tetra%tetradf)) then
	task=310
	call xsinit
	call tetcalccw
	call xsfinit
    endif
    
    task=320
    call xsinit
    call writepmatxs
    call xsfinit
    
    task=330
    call xsinit
    call writeemat
    call xsfinit
    
    if(input%xs%tddft%fxctypenumber.eq.7 .or. input%xs%tddft%fxctypenumber.eq.8) then
    	task=401
    	call xsinit
    	call scrgeneigvec
    	call xsfinit
	
    	task=420
    	call xsinit
    	call scrwritepmat
    	call xsfinit
	
    	if((input%xs%tetra%tetradf)) then
    	    task=410
    	    call xsinit
    	    call scrtetcalccw
    	    call xsfinit
    	endif
	
    	task=430
    	call xsinit
    	call screen
    	call xsfinit
	
    	task=440
    	call xsinit
    	call scrcoulint
    	call xsfinit
	
    	task=450
    	call xsinit
    	call kernxc_bse
    	call xsfinit
    endif

    task=340
    call xsinit
    call df
    call xsfinit
    
    task=350
    call xsinit
    call idf
    call xsfinit
    
else if(trim(input%xs%xstype).eq."BSE")then

    task=401
    call xsinit
    call scrgeneigvec
    call xsfinit
    
    if((input%xs%tetra%tetradf)) then
    	    task=410
    	    call xsinit
    	    call scrtetcalccw
            call xsfinit
    endif

    task=420
    call xsinit
    call scrwritepmat
    call xsfinit
    
    task=430
    call xsinit
    call screen
    call xsfinit
    
    task=440
    call xsinit
    call scrcoulint
    call xsfinit
    
    task=441
    call xsinit
    call exccoulint
    call xsfinit
    
    task=445
    call xsinit
    call bse
    call xsfinit
else
  write(*,*)"error xstasklauncher"
  write(*,*)trim(input%xs%xstype),"no valid xstype"
  stop
endif


end subroutine
