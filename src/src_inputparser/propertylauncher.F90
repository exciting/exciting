subroutine propertylauncher
  use modinput
  use modmain,only: task
  implicit none

  if(associated(input%properties%dos)) then
     task=10
     call dos()
  endif
  if(associated(input%properties%wfplot)) then
#define NOTSTM .false.
     call wfplot(NOTSTM)
  endif
  if(associated(input%properties%STM)) then
#define STM .true.
     call wfplot(STM)
  endif
  if(associated(input%properties%LSJ)) then
     call writelsj
  endif
  if(associated(input%properties%masstensor)) then
     call effmass
  endif
  if(associated(input%properties%chargedesityplot)) then
     call rhoplot
  endif

  if(associated(input%properties%exccplot)) then
     call potplot
  endif

  if(associated(input%properties%elfplot)) then
     call elfplot
  endif

  if(associated(input%properties%xcmvecfield)) then

     if(associated(input%properties%xcmvecfield%plot2d))then
    	task=82
    	call vecplot
     endif
     if(associated(input%properties%xcmvecfield%plot3d))then
    	task=83
    	call vecplot
     endif
  endif
  if(associated(input%properties%mvecfield) ) then
     if(associated(input%properties%mvecfield%plot2d))then
    	task=72
    	call vecplot
     endif
     if(associated(input%properties%mvecfield%plot3d))then
    	task=73
    	call vecplot
     endif
  endif

  if(associated(input%properties%electricfield ) )  then
     if(associated(input%properties%electricfield%plot2d))then
    	task=142
    	call vecplot
     endif
     if(associated(input%properties%electricfield%plot3d))then
    	task=143
    	call vecplot
     endif
  endif
  if(associated(input%properties%gradmvecfield)) then
     call dbxcplot
  endif

  if(associated(input%properties%EFG)) then
     call writeefg
  endif
  if(associated(input%properties%momentummatrix)) then
     task=120
     call writepmat
  endif
  if(associated(input%properties%dielectric)) then
     task= 121
     call dielectric
  endif
  if(associated(input%properties%mossbauer)) then
     call mossbauer
  endif
  if(associated(input%properties%expiqr)) then
     call  writeexpiqr
  endif
  if(associated(input%properties%elnes)) then
     call elnes
  endif
  if(associated(input%properties%eliashberg)) then
     call alpha2f
  endif
   if(associated(input%properties%fermisurfaceplot)) then
     if(input%properties%fermisurfaceplot%separate)then
        task=101
        call fermisurf
     else
        task=100
        call fermisurf
     endif
  endif
  if(associated(input%properties%bandstructure)) then
     task=20
     call bandstr
  endif
end subroutine propertylauncher
