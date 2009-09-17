subroutine propertylauncher
  use modinput
  use modmain,only: task
  use modmpi,only:rank
  implicit none

  if(associated(input%properties%dos) .and. rank .eq. 0) then
     task=10
     call dos()
  endif
  if(associated(input%properties%wfplot) .and. rank .eq. 0) then
#define NOTSTM .false.
     call wfplot(NOTSTM)
  endif
  if(associated(input%properties%STM) .and. rank .eq. 0) then
#define STM .true.
     call wfplot(STM)
  endif
  if(associated(input%properties%LSJ) .and. rank .eq. 0) then
     call writelsj
  endif
  if(associated(input%properties%masstensor) .and. rank .eq. 0) then
     call effmass
  endif
  if(associated(input%properties%chargedensityplot ).and. rank .eq. 0) then
     call rhoplot
  endif

  if(associated(input%properties%exccplot) .and. rank .eq. 0) then
     call potplot
  endif

  if(associated(input%properties%elfplot) .and. rank .eq. 0) then
     call elfplot
  endif

  if(associated(input%properties%xcmvecfield) .and. rank .eq. 0) then

     if(associated(input%properties%xcmvecfield%plot2d) .and. rank .eq. 0)then
    	task=82
    	call vecplot
     endif
     if(associated(input%properties%xcmvecfield%plot3d) .and. rank .eq. 0)then
    	task=83
    	call vecplot
     endif
  endif
  if(associated(input%properties%mvecfield) ) then
     if(associated(input%properties%mvecfield%plot2d) .and. rank .eq. 0)then
    	task=72
    	call vecplot
     endif
     if(associated(input%properties%mvecfield%plot3d) .and. rank .eq. 0)then
    	task=73
    	call vecplot
     endif
  endif

  if(associated(input%properties%electricfield ) )  then
     if(associated(input%properties%electricfield%plot2d) .and. rank .eq. 0)then
    	task=142
    	call vecplot
     endif
     if(associated(input%properties%electricfield%plot3d) .and. rank .eq. 0)then
    	task=143
    	call vecplot
     endif
  endif
  if(associated(input%properties%gradmvecfield) .and. rank .eq. 0) then
     call dbxcplot
  endif

  if(associated(input%properties%EFG) .and. rank .eq. 0) then
     call writeefg
  endif
  if(associated(input%properties%momentummatrix) .and. rank .eq. 0) then
     task=120
     call writepmat
  endif
  if(associated(input%properties%dielectric) .and. rank .eq. 0) then
     task= 121
     call dielectric
  endif
  if(associated(input%properties%mossbauer) .and. rank .eq. 0) then
     call mossbauer
  endif
  if(associated(input%properties%expiqr) .and. rank .eq. 0) then
     call  writeexpiqr
  endif
  if(associated(input%properties%elnes) .and. rank .eq. 0) then
     call elnes
  endif
  if(associated(input%properties%eliashberg) .and. rank .eq. 0) then
     call alpha2f
  endif
   if(associated(input%properties%fermisurfaceplot) .and. rank .eq. 0) then
     if(input%properties%fermisurfaceplot%separate)then
        task=101
        call fermisurf
     else
        task=100
        call fermisurf
     endif
  endif
  if(associated(input%properties%bandstructure) .and. rank .eq. 0) then
     task=20
     call bandstr
  endif
end subroutine propertylauncher
