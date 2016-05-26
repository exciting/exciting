
subroutine gwtasklauncher()
   
  use modinput
  use modmain, only: task
  use inputdom

  if (associated(input%gw)) then
     call rereadinput
     task = 1
     call gw_main
  else
     write (*,*) "error gwtasklauncher"
     stop
  end if

  return
end subroutine gwtasklauncher
