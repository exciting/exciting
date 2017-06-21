subroutine testingfun
  use modinput
  use modmpi
  use mod_ematgrids
  use mod_misc, only: task

  implicit none

  write(*,*) "Hello, testing fun here!"

  call init0
  call init1

  call init2



  call terminate

end subroutine testingfun
