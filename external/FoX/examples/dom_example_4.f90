program dom_example_4

  ! Example of error handling from the DOM interface
  !
  ! Note the need to check the value of the iostat as well
  ! as ex. If the iostat optional argument is present and
  ! the file cannot be found that is the variable that gets
  ! set in the parseFile function. If iostat is not present 
  ! in this case FoX will abort itself reporting that it 
  ! "Cannot open file". This kind of error does not get 
  ! attached to the DOMException.
  ! 
  ! If the file exists but is not well-formed XML the error
  ! should be reported in the DOMException. One oddity to note 
  ! is that if ex is not present but ios is, you may end up 
  ! with a FoX error code in ios. That only happens if the file
  ! exists but it not well-formed XML.

  use FoX_dom
  implicit none

  type(Node), pointer :: myDoc
  type(DOMException) :: ex
  integer :: ios

  myDoc => parseFile("h2o.xml", iostat=ios, ex=ex)
  if (ios.ne.0) then
      print*, "Problem reading file. iostat was ", ios
      stop

  elseif (inException(ex)) then
      print*,"DOM Parse error ", getExceptionCode(ex)
      stop

  else
      print*, "Not in exception and file read OK"

  endif

  ! Do something with myDoc

end program dom_example_4
