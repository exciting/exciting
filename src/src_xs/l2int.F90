integer function l2int(l)
  implicit none
  logical, intent(in) :: l
  if (l) then
     l2int=1
  else
     l2int=0
  end if
end function l2int
