module  modmixermsec
real(8),allocatable:: work(:),residual(:),last_outputp(:),last_inputp(:)
contains
subroutine initmixermsec(n)
integer,intent(in)::n
allocate(residual(n),last_outputp(n),last_inputp(n))

residual=0
last_outputp=0
last_inputp=0

end subroutine

 subroutine freearraysmixermsec()
 deallocate(residual,last_outputp,last_inputp)

end subroutine

end module


