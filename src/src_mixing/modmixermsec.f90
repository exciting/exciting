module  modmixermsec
real(8),allocatable:: work(:),residual(:),last_outputp(:),last_inputp(:)
real(8),allocatable::PWHIST(:),FHIST(:),CLMHIST(:),yhist(:)
integer:: record_of_last_iter
integer, parameter::noldstepsmax=8
contains
subroutine initmixermsec(n)
integer,intent(in)::n
allocate(residual(n),last_outputp(n),last_inputp(n))
allocate(PWHIST(noldstepsmax),FHIST(noldstepsmax),CLMHIST(noldstepsmax),yhist(noldstepsmax))
record_of_last_iter=0
residual=0
last_outputp=0
last_inputp=0

end subroutine

 subroutine freearraysmixermsec()
 deallocate(residual,last_outputp,last_inputp)
 deallocate(PWHIST,FHIST,CLMHIST)

end subroutine

end module


