module  modmixermsec
real(8),allocatable:: residual(:),last_outputp(:),last_inputp(:)
real(8),allocatable::PWHIST(:),FHIST(:),CLMHIST(:),yhist(:)
integer:: record_of_last_iter
integer, parameter::noldstepsmax=8,icond=1 !
real(8)::   scl_plane
contains
subroutine initmixermsec(n)
integer,intent(in)::n
allocate(residual(n),last_outputp(n),last_inputp(n))
allocate(PWHIST(noldstepsmax),FHIST(noldstepsmax),CLMHIST(noldstepsmax),yhist(noldstepsmax))
record_of_last_iter=0
residual=0
last_outputp=0
last_inputp=0
PWHIST=0
FHIST=0
CLMHIST=0
yhist=0
scl_plane=1
end subroutine

 subroutine freearraysmixermsec()
 deallocate(residual,last_outputp,last_inputp)
 deallocate(PWHIST,FHIST,CLMHIST)

end subroutine

end module
