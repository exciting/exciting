module  modmixermsec
real(8),allocatable:: residual(:),last_outputp(:),last_inputp(:)
real(8),allocatable::PWHIST(:),FHIST(:),CLMHIST(:),yhist(:)
integer:: record_of_last_iter,noldstepsin_file,noldsteps,MUSE,  IDSCALE
integer, parameter::icond=1,noldstepsmax=8
real(8)::   scl_plane,qmx,RedOld,RedPred,qmx_input,PM1,DIAG,dmix_last,dmixout(4)
real*8 MSECINFO(20)
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
RedOld=1
RedPred=1
qmx_input=.5
qmx=qmx_input
PM1=0
IDSCALE=0
DIAG =5D-4
noldstepsin_file=0
noldsteps=0

end subroutine

 subroutine freearraysmixermsec()
 deallocate(residual,last_outputp,last_inputp)
 deallocate(PWHIST,FHIST,CLMHIST)

end subroutine

end module
