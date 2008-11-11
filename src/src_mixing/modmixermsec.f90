module  modmixermsec
use modmain,only:beta0
real(8),allocatable:: residual(:),last_outputp(:),last_inputp(:),last_outputp2(:)
real(8),allocatable::PWHIST(:),FHIST(:),CLMHIST(:),yhist(:)
integer:: record_of_last_iter,noldstepsin_file,noldsteps,MUSE,  IDSCALE
integer, parameter::icond=1,noldstepsmax=8,dbase=0.005D0
real(8)::   scl_plane,qmx,RedOld,RedPred,qmx_input,PM1,DIAG,dmix_last,dmixout(4)
real(8):: MSECINFO(20),rtrap,SCHARGE,TCharge,splane,tplane,qtot
real(8)::dmix
contains
subroutine initmixermsec(n)
integer,intent(in)::n

integer::niter
allocate(residual(n),last_outputp(n),last_outputp2(n),last_inputp(n))
allocate(PWHIST(noldstepsmax),FHIST(noldstepsmax),CLMHIST(noldstepsmax),yhist(noldstepsmax))
record_of_last_iter=0
residual=0
last_outputp=0
last_inputp=0
last_outputp2=0
PWHIST=0
FHIST=0
CLMHIST=0
yhist=0
scl_plane=4
RedOld=1
RedPred=1
qmx_input=.2
qmx=qmx_input
PM1=1
IDSCALE=1
DIAG =5D-4
noldstepsin_file=0
noldsteps=0
rtrap =0.1
SCHARGE=chgir
TCharge=chgmttot
splane=.000001
tplane=.000001
MSECINFO=0

DMIX=.5
end subroutine

 subroutine freearraysmixermsec()
 deallocate(residual,last_outputp,last_inputp)
 deallocate(PWHIST,FHIST,CLMHIST)
 open(23,file="BROYDEN.OUT")
 close(23,STATUS='DELETE')
end subroutine

end module
