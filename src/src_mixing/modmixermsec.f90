

module	modmixermsec

real(8), allocatable:: residual(:), last_outputp(:), work2(:), work3(:)
real(8), allocatable::PWHIST(:), FHIST(:), CLMHIST(:), yhist(:)
integer:: record_of_last_iter, noldstepsin_file, noldsteps, MUSE,  IDSCALE
integer, parameter::icond=1, noldstepsmax=8, dbase=0.005D0
real(8)::   scl_plane, qmx, RedOld, RedPred, qmx_input, PM1, DIAG, dmix_last, dmixout(4)
real(8):: MSECINFO(20), rtrap, SCHARGE, TCharge, splane, tplane, qtot
real(8)::dmix
contains


subroutine initmixermsec(n)
use modmain, only:CHGIR, CHGMTTOT
integer, intent(in)::n

integer::niter
allocate(residual(n), last_outputp(n), work2(n), work3(n))

allocate(PWHIST(noldstepsmax), FHIST(noldstepsmax), CLMHIST(noldstepsmax), yhist(noldstepsmax))
record_of_last_iter=0
residual=0
last_outputp=0
work2=0
work3=0
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
MSECINFO=1

DMIX=.5
end subroutine


subroutine freearraysmixermsec()
 character(256), external:: outfilenamestring
	character(256)::filetag
	filetag="BROYDEN"
 deallocate(residual, last_outputp)
 if (allocated(work2)) deallocate(work2)
 if (allocated(work3)) deallocate(work3)
 deallocate(PWHIST, FHIST, CLMHIST, yhist)
 open(23, file=outfilenamestring(filetag, 1))
 close(23, STATUS='DELETE')
end subroutine

end module
