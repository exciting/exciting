




! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readinput
! !INTERFACE:


subroutine setdefault
! !USES:
use modinput
use modmain
#ifdef TETRA
use modtetra
#endif
#ifdef XS
use modmpi, only: rank
use modxs
#endif

use sclcontroll
! !DESCRIPTION:
!   Reads in the input parameters from the file {\tt exciting.in}. Also sets
!   default values for the input parameters.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!   Additional parmeters for excited states and tetrahedron method
!     2004-2008 (Sagmeister)
!EOP
!BOC
implicit none
! local variables
integer::is, js, ia, ja, ias
integer::i, l, iv, iostat
!!real(8) sc,sc1,sc2,sc3
!real(8) vacuum,v(3),t1,t2
character(256)::bname

#ifdef XS
character(256) :: fname
logical, parameter :: dumpmain=.true.
logical, parameter :: dumpadd=.true.
logical, parameter :: dumptetra=.true.
logical, parameter :: dumpmpiiter=.true.
logical, parameter :: dumpxs=.true.
#endif

!------------------------!
!     default values     !
!------------------------!
ntasks=1
tasks(1)=-1
!avec(:,:)=0.d0
!avec(1,1)=1.d0
!avec(2,2)=1.d0
!avec(3,3)=1.d0
!sc=1.d0
!sc1=1.d0
!sc2=1.d0
!sc3=1.d0
!epslat=1.d-6
!primcell=.false.
tshift=.true.
!ngridk(:)=1
!vkloff(:)=0.d0
!autokpt=.false.
!radkpt=40.0
!reducek=.true.
!ngridq(:)=1
!reduceq=.true.
!rgkmax=7.d0
!gmaxvr=12.d0
!lmaxapw=8
!lmaxvr=7
!lmaxmat=5
!lmaxinr=2
!fracinr=0.25d0
!npsden=9
!xctype=3
!input%groundstate%stypenumber=0
!input%groundstate%swidth=0.01d0
!epsocc=1.d-8
!epschg=1.d-3
!nempty=5
!maxscl=200
!mixtype=1
beta0=0.4d0
betainc=1.1d0
betadec=0.6d0
!epspotatom=1.d-6
epsengy=1.d-7
!epsforce=5.d-4
cfdamp=0.d0
!molecule=.false.
!vacuum=10.d0
nspecies=0
epsedirac=1.d-11
epspotatom=1.d-6
!atposc(:,:,:)=0.d0
!bfcmt(:,:,:)=0.d0
bflmt(:, :, :)=0.d0
!sppath='./'
scrpath='./'
nvp1d=2
!iterativetype
 tarpack=.false.
 tlapack=.true.
 tdiis=.false.
 tjdqz=.false.
 diisfirstscl=10
lowesteval=-1.d0
packedmatrixstorage=.true.
epsarpack=1e-8
!epsresid=1e-12
maxncv=200
if (allocated(vvlp1d)) deallocate(vvlp1d)
allocate(vvlp1d(3, nvp1d))
vvlp1d(:, 1)=0.d0
vvlp1d(:, 2)=1.d0
npp1d=200
vclp2d(:, :)=0.d0
vclp2d(1, 2)=1.d0
vclp2d(2, 3)=1.d0
np2d(:)=40
vclp3d(:, :)=0.d0
vclp3d(1, 2)=1.d0
vclp3d(2, 3)=1.d0
vclp3d(3, 4)=1.d0
np3d(:)=20
!nwdos=500
!ngrdos=100
!nsmdos=0
wdos(1)=-0.5d0
wdos(2)=0.5d0
!lmirep=.false.
!spinpol=.false.
!spinorb=.false.
!tau0atm=0.2d0
!nstfsp=6
!lradstp=4
!chgexs=0.d0
!nprad=4
!!!scissor=0.d0
!noptcomp=1
!optcomp(:,1)=1
!<sag>
!optswidth=0.d0
!</sag>
usegdft=.false.
intraband=.false.
evaltol=1.d-8
!evalmin=-4.5d0
!deband=0.0025d0
!bfieldc(:)=0.d0
!fixspin=0
!momfix(:)=0.d0
mommtfix(:, :, :)=0.d0
!taufsm=0.01d0
!autormt=.false.
!rmtapm(1)=0.25d0
!rmtapm(2)=0.95d0
!isgkmax=-1
!nosym=.false.
!deltaph=0.03d0
nphwrt=1
if (allocated(vqlwrt)) deallocate(vqlwrt)
allocate(vqlwrt(3, nphwrt))
vqlwrt(:, :)=0.d0
notelns=0
!!tforce=.false.
!!tfibs=.true.
maxitoep=120
tauoep(1)=1.d0
tauoep(2)=0.2d0
tauoep(3)=1.5d0
nkstlist=1
kstlist(:, 1)=1
!vklem(:)=0.d0
!deltaem=0.025d0
!ndspem=1
nosource=.false.
!spinsprl=.false.
!vqlss(:)=0.d0
nwrite=0
tevecsv=.false.
ldapu=0
llu(:)=-1
ujlu(:, :)=0.d0
rdmxctype=2
rdmmaxscl=1
maxitn=250
maxitc=10
taurdmn=1.d0
taurdmc=0.5d0
rdmalpha=0.7d0
rdmtemp=0.d0
!reducebf=1.d0
ptnucl=.true.
tseqit=.false.
nseqit=40
tauseq=0.1d0
vecql(:)=0.d0
mustar=0.15d0
sqados(1:2)=0.d0
sqados(3)=1.d0
#ifdef TETRA
! tetrahedron method variables
tetraocc=.false.
tetraopt=.false.
tetradf=.false.
tetrakordexc=.false.
#endif
#ifdef XS
! TDDFT variables
imbandstr=.false.
nqptmt=1
if (allocated(vgqlmt)) deallocate(vgqlmt)
allocate(vgqlmt(3, nqptmt))
vgqlmt(:, :)=0.d0
mdfqtype=0
vqloff(:)=0.d0
tq0ev=.true.
gqmax=0.d0
lmaxapwwf=-1
fastpmat=.true.
fastemat=.true.
emattype=1
lmaxemat=3
rsptype='reta'
acont=.false.
nwacont=0
broad=0.01d0
aresdf=.true.
epsdfde=1.d-8
emaxdf=1.d10
symwings=.false.
lfediag=.false.
!fxctype=0
nexcitmax=100
alphalrc=0.d0
alphalrcdyn=0.d0
betalrcdyn=0.d0
ndftrans=1
if (allocated(dftrans)) deallocate(dftrans)
allocate(dftrans(3, ndftrans))
dftrans(:, :)=0
tetraqweights=1
gather=.false.
symmorph=.false.
tevout=.false.
tappinfo=.false.
dbglev=0
! screening variables
screentype='full'
nosymscr=.false.
reducekscr=.true.
ngridkscr(:)=-1
vkloffscr(:)=-1.d0
rgkmaxscr=-1.d0
nemptyscr=-1
scrherm=0
fnevecfvscr='EVECFV_SCR.OUT'
fnevalsvscr='EVALSV_SCR.OUT'
fnoccsvscr='OCCSV_SCR.OUT'
! BSE (-kernel) variables
bsetype='ip'
nosymbse=.false.
reducekbse=.true.
vkloffbse(:)=-1.d0
bsediagweight=1
bsediagsym=0
fnevecfvbse='EVECFV_BSE.OUT'
fnevalsvbse='EVALSV_BSE.OUT'
fnoccsvbse='OCCSV_BSE.OUT'
nbfce=-1
nafce=-1
nbfbse=-1
nafbse=-1
! dump default parameters
if (rank.eq.0) then
   fname='PARAMS_DEFAULT.OUT'
   if (dumpmain) call dumpparams(trim(fname),  main parameters:', input%structure%speciespath, &
    &input%structure%crystal%scale, input%structure%crystal%stretch, input%structure%crystal%stretch, &
    &input%structure%crystal%stretch, input%structure%vacuum)
   if (dumpadd) call dumpparams_add(trim(fname),       '! additional parameters:')
   if (dumpmpiiter) call dumpparams_mpiiter(trim(fname),	'! MPI parallelization and iterative&
    &input%groundstate%solvernumber parameters:')
   if (dumptetra) call dumpparams_tetra(trim(fname),	     '! tetrahedron method parameters:')
   if (dumpxs) call dumpparams_xs(trim(fname),	       '! excited states parameters:')
end if
#endif

end subroutine
!EOC
