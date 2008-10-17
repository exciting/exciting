
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readinput
! !INTERFACE:
subroutine readinput
! !USES:
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
!   Reads in the input parameters from the file {\tt exciting.in} as well as
!   from the species files. Also sets default values for the input parameters.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!   Additional parmeters for excited states and tetrahedron method
!     2004-2008 (Sagmeister)
!EOP
!BOC
implicit none
! local variables
integer is,js,ia,ja,ias
integer i,l,iv,iostat
integer ist,io,nlx,ilx,lx,ilo
real(8) sc,sc1,sc2,sc3
real(8) vacuum,v(3),t1,t2
character(256) str
character(256) bname
character(256) sppath
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
avec(:,:)=0.d0
avec(1,1)=1.d0
avec(2,2)=1.d0
avec(3,3)=1.d0
sc=1.d0
sc1=1.d0
sc2=1.d0
sc3=1.d0
epslat=1.d-6
primcell=.false.
tshift=.true.
ngridk(:)=1
vkloff(:)=0.d0
autokpt=.false.
radkpt=40.0
reducek=.true.
ngridq(:)=1
reduceq=.true.
rgkmax=7.d0
gmaxvr=12.d0
lmaxapw=8
lmaxvr=7
lmaxmat=5
lmaxinr=2
fracinr=0.25d0
npsden=9
xctype=3
stype=0
swidth=0.01d0
epsocc=1.d-8
epschg=1.d-3
nempty=5
maxscl=200
mixtype=1
beta0=0.4d0
betainc=1.1d0
betadec=0.6d0
epspot=1.d-6
epsengy=1.d-7
epsforce=5.d-4
cfdamp=0.d0
molecule=.false.
vacuum=10.d0
nspecies=0
natoms(:)=0
atposl(:,:,:)=0.d0
atposc(:,:,:)=0.d0
bfcmt(:,:,:)=0.d0
bflmt(:,:,:)=0.d0
sppath='./'
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
epsresid=1e-12
maxncv=200
if (allocated(vvlp1d)) deallocate(vvlp1d)
allocate(vvlp1d(3,nvp1d))
vvlp1d(:,1)=0.d0
vvlp1d(:,2)=1.d0
npp1d=200
vclp2d(:,:)=0.d0
vclp2d(1,2)=1.d0
vclp2d(2,3)=1.d0
np2d(:)=40
nup3d(:)=1
np3d(:)=20
nwdos=500
ngrdos=100
nsmdos=0
wdos(1)=-0.5d0
wdos(2)=0.5d0
lmirep=.false.
spinpol=.false.
spinorb=.false.
tau0atm=0.2d0
nstfsp=6
lradstp=4
chgexs=0.d0
nprad=4
scissor=0.d0
noptcomp=1
optcomp(:,1)=1
!<sag>
optswidth=0.d0
!</sag>
usegdft=.false.
intraband=.false.
evaltol=1.d-8
evalmin=-4.5d0
deband=0.0025d0
bfieldc(:)=0.d0
fixspin=0
momfix(:)=0.d0
mommtfix(:,:,:)=0.d0
taufsm=0.01d0
autormt=.false.
rmtapm(1)=0.25d0
rmtapm(2)=0.95d0
isgkmax=-1
nosym=.false.
deltaph=0.03d0
nphwrt=1
if (allocated(vqlwrt)) deallocate(vqlwrt)
allocate(vqlwrt(3,nphwrt))
vqlwrt(:,:)=0.d0
notelns=0
tforce=.false.
tfibs=.true.
maxitoep=120
tauoep(1)=1.d0
tauoep(2)=0.2d0
tauoep(3)=1.5d0
nkstlist=1
kstlist(:,1)=1
vklem(:)=0.d0
deltaem=0.025d0
ndspem=1
nosource=.false.
spinsprl=.false.
vqlss(:)=0.d0
nwrite=0
tevecsv=.false.
ldapu=0
llu(:)=-1
ujlu(:,:)=0.d0
rdmxctype=2
rdmmaxscl=1
maxitn=250
maxitc=10
taurdmn=1.d0
taurdmc=0.5d0
rdmalpha=0.7d0
rdmtemp=0.d0
reducebf=1.d0
ptnucl=.true.
tseqit=.false.
nseqit=6
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
allocate(vgqlmt(3,nqptmt))
vgqlmt(:,:)=0.d0
mdfqtype=0
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
dfoffdiag=.false.
kerndiag=.false.
fxctype=0
nexcitmax=100
alphalrc=0.d0
alphalrcdyn=0.d0
betalrcdyn=0.d0
ndftrans=1
if (allocated(dftrans)) deallocate(dftrans)
allocate(dftrans(3,ndftrans))
dftrans(:,:)=0
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
ngridkscr(:)=0
vkloffscr(:)=-1.d0
rgkmaxscr=-1.d0
nemptyscr=0
scrherm=0
! BSE (-kernel) variables
bsetype='ip'
nosymbse=.false.
reducekbse=.true.
vkloffbse(:)=-1.d0
bsediagweight=1
bsediagsym=0
nbfce=0
nafce=0
nbfbse=0
nafbse=0
! dump default parameters
if (rank.eq.0) then
   fname='PARAMS_DEFAULT.OUT'
   if (dumpmain) call dumpparams(trim(fname), &
        '! main parameters:',sppath,sc,sc1,sc2,sc3,vacuum)
   if (dumpadd) call dumpparams_add(trim(fname), &
        '! additional parameters:')
   if (dumpmpiiter) call dumpparams_mpiiter(trim(fname), &
        '! MPI parallelization and iterative solver parameters:')
   if (dumptetra) call dumpparams_tetra(trim(fname), &
        '! tetrahedron method parameters:')
   if (dumpxs) call dumpparams_xs(trim(fname), &
        '! excited states parameters:')
end if
#endif

!-------------------------------!
!     read from exciting.in     !
!-------------------------------!
open(50,file='exciting.in',action='READ',status='OLD',form='FORMATTED', &
 iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readinput): error opening exciting.in")')
  write(*,*)
  stop
end if
10 continue
read(50,*,end=30) bname
! check for a comment
if ((scan(trim(bname),'!').eq.1).or.(scan(trim(bname),'#').eq.1)) goto 10
select case(trim(bname))
case('tasks')
  do i=1,maxtasks
    read(50,'(A80)') str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): no tasks to perform")')
        write(*,*)
        stop
      end if
      ntasks=i-1
      goto 10
    end if
    read(str,*,iostat=iostat) tasks(i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading tasks")')
      write(*,'("(blank line required after tasks block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): too many tasks")')
  write(*,*)
  stop
case('avec')
  read(50,*,err=20) avec(:,1)
  read(50,*,err=20) avec(:,2)
  read(50,*,err=20) avec(:,3)
case('scale')
  read(50,*,err=20) sc
case('scale1')
  read(50,*,err=20) sc1
case('scale2')
  read(50,*,err=20) sc2
case('scale3')
  read(50,*,err=20) sc3
case('epslat')
  read(50,*,err=20) epslat
  if (epslat.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epslat <= 0 : ",G18.10)') epslat
    write(*,*)
    stop
  end if
case('primcell')
  read(50,*,err=20) primcell
case('tshift')
  read(50,*,err=20) tshift
case('autokpt')
  read(50,*,err=20) autokpt
case('radkpt')
  read(50,*,err=20) radkpt
  if (radkpt.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): radkpt <= 0 : ",G18.10)') radkpt
    write(*,*)
    stop
  end if
case('ngridk')
  read(50,*,err=20) ngridk(:)
  if ((ngridk(1).le.0).or.(ngridk(2).le.0).or.(ngridk(3).le.0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid ngridk : ",3I8)') ngridk
    write(*,*)
    stop
  end if
case('vkloff')
  read(50,*,err=20) vkloff(:)
case('reducek')
  read(50,*,err=20) reducek
case('ngridq')
  read(50,*,err=20) ngridq(:)
  if ((ngridq(1).le.0).or.(ngridq(2).le.0).or.(ngridq(3).le.0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid ngridq : ",3I8)') ngridq
    write(*,*)
    stop
  end if
case('reduceq')
  read(50,*,err=20) reduceq
case('rgkmax')
  read(50,*,err=20) rgkmax
  if (rgkmax.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rgkmax <= 0 : ",G18.10)') rgkmax
    write(*,*)
    stop
  end if
case('gmaxvr')
  read(50,*,err=20) gmaxvr
case('lmaxapw')
  read(50,*,err=20) lmaxapw
  if (lmaxapw.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxapw < 0 : ",I8)') lmaxapw
    write(*,*)
    stop
  end if
  if (lmaxapw.ge.maxlapw) then
    write(*,*)
    write(*,'("Error(readinput): lmaxapw too large : ",I8)') lmaxapw
    write(*,'("Adjust maxlapw in modmain and recompile code")')
    write(*,*)
    stop
  end if
case('lmaxvr')
  read(50,*,err=20) lmaxvr
  if (lmaxvr.lt.3) then
    write(*,*)
    write(*,'("Error(readinput): lmaxvr < 3 : ",I8)') lmaxvr
    write(*,*)
    stop
  end if
case('lmaxmat')
  read(50,*,err=20) lmaxmat
  if (lmaxmat.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxmat < 0 : ",I8)') lmaxmat
    write(*,*)
    stop
  end if
case('lmaxinr')
  read(50,*,err=20) lmaxinr
  if (lmaxinr.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxinr < 0 : ",I8)') lmaxinr
    write(*,*)
    stop
  end if
case('fracinr')
  read(50,*,err=20) fracinr
case('npsden')
  read(50,*,err=20) npsden
  if (npsden.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): npsden < 2 : ",I8)') npsden
    write(*,*)
    stop
  end if
case('spinpol')
  read(50,*,err=20) spinpol
case('spinorb')
  read(50,*,err=20) spinorb
case('xctype')
  read(50,*,err=20) xctype
case('stype')
  read(50,*,err=20) stype
case('tarpack')
  read(50,*,err=20)tarpack
case('tlapack') 
  read(50,*,err=20)tlapack
case('tdiis')
  read(50,*,err=20)tdiis
case('tjdqz')
  read(50,*,err=20)tjdqz
case('diisfirstscl') 
  read(50,*,err=20)diisfirstscl 
!case('iterativetype')
 ! read(50,*) iterativetype
case('packedmatrixstorage')
  read(50,*,err=20) packedmatrixstorage
case('epsarpack')
  read(50,*,err=20) epsarpack
case('epsresid')
  read(50,*,err=20) epsresid
case ('maxncv')
  read(50,*,err=20) maxncv
case('lowesteval')
  read(50,*,err=20) lowesteval
case('swidth')
  read(50,*,err=20) swidth
  if (swidth.lt.1.d-9) then
    write(*,*)
    write(*,'("Error(readinput): swidth too small or negative : ",G18.10)') &
     swidth
    write(*,*)
    stop
  end if
case('epsocc')
  read(50,*,err=20) epsocc
  if (epsocc.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsocc <= 0 : ",G18.10)') epsocc
    write(*,*)
    stop
  end if
case('epschg')
  read(50,*,err=20) epschg
  if (epschg.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epschg <= 0 : ",G18.10)') epschg
    write(*,*)
    stop
  end if
case('nempty')
  read(50,*,err=20) nempty
  if (nempty.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nempty <= 0 : ",I8)') nempty
    write(*,*)
    stop
  end if
case('mixtype')
  read(50,*,err=20) mixtype
case('beta0')
  read(50,*,err=20) beta0
  if (beta0.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): beta0 < 0 : ",G18.10)') beta0
    write(*,*)
    stop
  end if
case('betainc')
  read(50,*,err=20) betainc
  if (betainc.lt.1.d0) then
    write(*,*)
    write(*,'("Error(readinput): betainc < 1 : ",G18.10)') betainc
    write(*,*)
    stop
  end if
case('betadec')
  read(50,*,err=20) betadec
  if ((betadec.le.0.d0).or.(betadec.gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): betadec should be in (0,1] : ",G18.10)') &
     betadec
    write(*,*)
    stop
  end if
case('maxscl')
  read(50,*,err=20) maxscl
  if (maxscl.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): maxscl < 0 : ",I8)') maxscl
    write(*,*)
    stop
  end if
case('epspot')
  read(50,*,err=20) epspot
case('epsengy')
  read(50,*,err=20) epsengy
case('epsforce')
  read(50,*,err=20) epsforce
case('cfdamp')
  read(50,*,err=20) cfdamp
  if (cfdamp.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): cfdamp < 0 : ",G18.10)') cfdamp
    write(*,*)
    stop
  end if
case('sppath')
  read(50,*,err=20) sppath
  sppath=adjustl(sppath)
case('scrpath')
  read(50,*,err=20) scrpath
case('molecule')
  read(50,*,err=20) molecule
case('vacuum')
  read(50,*,err=20) vacuum
  if (vacuum.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): vacuum < 0 : ",G18.10)') vacuum
    write(*,*)
    stop
  end if
case('atoms')
  read(50,*,err=20) nspecies
  if (nspecies.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nspecies <= 0 : ",I8)') nspecies
    write(*,*)
    stop
  end if
  if (nspecies.gt.maxspecies) then
    write(*,*)
    write(*,'("Error(readinput): nspecies too large : ",I8)') nspecies
    write(*,'("Adjust maxspecies in modmain and recompile code")')
    write(*,*)
    stop
  end if
  do is=1,nspecies
    read(50,*,err=20) spfname(is)
    spfname(is)=adjustl(spfname(is))
    read(50,*,err=20) natoms(is)
    if (natoms(is).le.0) then
      write(*,*)
      write(*,'("Error(readinput): natoms <= 0 : ",I8)') natoms(is)
      write(*,'(" for species ",I4)') is
      write(*,*)
      stop
    end if
    if (natoms(is).gt.maxatoms) then
      write(*,*)
      write(*,'("Error(readinput): natoms too large : ",I8)') natoms(is)
      write(*,'(" for species ",I4)') is
      write(*,'("Adjust maxatoms in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do ia=1,natoms(is)
      read(50,*,err=20) atposl(:,ia,is),bfcmt(:,ia,is)
    end do
  end do
case('plot1d')
  read(50,*,err=20) nvp1d,npp1d
  if (nvp1d.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): nvp1d < 1 : ",I8)') nvp1d
    write(*,*)
    stop
  end if
  if (npp1d.lt.nvp1d) then
    write(*,*)
    write(*,'("Error(readinput): npp1d < nvp1d : ",2I8)') npp1d,nvp1d
    write(*,*)
    stop
  end if
  if (allocated(vvlp1d)) deallocate(vvlp1d)
  allocate(vvlp1d(3,nvp1d))
  do iv=1,nvp1d
    read(50,*,err=20) vvlp1d(:,iv)
  end do
case('plot2d')
  read(50,*,err=20) vclp2d(:,1)
  read(50,*,err=20) vclp2d(:,2)
  read(50,*,err=20) vclp2d(:,3)
  read(50,*,err=20) np2d(:)
  if ((np2d(1).lt.1).or.(np2d(2).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): np2d < 1 : ",2I8)') np2d
    write(*,*)
    stop
  end if
case('plot3d')
  read(50,*,err=20) nup3d(:)
  if ((nup3d(1).lt.1).or.(nup3d(2).lt.1).or.(nup3d(3).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): nup3d < 1 : ",3I8)') nup3d
    write(*,*)
    stop
  end if
  read(50,*,err=20) np3d(:)
  if ((np3d(1).lt.1).or.(np3d(2).lt.1).or.(np3d(3).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): np3d < 1 : ",3I8)') np3d
    write(*,*)
    stop
  end if
case('dos')
  read(50,*,err=20) nwdos,ngrdos,nsmdos
  if (nwdos.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): nwdos < 2 : ",I8)') nwdos
    write(*,*)
    stop
  end if
  if (ngrdos.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): ngrdos < 1 : ",I8)') ngrdos
    write(*,*)
    stop
  end if
  if (nsmdos.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nsmdos < 0 : ",I8)') nsmdos
    write(*,*)
    stop
  end if
  read(50,*,err=20) wdos(:)
  if (wdos(1).ge.wdos(2)) then
    write(*,*)
    write(*,'("Error(readinput): wdos(1) >= wdos(2) : ",2G18.10)') wdos
    write(*,*)
    stop
  end if
case('lmirep')
  read(50,*,err=20) lmirep
case('tau0atm')
  read(50,*,err=20) tau0atm
case('nstfsp')
  read(50,*,err=20) nstfsp
  if (nstfsp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nstfsp <= 0 : ",I8)') nstfsp
    write(*,*)
    stop
  end if
case('lradstp')
  read(50,*,err=20) lradstp
  if (lradstp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): lradstp <= 0 : ",I8)') lradstp
    write(*,*)
    stop
  end if
case('chgexs')
  read(50,*,err=20) chgexs
case('nprad')
  read(50,*,err=20) nprad
  if (nprad.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): nprad < 2 : ",I8)') nprad
    write(*,*)
    stop
  end if
case('scissor')
  read(50,*,err=20) scissor
case('optcomp')
  do i=1,27
    read(50,'(A80)') str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty optical component list")')
        write(*,*)
        stop
      end if
      noptcomp=i-1
      goto 10
    end if
    str=trim(str)//' 1'
    read(str,*,iostat=iostat) optcomp(:,i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading optical component list")')
      write(*,'("(blank line required after optcomp block)")')
      write(*,*)
      stop
    end if
    if ((optcomp(1,i).lt.1).or.(optcomp(1,i).gt.3).or. &
        (optcomp(2,i).lt.1).or.(optcomp(2,i).gt.3).or. &
        (optcomp(3,i).lt.1).or.(optcomp(3,i).gt.3)) then
      write(*,*)
      write(*,'("Error(readinput): invalid optcomp : ",3I8)') optcomp
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): optical component list too long")')
  write(*,*)
  stop
!</sag>
case('optswidth')
  read(50,*,err=20) optswidth
  if (optswidth.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): optswidth < 0 : ",G18.10)') optswidth
    write(*,*)
    stop
  end if
!</sag>
case('usegdft')
  read(50,*,err=20) usegdft
case('intraband')
  read(50,*,err=20) intraband
case('evaltol')
  read(50,*,err=20) evaltol
  if (evaltol.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): evaltol <= 0 : ",G18.10)') evaltol
    write(*,*)
    stop
  end if
case('evalmin')
  read(50,*,err=20) evalmin
case('deband')
  read(50,*,err=20) deband
  if (deband.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): deband < 0 : ",G18.10)') deband
    write(*,*)
    stop
  end if
case('bfieldc')
  read(50,*,err=20) bfieldc
case('fixspin')
  read(50,*,err=20) fixspin
case('momfix')
  read(50,*,err=20) momfix
case('mommtfix')
  do ias=1,maxspecies*maxatoms
    read(50,'(A80)') str
    if (trim(str).eq.'') goto 10
    read(str,*,iostat=iostat) is,ia,mommtfix(:,ia,is)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading muffin-tin fixed spin &
       &moments")')
      write(*,'("(blank line required after mommtfix block")')
      write(*,*)
      stop
    end if
  end do
case('taufsm')
  read(50,*,err=20) taufsm
  if (taufsm.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taufsm < 0 : ",G18.10)') taufsm
    write(*,*)
    stop
  end if
case('autormt')
  read(50,*,err=20) autormt
case('rmtapm')
  read(50,*,err=20) rmtapm(:)
  if (rmtapm(1).lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rmtapm(1) < 0 : ",G18.10)') rmtapm(1)
    write(*,*)
    stop
  end if
  if ((rmtapm(2).le.0.d0).or.(rmtapm(2).gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): rmtapm(2) not in (0,1] : ",G18.10)') rmtapm(2)
    write(*,*)
    stop
  end if
case('isgkmax')
  read(50,*,err=20) isgkmax
case('nosym')
  read(50,*,err=20) nosym
case('deltaph')
  read(50,*,err=20) deltaph
case('phwrite')
  read(50,*,err=20) nphwrt
  if (nphwrt.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nphwrt <= 0 : ",I8)') nphwrt
    write(*,*)
    stop
  end if
  if (allocated(vqlwrt)) deallocate(vqlwrt)
  allocate(vqlwrt(3,nphwrt))
  do i=1,nphwrt
    read(50,*,err=20) vqlwrt(:,i)
  end do
case('notes')
  do i=1,maxnlns
    read(50,'(A80)') notes(i)
    if (trim(notes(i)).eq.'') then
      notelns=i-1
      goto 10
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): too many note lines")')
  write(*,*)
  stop
case('tforce')
  read(50,*,err=20) tforce
case('tfibs')
  read(50,*,err=20) tfibs
case('maxitoep')
  read(50,*,err=20) maxitoep
  if (maxitoep.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): maxitoep < 1 : ",I8)') maxitoep
    write(*,*)
    stop
  end if
case('tauoep')
  read(50,*,err=20) tauoep(:)
  if ((tauoep(1).lt.0.d0).or.(tauoep(2).lt.0.d0).or.(tauoep(3).lt.0.d0)) then
    write(*,*)
    write(*,'("Error(readinput): tauoep < 0 : ",3G18.10)') tauoep
    write(*,*)
    stop
  end if
case('kstlist')
  do i=1,maxkst
    read(50,'(A80)') str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty k-point and state list")')
        write(*,*)
        stop
      end if
      nkstlist=i-1
      goto 10
    end if
    str=trim(str)//' 1 1'
    read(str,*,iostat=iostat) kstlist(:,i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading k-point and state list")')
      write(*,'("(blank line required after kstlist block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): k-point and state list too long")')
  write(*,*)
  stop
case('vklem')
  read(50,*,err=20) vklem
case('deltaem')
  read(50,*,err=20) deltaem
case('ndspem')
  read(50,*,err=20) ndspem
  if ((ndspem.lt.1).or.(ndspem.gt.3)) then
    write(*,*)
    write(*,'("Error(readinput): ndspem out of range : ",I8)') ndspem
    write(*,*)
    stop
  end if
case('nosource')
  read(50,*,err=20) nosource
case('spinsprl')
  read(50,*,err=20) spinsprl
case('vqlss')
  read(50,*,err=20) vqlss
case('nwrite')
  read(50,*,err=20) nwrite
case('tevecsv')
  read(50,*,err=20) tevecsv
case('lda+u')
  read(50,*) ldapu
  do is=1,maxspecies
    read(50,'(A80)') str
    if (trim(str).eq.'') goto 10
    read(str,*,iostat=iostat) js,l,t1,t2
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading LDA+U parameters")')
      write(*,'("(blank line required after lda+u block)")')
      write(*,*)
      stop
    end if
    if ((js.le.0).or.(js.ge.maxspecies)) then
      write(*,*)
      write(*,'("Error(readinput): invalid species number in lda+u block : ",&
       &I8)') js
      write(*,*)
      stop
    end if
    if (l.gt.lmaxlu) then
      write(*,*)
      write(*,'("Error(readinput): l > lmaxlu in lda+u block : ",2I8)') l,lmaxlu
      write(*,*)
      stop
    end if
    llu(js)=l
    ujlu(1,js)=t1
    ujlu(2,js)=t2
  end do
case('rdmxctype')
  read(50,*,err=20) rdmxctype
case('rdmmaxscl')
  read(50,*,err=20) rdmmaxscl
  if (rdmmaxscl.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): rdmmaxscl < 0 : ",I8)') rdmmaxscl
    write(*,*)
  end if
case('maxitn')
  read(50,*,err=20) maxitn
case('maxitc')
  read(50,*,err=20) maxitc
case('taurdmn')
  read(50,*,err=20) taurdmn
  if (taurdmn.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taurdmn < 0 : ",G18.10)') taurdmn
    write(*,*)
    stop
  end if
case('taurdmc')
  read(50,*,err=20) taurdmc
  if (taurdmc.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taurdmc < 0 : ",G18.10)') taurdmc
    write(*,*)
    stop
  end if
case('rdmalpha')
  read(50,*,err=20) rdmalpha
  if ((rdmalpha.le.0.d0).or.(rdmalpha.ge.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): rdmalpha not in (0,1) : ",G18.10)') rdmalpha
    write(*,*)
    stop
  end if
case('rdmtemp')
  read(50,*,err=20) rdmtemp
  if (rdmtemp.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rdmtemp < 0 : ",G18.10)') rdmtemp
    write(*,*)
    stop
  end if
case('reducebf')
  read(50,*,err=20) reducebf
  if ((reducebf.lt.0.d0).or.(reducebf.gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): reducebf not in [0,1] : ",G18.10)') reducebf
    write(*,*)
    stop
  end if
case('ptnucl')
  read(50,*,err=20) ptnucl
case('tseqit')
  read(50,*,err=20) tseqit
case('nseqit')
  read(50,*,err=20) nseqit
  if (nseqit.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): nseqit < 1 : ",I8)') nseqit
    write(*,*)
    stop
  end if
case('vecql')
  read(50,*,err=20) vecql(:)
case('mustar')
  read(50,*,err=20) mustar
case('sqados')
  read(50,*,err=20) sqados(:)
#ifdef TETRA
!  tetrahedron method variables
case('tetraocc')
  read(50,*,err=20) tetraocc
case('tetraopt')
  read(50,*,err=20) tetraopt
case('tetradf')
  read(50,*,err=20) tetradf
case('tetrakordexc')
  read(50,*,err=20) tetrakordexc
#endif
#ifdef XS
! TDDFT variables
case('imbandstr')
  read(50,*,err=20) imbandstr
case('vgqlmt')
  read(50,*,err=20) nqptmt
  if (nqptmt.le.0) then
    write(*,*)
    write(*,'("Error(readinput/xs): nqptmt <= 0 : ",I8)') nqptmt
    write(*,*)
    stop
  end if
  if (allocated(vgqlmt)) deallocate(vgqlmt)
  allocate(vgqlmt(3,nqptmt))
  do i=1,nqptmt
    read(50,*,err=20) vgqlmt(:,i)
  end do
case('mdfqtype')
  read(50,*,err=20) mdfqtype
  if ((mdfqtype.lt.0).or.(mdfqtype.gt.1)) then
    write(*,*)
    write(*,'("Error(readinput/xs): mdfqtype not in {0,1} : ",I8)') mdfqtype
    write(*,*)
    stop
  end if
case('tq0ev')
   read(50,*,err=20) tq0ev
case('gqmax')
  read(50,*,err=20) gqmax
case('lmaxapwwf')
  read(50,*,err=20) lmaxapwwf
  if (lmaxapwwf.lt.0) then
    write(*,*)
    write(*,'("Error(readinput)/xs: lmaxapwwf < 0 : ",I8)') lmaxapwwf
    write(*,*)
    stop
  end if
case('fastpmat')
  read(50,*,err=20) fastpmat
case('fastemat')
  read(50,*,err=20) fastemat
case('emattype')
  read(50,*,err=20) emattype
case('lmaxemat')
  read(50,*,err=20) lmaxemat
  if (lmaxemat.lt.0) then
    write(*,*)
    write(*,'("Error(readinput/xs): lmaxemat < 0 : ",I8)') lmaxemat
    write(*,*)
    stop
  end if
case('rsptype')
  read(50,*,err=20) rsptype
  if ((rsptype.ne.'tord').and.(rsptype.ne.'reta')) then
    write(*,*)
    write(*,'("Error(readinput/xs): invalid rsptype : ",a)') rsptype
    write(*,*)
    stop
  end if
case('acont')
  read(50,*,err=20) acont
case('nwacont')
  read(50,*,err=20) nwacont
  if (broad.le.0) then
    write(*,*)
    write(*,'("Error(readinput/xs): nwacont <= 0 : ",g18.10)') nwacont
    write(*,*)
    stop
  end if
case('broad')
  read(50,*,err=20) broad
  if (broad.le.0) then
    write(*,*)
    write(*,'("Warning(readinput/xs): broad <= 0 : ",g18.10)') broad
    write(*,*)
  end if
case('aresdf')
  read(50,*,err=20) aresdf
case('epsdfde')
  read(50,*,err=20) epsdfde
  if (broad.le.0) then
    write(*,*)
    write(*,'("Warning(readinput/xs): epsdfde <= 0 : ",g18.10)') epsdfde
    write(*,*)
  end if
case('emaxdf')
  read(50,*,err=20) emaxdf
case('dfoffdiag')
  read(50,*,err=20) dfoffdiag
case('kerndiag')
  read(50,*,err=20) kerndiag
case('fxctype')
  read(50,*,err=20) fxctype
case('nexcitmax')
  read(50,*,err=20) nexcitmax
case('alphalrc')
  read(50,*,err=20) alphalrc
case('alphalrcdyn')
  read(50,*,err=20) alphalrcdyn
case('betalrcdyn')
  read(50,*,err=20) betalrcdyn
case('dftrans')
  read(50,*,err=20) ndftrans
  if (allocated(dftrans)) deallocate(dftrans)
  allocate(dftrans(3,ndftrans))
  do i=1,ndftrans
    read(50,'(A80)') str
    if (trim(str).eq.'') then
      write(*,*)
      write(*,'("Error(readinput/xs): missing k-point and state in list for &
           &transition analysis")')
      write(*,*)
      stop
    end if
    read(str,*,iostat=iostat) dftrans(:,i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput/xs): error reading k-point and state list for&
           &transition analysis")')
      write(*,'("(blank line required after dftrans block)")')
      write(*,*)
      stop
    end if
  end do
case('tetraqweights')
  read(50,*,err=20) tetraqweights
case('gather')
  read(50,*,err=20) gather
case('symmorph')
  read(50,*,err=20) symmorph
case('tevout')
  read(50,*,err=20) tevout
case('appinfo')
  read(50,*,err=20) tappinfo
case('dbglev')
  read(50,*,err=20) dbglev
  ! screening variables
case('screentype')
   read(50,*,err=20) screentype
   select case(trim(screentype))
   case('full','diag','noinvdiag','longrange')
   case default
      write(*,*)
      write(*,'("Error(readinput/xs): unknown screening type: ",a)') &
           trim(screentype)
      write(*,*)
      stop
   end select
case('nosymscr')
  read(50,*,err=20) nosymscr
case('reducekscr')
  read(50,*,err=20) reducekscr
case('ngridkscr')
  read(50,*,err=20) ngridkscr(1),ngridkscr(2),ngridkscr(3)
  if ((ngridkscr(1).le.0).or.(ngridkscr(2).le.0).or.(ngridkscr(3).le.0)) then
    write(*,*)
    write(*,'("Error(readinput/xs): invalid ngridkscr : ",3I8)') ngridkscr
    write(*,*)
    stop
  end if
case('vkloffscr')
  read(50,*,err=20) vkloffscr(1),vkloffscr(2),vkloffscr(3)
case('rgkmaxscr')
  read(50,*,err=20) rgkmaxscr
  if (rgkmaxscr.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput/xs): rgkmaxscr <= 0 : ",G18.10)') rgkmaxscr
    write(*,*)
    stop
  end if
case('nemptyscr')
  read(50,*,err=20) nemptyscr
  if (nemptyscr.le.0) then
    write(*,*)
    write(*,'("Error(readinput/xs): nemptyscr <= 0 : ",I8)') nemptyscr
    write(*,*)
    stop
  end if
case('scrherm')
  read(50,*,err=20) scrherm
  ! BSE (-kernel) variables
case('bsetype')
   read(50,*,err=20) bsetype
   select case(trim(bsetype))
   case('ip','rpa','singlet','triplet')
   case default
      write(*,*)
      write(*,'("Error(readinput/xs): unknown BSE-type: ",a)') &
           trim(bsetype)
      write(*,*)
      stop
   end select
case('nosymbse')
  read(50,*,err=20) nosymbse
case('reducekbse')
  read(50,*,err=20) reducekbse
case('vkloffbse')
  read(50,*,err=20) vkloffbse(1),vkloffbse(2),vkloffbse(3)
case('bsediagweight')
  read(50,*,err=20) bsediagweight
case('bsediagsym')
  read(50,*,err=20) bsediagsym
case('nstlce')
  read(50,*,err=20) nbfce,nafce
  if ((nbfce.le.0).or.(nafce.le.0)) then
    write(*,*)
    write(*,'("Error(readinput/xs): nbfce or nafce <= 0 : ",I8)') nbfce,nafce
    write(*,*)
    stop
  end if
case('nstlbse')
  read(50,*,err=20) nbfbse,nafbse
  if ((nbfbse.le.0).or.(nafbse.le.0)) then
    write(*,*)
    write(*,'("Error(readinput/xs): nbfbse or nafbse <= 0 : ",I8)') nbfbse, &
         nafbse
    write(*,*)
    stop
  end if
#endif
case('')
  goto 10
case default
  write(*,*)
  write(*,'("Error(readinput): invalid block name : ",A)') trim(bname)
  write(*,*)
  stop
end select
goto 10
20 continue
write(*,*)
write(*,'("Error(readinput): error reading from exciting.in")')
write(*,'("Problem occurred in ''",A,"'' block")') trim(bname)
write(*,*)
stop
30 continue
close(50)
! scale the lattice vectors (scaling not referenced again in code)
avec(:,1)=sc1*avec(:,1)
avec(:,2)=sc2*avec(:,2)
avec(:,3)=sc3*avec(:,3)
avec(:,:)=sc*avec(:,:)
! check if system is an isolated molecule
if (molecule) then
! set up cubic unit cell with vacuum region around molecule
  avec(:,:)=0.d0
  do is=1,nspecies
    do ia=1,natoms(is)
      do js=1,nspecies
        do ja=1,natoms(is)
          do i=1,3
            t1=abs(atposl(i,ia,is)-atposl(i,ja,js))
            if (t1.gt.avec(i,i)) avec(i,i)=t1
          end do
        end do
      end do
    end do
  end do
  do i=1,3
    avec(i,i)=avec(i,i)+vacuum
  end do
! convert atomic positions from Cartesian to lattice coordinates
  call r3minv(avec,ainv)
  do is=1,nspecies
    do ia=1,natoms(is)
      call r3mv(ainv,atposl(:,ia,is),v)
      atposl(:,ia,is)=v(:)
    end do
  end do
end if
#ifdef XS
! dump default parameters (partially) redifined in input file
if (rank.eq.0) then
   fname='PARAMS.OUT'
   if (dumpmain) call dumpparams(trim(fname), &
        '! main parameters:',sppath,sc,sc1,sc2,sc3,vacuum)
   if (dumpadd) call dumpparams_add(trim(fname), &
        '! additional parameters:')
   if (dumpmpiiter) call dumpparams_mpiiter(trim(fname), &
        '! MPI parallelization and iterative solver parameters:')
   if (dumptetra) call dumpparams_tetra(trim(fname), &
        '! tetrahedron method parameters:')
   if (dumpxs) call dumpparams_xs(trim(fname), &
        '! excited states parameters:')
end if
#endif

!---------------------------------------------!
!     read from atomic species data files     !
!---------------------------------------------!
do is=1,nspecies
  open(50,file=trim(sppath)//trim(spfname(is)),action='READ',status='OLD', &
   form='FORMATTED',iostat=iostat)
  if (iostat.ne.0) then
    write(*,*)
    write(*,'("Error(readinput): error opening species file ",A)') &
     trim(sppath)//trim(spfname(is))
    write(*,*)
    stop
  end if
  read(50,*) spsymb(is)
  read(50,*) spname(is)
  read(50,*) spzn(is)
  read(50,*) spmass(is)
  read(50,*) sprmin(is),rmt(is),sprmax(is),nrmt(is)
  if (sprmin(is).le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): sprmin <= 0 : ",G18.10)') sprmin(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (rmt(is).le.sprmin(is)) then
    write(*,*)
    write(*,'("Error(readinput): rmt <= sprmin : ",2G18.10)') rmt(is),sprmin(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (sprmax(is).lt.rmt(is)) then
    write(*,*)
    write(*,'("Error(readinput): sprmax < rmt : ",2G18.10)') sprmax(is),rmt(is)
    write(*,*)
    stop
  end if
  if (nrmt(is).lt.20) then
    write(*,*)
    write(*,'("Error(readinput): nrmt too small : ",I8)') nrmt(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  read(50,*) spnst(is)
  if (spnst(is).le.0) then
    write(*,*)
    write(*,'("Error(readinput): invalid spnst : ",I8)') spnst(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (spnst(is).gt.maxspst) then
    write(*,*)
    write(*,'("Error(readinput): too many states for species ",I8)') is
    write(*,*)
    stop
  end if
  do ist=1,spnst(is)
    read(50,*) spn(ist,is),spl(ist,is),spk(ist,is),spocc(ist,is),spcore(ist,is)
    if (spn(ist,is).lt.1) then
      write(*,*)
      write(*,'("Error(readinput): spn < 1 : ",I8)') spn(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (spl(ist,is).lt.0) then
      write(*,*)
      write(*,'("Error(readinput): spl < 0 : ",I8)') spl(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (spk(ist,is).lt.1) then
      write(*,*)
      write(*,'("Error(readinput): spk < 1 : ",I8)') spk(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (spocc(ist,is).lt.0.d0) then
      write(*,*)
      write(*,'("Error(readinput): spocc < 0 : ",G18.10)') spocc(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
  end do
  read(50,*) apword(0,is)
  if (apword(0,is).le.0) then
    write(*,*)
    write(*,'("Error(readinput): apword <= 0 : ",I8)') apword(0,is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (apword(0,is).gt.maxapword) then
    write(*,*)
    write(*,'("Error(readinput): apword too large : ",I8)') apword(0,is)
    write(*,'(" for species ",I4)') is
    write(*,'("Adjust maxapword in modmain and recompile code")')
    write(*,*)
    stop
  end if
! set the APW orders for l>0
  apword(1:lmaxapw,is)=apword(0,is)
  do io=1,apword(0,is)
    read(50,*) apwe0(io,0,is),apwdm(io,0,is),apwve(io,0,is)
    if (apwdm(io,0,is).lt.0) then
      write(*,*)
      write(*,'("Error(readinput): apwdm < 0 : ",I8)') apwdm(io,0,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and order ",I4)') io
      write(*,*)
      stop
    end if
! set the APW linearisation energies, derivative orders and variability for l>0
    apwe0(io,1:lmaxapw,is)=apwe0(io,0,is)
    apwdm(io,1:lmaxapw,is)=apwdm(io,0,is)
    apwve(io,1:lmaxapw,is)=apwve(io,0,is)
  end do
  read(50,*) nlx
  if (nlx.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nlx < 0 : ",I8)') nlx
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  do ilx=1,nlx
    read(50,*) lx,io
    if (lx.lt.0) then
      write(*,*)
      write(*,'("Error(readinput): lx < 0 : ",I8)') lx
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    if (lx.gt.lmaxapw) then
      write(*,*)
      write(*,'("Error(readinput): lx > lmaxapw : ",I8)') lx
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    apword(lx,is)=io
    if (apword(lx,is).le.0) then
      write(*,*)
      write(*,'("Error(readinput): apword <= 0 : ",I8)') apword(lx,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    if (apword(lx,is).gt.maxapword) then
      write(*,*)
      write(*,'("Error(readinput): apword too large : ",I8)') apword(lx,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,'("Adjust maxapword in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do io=1,apword(lx,is)
      read(50,*) apwe0(io,lx,is),apwdm(io,lx,is),apwve(io,lx,is)
      if (apwdm(io,lx,is).lt.0) then
        write(*,*)
        write(*,'("Error(readinput): apwdm < 0 : ",I8)') apwdm(io,lx,is)
        write(*,'(" for species ",I4)') is
        write(*,'(" exception number ",I4)') ilx
        write(*,'(" and order ",I4)') io
        write(*,*)
        stop
      end if
    end do
  end do
  read(50,*) nlorb(is)
  if (nlorb(is).lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nlorb < 0 : ",I8)') nlorb(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (nlorb(is).gt.maxlorb) then
    write(*,*)
    write(*,'("Error(readinput): nlorb too large : ",I8)') nlorb(is)
    write(*,'(" for species ",I4)') is
    write(*,'("Adjust maxlorb in modmain and recompile code")')
    write(*,*)
    stop
  end if
  do ilo=1,nlorb(is)
    read(50,*) lorbl(ilo,is),lorbord(ilo,is)
    if (lorbl(ilo,is).lt.0) then
      write(*,*)
      write(*,'("Error(readinput): lorbl < 0 : ",I8)') lorbl(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbl(ilo,is).gt.lmaxmat) then
      write(*,*)
      write(*,'("Error(readinput): lorbl > lmaxmat : ",2I8)') lorbl(ilo,is), &
       lmaxmat
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbord(ilo,is).lt.2) then
      write(*,*)
      write(*,'("Error(readinput): lorbord < 2 : ",I8)') lorbord(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbord(ilo,is).gt.maxlorbord) then
      write(*,*)
      write(*,'("Error(readinput): lorbord too large : ",I8)') lorbord(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,'("Adjust maxlorbord in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do io=1,lorbord(ilo,is)
      read(50,*) lorbe0(io,ilo,is),lorbdm(io,ilo,is),lorbve(io,ilo,is)
      if (lorbdm(io,ilo,is).lt.0) then
        write(*,*)
        write(*,'("Error(readinput): lorbdm < 0 : ",I8)') lorbdm(io,ilo,is)
        write(*,'(" for species ",I4)') is
        write(*,'(" local-orbital ",I4)') ilo
        write(*,'(" and order ",I4)') io
        write(*,*)
        stop
      end if
    end do
  end do
  close(50)
end do
end subroutine
!EOC
