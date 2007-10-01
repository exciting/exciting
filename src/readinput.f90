
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readinput
! !INTERFACE:
subroutine readinput
! !USES:
use modmain
! !DESCRIPTION:
!   Reads in the input parameters from the file {\tt exciting.in} as well as
!   from the species files. Also sets default values for the input parameters.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ia,is,ja,js,iv,i,iostat
integer ist,io,nlx,ilx,lx,ilo
real(8) sc,sc1,sc2,sc3
real(8) vacuum,t1
character(256) str
character(256) bname
character(256) sppath

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
rlambda=30.d0
autokpt=.false.
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
beta0=0.1d0
betamax=1.d0
maxscl=200
epspot=1.d-6
epsengy=1.d-7
epsforce=5.d-4
cfdamp=0.d0
molecule=.false.
vacuum=10.d0
nspecies=0
natoms(:)=0
sppath='./'
scrpath='./'
nvp1d=2
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
bcsym=.true.
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
usegdft=.false.
intraband=.false.
evalmin=-4.5d0
deband=0.0025d0
bfieldc(:)=0.d0
fixspin=.false.
momfix(:)=0.d0
taufsm=0.01d0
autormt=.false.
rmtapm(1)=0.25d0
rmtapm(2)=0.95d0
nosym=.false.
deltaph=0.03d0
nphwrt=1
if (allocated(vqlwrt)) deallocate(vqlwrt)
allocate(vqlwrt(3,nphwrt))
vqlwrt(:,:)=0.d0
notelns=0
tforce=.false.
tfibs=.true.
maxitoep=150
tau0oep=0.5d0
dtauoep=0.5d0
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

!-------------------------------!
!     read from exciting.in     !
!-------------------------------!
open(50,file='exciting.in',action='READ',status='OLD',form='FORMATTED')
10 continue
read(50,*,end=20) bname
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
      write(*,'(" (blank line required after tasks block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): too many tasks")')
  write(*,*)
  stop
case('avec')
  read(50,*) avec(1,1),avec(2,1),avec(3,1)
  read(50,*) avec(1,2),avec(2,2),avec(3,2)
  read(50,*) avec(1,3),avec(2,3),avec(3,3)
case('scale')
  read(50,*) sc
case('scale1')
  read(50,*) sc1
case('scale2')
  read(50,*) sc2
case('scale3')
  read(50,*) sc3
case('epslat')
  read(50,*) epslat
  if (epslat.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epslat <= 0 : ",G18.10)') epslat
    write(*,*)
    stop
  end if
case('primcell')
  read(50,*) primcell
case('tshift')
  read(50,*) tshift
case('rlambda')
  read(50,*) rlambda
  if (rlambda.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rlambda <= 0 : ",G18.10)') rlambda
    write(*,*)
    stop
  end if
case('autokpt')
  read(50,*) autokpt
case('ngridk')
  read(50,*) ngridk(1),ngridk(2),ngridk(3)
  if ((ngridk(1).le.0).or.(ngridk(2).le.0).or.(ngridk(3).le.0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid ngridk : ",3I8)') ngridk
    write(*,*)
    stop
  end if
case('vkloff')
  read(50,*) vkloff(1),vkloff(2),vkloff(3)
case('reducek')
  read(50,*) reducek
case('ngridq')
  read(50,*) ngridq(1),ngridq(2),ngridq(3)
  if ((ngridq(1).le.0).or.(ngridq(2).le.0).or.(ngridq(3).le.0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid ngridq : ",3I8)') ngridq
    write(*,*)
    stop
  end if
case('reduceq')
  read(50,*) reduceq
case('rgkmax')
  read(50,*) rgkmax
  if (rgkmax.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rgkmax <= 0 : ",G18.10)') rgkmax
    write(*,*)
    stop
  end if
case('gmaxvr')
  read(50,*) gmaxvr
case('lmaxapw')
  read(50,*) lmaxapw
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
  read(50,*) lmaxvr
  if (lmaxvr.lt.3) then
    write(*,*)
    write(*,'("Error(readinput): lmaxvr < 3 : ",I8)') lmaxvr
    write(*,*)
    stop
  end if
case('lmaxmat')
  read(50,*) lmaxmat
  if (lmaxmat.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxmat < 0 : ",I8)') lmaxmat
    write(*,*)
    stop
  end if
case('lmaxinr')
  read(50,*) lmaxinr
  if (lmaxinr.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxinr < 0 : ",I8)') lmaxinr
    write(*,*)
    stop
  end if
case('fracinr')
  read(50,*) fracinr
case('npsden')
  read(50,*) npsden
  if (npsden.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): npsden < 2 : ",I8)') npsden
    write(*,*)
    stop
  end if
case('spinpol')
  read(50,*) spinpol
case('spinorb')
  read(50,*) spinorb
case('xctype')
  read(50,*) xctype
case('stype')
  read(50,*) stype
case('swidth')
  read(50,*) swidth
  if (swidth.lt.1.d-9) then
    write(*,*)
    write(*,'("Error(readinput): swidth too small or negative : ",G18.10)') &
     swidth
    write(*,*)
    stop
  end if
case('epsocc')
  read(50,*) epsocc
  if (epsocc.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsocc <= 0 : ",G18.10)') epsocc
    write(*,*)
    stop
  end if
case('epschg')
  read(50,*) epschg
  if (epschg.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epschg <= 0 : ",G18.10)') epschg
    write(*,*)
    stop
  end if
case('nempty')
  read(50,*) nempty
  if (nempty.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nempty <= 0 : ",I8)') nempty
    write(*,*)
    stop
  end if
case('beta0')
  read(50,*) beta0
  if (beta0.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): beta0 < 0 : ",G18.10)') beta0
    write(*,*)
    stop
  end if
case('betamax')
  read(50,*) betamax
  if (betamax.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): betamax < 0 : ",G18.10)') betamax
    write(*,*)
    stop
  end if
case('maxscl')
  read(50,*) maxscl
  if (maxscl.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): maxscl < 0 : ",I8)') maxscl
    write(*,*)
    stop
  end if
case('epspot')
  read(50,*) epspot
case('epsengy')
  read(50,*) epsengy
case('epsforce')
  read(50,*) epsforce
case('cfdamp')
  read(50,*) cfdamp
  if (cfdamp.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): cfdamp < 0 : ",G18.10)') cfdamp
    write(*,*)
    stop
  end if
case('sppath')
  read(50,*) sppath
case('scrpath')
  read(50,*) scrpath
case('molecule')
  read(50,*) molecule
case('vacuum')
  read(50,*) vacuum
  if (vacuum.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): vacuum < 0 : ",G18.10)') vacuum
    write(*,*)
    stop
  end if
case('atoms')
  read(50,*) nspecies
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
    read(50,*) spfname(is)
    read(50,*) natoms(is)
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
      read(50,*) atposl(1,ia,is),atposl(2,ia,is),atposl(3,ia,is), &
       bfcmt(1,ia,is),bfcmt(2,ia,is),bfcmt(3,ia,is)
    end do
  end do
case('plot1d')
  read(50,*) nvp1d,npp1d
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
    read(50,*) vvlp1d(1,iv),vvlp1d(2,iv),vvlp1d(3,iv)
  end do
case('plot2d')
  read(50,*) vclp2d(1,1),vclp2d(2,1),vclp2d(3,1)
  read(50,*) vclp2d(1,2),vclp2d(2,2),vclp2d(3,2)
  read(50,*) vclp2d(1,3),vclp2d(2,3),vclp2d(3,3)
  read(50,*) np2d(1),np2d(2)
  if ((np2d(1).lt.1).or.(np2d(2).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): np2d < 1 : ",2I8)') np2d
    write(*,*)
    stop
  end if
case('plot3d')
  read(50,*) nup3d(1),nup3d(2),nup3d(3)
  if ((nup3d(1).lt.1).or.(nup3d(2).lt.1).or.(nup3d(3).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): nup3d < 1 : ",3I8)') nup3d
    write(*,*)
    stop
  end if
  read(50,*) np3d(1),np3d(2),np3d(3)
  if ((np3d(1).lt.1).or.(np3d(2).lt.1).or.(np3d(3).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): np3d < 1 : ",3I8)') np3d
    write(*,*)
    stop
  end if
case('dos')
  read(50,*) nwdos,ngrdos,nsmdos
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
  read(50,*) wdos(1),wdos(2)
  if (wdos(1).ge.wdos(2)) then
    write(*,*)
    write(*,'("Error(readinput): wdos(1) >= wdos(2) : ",2G18.10)') wdos
    write(*,*)
    stop
  end if
case('bcsym')
  read(50,*) bcsym
case('tau0atm')
  read(50,*) tau0atm
case('nstfsp')
  read(50,*) nstfsp
  if (nstfsp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nstfsp <= 0 : ",I8)') nstfsp
    write(*,*)
    stop
  end if
case('lradstp')
  read(50,*) lradstp
  if (lradstp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): lradstp <= 0 : ",I8)') lradstp
    write(*,*)
    stop
  end if
case('chgexs')
  read(50,*) chgexs
case('nprad')
  read(50,*) nprad
  if (nprad.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): nprad < 2 : ",I8)') nprad
    write(*,*)
    stop
  end if
case('scissor')
  read(50,*) scissor
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
    read(str,*,iostat=iostat) optcomp(:,i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading optical component list")')
      write(*,'(" (blank line required after optcomp block)")')
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
case('usegdft')
  read(50,*) usegdft
case('intraband')
  read(50,*) intraband
case('evalmin')
  read(50,*) evalmin
case('deband')
  read(50,*) deband
  if (deband.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): deband < 0 : ",G18.10)') deband
    write(*,*)
    stop
  end if
case('bfieldc')
  read(50,*) bfieldc
case('fixspin')
  read(50,*) fixspin
case('momfix')
  read(50,*) momfix
case('taufsm')
  read(50,*) taufsm
  if (taufsm.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taufsm < 0 : ",G18.10)') taufsm
    write(*,*)
    stop
  end if
case('autormt')
  read(50,*) autormt
case('rmtapm')
  read(50,*) rmtapm
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
case('nosym')
  read(50,*) nosym
case('deltaph')
  read(50,*) deltaph
case('phwrite')
  read(50,*) nphwrt
  if (nphwrt.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nphwrt <= 0 : ",I8)') nphwrt
    write(*,*)
    stop
  end if
  if (allocated(vqlwrt)) deallocate(vqlwrt)
  allocate(vqlwrt(3,nphwrt))
  do i=1,nphwrt
    read(50,*) vqlwrt(:,i)
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
  read(50,*) tforce
case('tfibs')
  read(50,*) tfibs
case('maxitoep')
  read(50,*) maxitoep
  if (maxitoep.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): maxitoep < 1 : ",I8)') maxitoep
    write(*,*)
    stop
  end if
case('tau0oep')
  read(50,*) tau0oep
case('dtauoep')
  read(50,*) dtauoep
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
    read(str,*,iostat=iostat) kstlist(:,i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading k-point and state list")')
      write(*,'(" (blank line required after kstlist block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): k-point and state list too long")')
  write(*,*)
  stop
case('vklem')
  read(50,*) vklem
case('deltaem')
  read(50,*) deltaem
case('ndspem')
  read(50,*) ndspem
  if ((ndspem.lt.1).or.(ndspem.gt.3)) then
    write(*,*)
    write(*,'("Error(readinput): ndspem out of range : ",I8)') ndspem
    write(*,*)
    stop
  end if
case('nosource')
  read(50,*) nosource
case('spinsprl')
  read(50,*) spinsprl
case('vqlss')
  read(50,*) vqlss
case('nwrite')
  read(50,*) nwrite
case('tevecsv')
  read(50,*) tevecsv
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
      call r3mv(ainv,atposl(1,ia,is),atposl(1,ia,is))
    end do
  end do
end if

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
