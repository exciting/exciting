
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: init0
! !INTERFACE:
subroutine init0
! !USES:
use modmain
use modxcifc
! !DESCRIPTION:
!   Performs basic consistency checks as well as allocating and initialising
!   global variables not dependent on the $k$-point set.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ist
integer l,m,lm,iv(3)
real(8) cs,sn,r,t1
real(8) cpu0,cpu1
! external functions
real(8) dlamch
external dlamch

!-------------------------------!
!     zero timing variables     !
!-------------------------------!
timeinit=0.d0
timemat=0.d0
timefv=0.d0
timesv=0.d0
timerho=0.d0
timepot=0.d0
timefor=0.d0
call cpu_time(cpu0)

!------------------------------------!
!     angular momentum variables     !
!------------------------------------!
lmmaxvr=(lmaxvr+1)**2
lmmaxapw=(lmaxapw+1)**2
lmmaxmat=(lmaxmat+1)**2
lmmaxinr=(lmaxinr+1)**2
if (lmaxvr.gt.lmaxapw) then
  write(*,*)
  write(*,'("Error(init0): lmaxvr > lmaxapw : ",2I8)') lmaxvr,lmaxapw
  write(*,*)
  stop
end if
if (lmaxmat.gt.lmaxapw) then
  write(*,*)
  write(*,'("Error(init0): lmaxmat > lmaxapw : ",2I8)') lmaxmat,lmaxapw
  write(*,*)
  stop
end if
! index to (l,m) pairs
if (allocated(idxlm)) deallocate(idxlm)
allocate(idxlm(0:lmaxapw,-lmaxapw:lmaxapw))
lm=0
do l=0,lmaxapw
  do m=-l,l
    lm=lm+1
    idxlm(l,m)=lm
  end do
end do
! array of i**l values
if (allocated(zil)) deallocate(zil)
allocate(zil(0:lmaxapw))
do l=0,lmaxapw
  zil(l)=zi**l
end do

!------------------------------------!
!     index to atoms and species     !
!------------------------------------!
! check if the system is an isolated molecule
if (molecule) then
  primcell=.false.
  tshift=.false.
end if
! find primitive cell if required
if (primcell) call findprim
natmmax=0
ias=0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=ias+1
    idxas(ia,is)=ias
  end do
! maximum number of atoms over all species
  natmmax=max(natmmax,natoms(is))
end do
! total number of atoms
natmtot=ias

!------------------------!
!     spin variables     !
!------------------------!
if (spinsprl) then
  select case(task)
  case(2,3,15,51,52,53,61,62,63,120,121)
    write(*,*)
    write(*,'("Error(init0): spin-spirals do not work with task ",I4)') task
    write(*,*)
    stop
  end select
  if (xctype.lt.0) then
    write(*,*)
    write(*,'("Error(init0): spin-spirals do not work with the OEP method")')
    write(*,*)
    stop
  end if
end if
! spin-orbit coupling or fixed spin moment implies spin-polarised calculation
if ((spinorb).or.(fixspin).or.(spinsprl)) spinpol=.true.
! number of spinor components
if (spinpol) then
  nspinor=2
else
  nspinor=1
end if
! number of spin-dependent first-variational functions per state
if (spinsprl) then
  nspnfv=2
else
  nspnfv=1
end if
! spin-polarised calculations require second-variational eigenvectors
if (spinpol) tevecsv=.true.
! Hartree-Fock/RDMFT requires second-variational eigenvectors
if ((task.eq.5).or.(task.eq.300)) tevecsv=.true.
! get exchange-correlation functional data
call getxcdata(xctype,xcdescr,xcspin,xcgrad)
if ((spinpol).and.(xcspin.eq.0)) then
  write(*,*)
  write(*,'("Error(init0): requested spin-polarised run with &
   &spin-unpolarised")')
  write(*,'(" exchange-correlation functional")')
  write(*,*)
  stop
end if
! check for collinearity in the z-direction and set the dimension of the
! magnetisation and exchange-correlation vector fields
if (spinpol) then
  ndmag=1
  if ((abs(bfieldc(1)).gt.epslat).or.(abs(bfieldc(2)).gt.epslat)) ndmag=3
  do is=1,nspecies
    do ia=1,natoms(is)
      if ((abs(bfcmt(1,ia,is)).gt.epslat).or.(abs(bfcmt(2,ia,is)).gt.epslat)) &
       ndmag=3
    end do
  end do
! source-free fields and spin-spirals are non-collinear in general
  if ((nosource).or.(spinsprl)) ndmag=3
else
  ndmag=0
end if
! set fixed spin moment effective field to zero
bfsmc(:)=0.d0

!-------------------------------------!
!     lattice and symmetry set up     !
!-------------------------------------!
! generate the reciprocal lattice vectors and unit cell volume
call reciplat
! compute the inverse of the lattice vector matrix
call r3minv(avec,ainv)
! compute the inverse of the reciprocal vector matrix
call r3minv(bvec,binv)
do is=1,nspecies
  do ia=1,natoms(is)
! map atomic lattice coordinates to [0,1) if not in molecule mode
    if (.not.molecule) call r3frac(epslat,atposl(1,ia,is),iv)
! determine atomic Cartesian coordinates
    call r3mv(avec,atposl(1,ia,is),atposc(1,ia,is))
! lattice coordinates of the muffin-tin magnetic fields
    call r3mv(ainv,bfcmt(1,ia,is),bflmt(1,ia,is))
  end do
end do
! lattice coordinates of the global magnetic field
call r3mv(ainv,bfieldc,bfieldl)
! Cartesian coordinates of the spin-spiral vector
call r3mv(bvec,vqlss,vqcss)
! find Bravais lattice symmetries
call findsymlat
! use only the identity if required
if (nosym) nsymlat=1
! find the crystal symmetries and shift atomic positions if required
call findsymcrys
! find the site symmetries
call findsymsite
! automatically determine the muffin-tin radii if required
if (autormt) call autoradmt
! check for overlapping muffin-tins
call checkmt

!-----------------------!
!     radial meshes     !
!-----------------------!
nrmtmax=1
nrcmtmax=1
spnrmax=1
if (nspecies.eq.0) then
  rmtmin=2.d0
  rmtmax=2.d0
else
  rmtmin=rmt(1)
  rmtmax=rmt(1)
end if
do is=1,nspecies
! make the muffin-tin mesh commensurate with lradstp
  nrmt(is)=nrmt(is)-mod(nrmt(is)-1,lradstp)
  nrmtmax=max(nrmtmax,nrmt(is))
! number of coarse radial mesh points
  nrcmt(is)=(nrmt(is)-1)/lradstp+1
  nrcmtmax=max(nrcmtmax,nrcmt(is))
! estimate the number of radial mesh points to infinity
  t1=dble(nrmt(is))*log(sprmax(is)/sprmin(is))/log(rmt(is)/sprmin(is))
  spnr(is)=max(nint(t1),nrmt(is))
  spnrmax=max(spnrmax,spnr(is))
! smallest and largest muffin-tin radii
  rmtmin=min(rmtmin,rmt(is))
  rmtmax=max(rmtmax,rmt(is))
end do
! set up atomic and muffin-tin radial meshes
call genrmesh

!--------------------------------------!
!     charges and number of states     !
!--------------------------------------!
chgzn=0.d0
chgcr=0.d0
chgval=0.d0
spnstmax=0
do is=1,nspecies
! nuclear charge
  chgzn=chgzn+spzn(is)*dble(natoms(is))
! find the maximum number of atomic states
  spnstmax=max(spnstmax,spnst(is))
! compute the electronic charge for each species, as well as the total core and
! valence charge
  spze(is)=0.d0
  do ist=1,spnst(is)
    spze(is)=spze(is)+spocc(ist,is)
    if (spcore(ist,is)) then
      chgcr=chgcr+dble(natoms(is))*spocc(ist,is)
    else
      chgval=chgval+dble(natoms(is))*spocc(ist,is)
    end if
  end do
end do
! add excess charge
chgval=chgval+chgexs
! total charge
chgtot=chgcr+chgval
if (chgtot.lt.1.d-8) then
  write(*,*)
  write(*,'("Error(init0): zero total charge")')
  write(*,*)
  stop
end if
! number of first-variational states
nstfv=int(chgval/2.d0)+nempty+1
! number of second-variational states
nstsv=nstfv*nspinor

!-------------------------!
!     G-vector arrays     !
!-------------------------!
! find the G-vector grid sizes
call gridsize
! generate the G-vectors
call gengvec
! generate the spherical harmonics of the G-vectors
call genylmg
! allocate structure factor array for G-vectors
if (allocated(sfacg)) deallocate(sfacg)
allocate(sfacg(ngvec,natmtot))
! generate structure factors for G-vectors
call gensfacgp(ngvec,vgc,ngvec,sfacg)
! generate the characteristic function
call gencfun

!-------------------------!
!     atoms and cores     !
!-------------------------!
! solve the Kohn-Sham-Dirac equations for all atoms
call allatoms
! allocate core state eigenvalue array and set to default
if (allocated(evalcr)) deallocate(evalcr)
allocate(evalcr(spnstmax,natmtot))
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist=1,spnst(is)
      evalcr(ist,ias)=speval(ist,is)
    end do
  end do
end do
! allocate core state radial wavefunction array
if (allocated(rwfcr)) deallocate(rwfcr)
allocate(rwfcr(spnrmax,2,spnstmax,natmtot))
! allocate core state charge density array
if (allocated(rhocr)) deallocate(rhocr)
allocate(rhocr(spnrmax,natmtot))

!---------------------------------------!
!     charge density and potentials     !
!---------------------------------------!
! allocate charge density arrays
if (allocated(rhomt)) deallocate(rhomt)
allocate(rhomt(lmmaxvr,nrmtmax,natmtot))
if (allocated(rhoir)) deallocate(rhoir)
allocate(rhoir(ngrtot))
! allocate magnetisation arrays
if (allocated(magmt)) deallocate(magmt)
if (allocated(magir)) deallocate(magir)
if (spinpol) then
  allocate(magmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(magir(ngrtot,ndmag))
end if
! Coulomb potential
if (allocated(vclmt)) deallocate(vclmt)
allocate(vclmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(vclir)) deallocate(vclir)
allocate(vclir(ngrtot))
! exchange-correlation potential
if (allocated(vxcmt)) deallocate(vxcmt)
allocate(vxcmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(vxcir)) deallocate(vxcir)
allocate(vxcir(ngrtot))
! exchange-correlation magnetic field
if (allocated(bxcmt)) deallocate(bxcmt)
if (allocated(bxcir)) deallocate(bxcir)
if (spinpol) then
  allocate(bxcmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(bxcir(ngrtot,ndmag))
end if
! exchange energy density
if (allocated(exmt)) deallocate(exmt)
allocate(exmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(exir)) deallocate(exir)
allocate(exir(ngrtot))
! correlation energy density
if (allocated(ecmt)) deallocate(ecmt)
allocate(ecmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(ecir)) deallocate(ecir)
allocate(ecir(ngrtot))
! effective potential
if (allocated(veffmt)) deallocate(veffmt)
allocate(veffmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(veffir)) deallocate(veffir)
allocate(veffir(ngrtot))
if (allocated(veffig)) deallocate(veffig)
allocate(veffig(ngvec))
! allocate muffin-tin charge and moment arrays
if (allocated(chgmt)) deallocate(chgmt)
allocate(chgmt(natmtot))
if (allocated(mommt)) deallocate(mommt)
allocate(mommt(3,natmtot))

!--------------------------------------------!
!     forces and structural optimisation     !
!--------------------------------------------!
if (allocated(forcehf)) deallocate(forcehf)
allocate(forcehf(3,natmtot))
if (allocated(forcecr)) deallocate(forcecr)
allocate(forcecr(3,natmtot))
if (allocated(forceibs)) deallocate(forceibs)
allocate(forceibs(3,natmtot))
if (allocated(forcetot)) deallocate(forcetot)
allocate(forcetot(3,natmtot))
if (allocated(forcetp)) deallocate(forcetp)
allocate(forcetp(3,natmtot))
if (allocated(tauatm)) deallocate(tauatm)
allocate(tauatm(natmtot))
! initialise the previous force
forcetp(:,:)=0.d0
! initial step sizes
tauatm(:)=tau0atm

!-----------------------!
!     miscellaneous     !
!-----------------------!
! determine the nuclear-nuclear energy
call energynn
! call LAPACK routines which contain the SAVE attribute
call dlartg(1.d0,0.d0,cs,sn,r)
t1=dlamch('E')
! get smearing function data
call getsdata(stype,sdescr)
! generate the spherical harmonic transform (SHT) matrices
call genshtmat
! allocate 1D plotting arrays
if (allocated(dvp1d)) deallocate(dvp1d)
allocate(dvp1d(nvp1d))
if (allocated(vplp1d)) deallocate(vplp1d)
allocate(vplp1d(3,npp1d))
if (allocated(dpp1d)) deallocate(dpp1d)
allocate(dpp1d(npp1d))
! zero self-consistent loop number
iscl=0

call cpu_time(cpu1)
timeinit=timeinit+cpu1-cpu0

return
end subroutine
!EOC

