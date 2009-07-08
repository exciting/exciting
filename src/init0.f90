

! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: init0
! !INTERFACE:


subroutine init0
! !USES:
use modinput
use modmain
use modxcifc
#ifdef XS
use modxs
#endif
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
integer::is, js, ia, ias
integer::ist, l, m, lm, iv(3)
real(8)::ts0, ts1 , tv3(3)

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
call timesec(ts0)

!------------------------------------!
!     angular momentum variables     !
!------------------------------------!
lmmaxvr=(input%groundstate%lmaxvr+1)**2
lmmaxapw=(input%groundstate%lmaxapw+1)**2
lmmaxmat=(input%groundstate%lmaxmat+1)**2
lmmaxinr=(input%groundstate%lmaxinr+1)**2
if (input%groundstate%lmaxvr.gt.input%groundstate%lmaxapw) then
  write(*, *)
  write(*, '("Error(init0): lmaxvr > lmaxapw : ", 2I8)') input%groundstate%lmaxvr, input%groundstate%lmaxapw
  write(*, *)
  stop
end if
if (input%groundstate%lmaxmat.gt.input%groundstate%lmaxapw) then
  write(*, *)
  write(*, '("Error(init0): lmaxmat > lmaxapw : ", 2I8)') input%groundstate%lmaxmat, input%groundstate%lmaxapw
  write(*, *)
  stop
end if
! index to (l,m) pairs
if (allocated(idxlm)) deallocate(idxlm)
allocate(idxlm(0:input%groundstate%lmaxapw, -input%groundstate%lmaxapw:input%groundstate%lmaxapw))
lm=0
do l=0, input%groundstate%lmaxapw
  do m=-l, l
    lm=lm+1
    idxlm(l, m)=lm
  end do
end do
! array of i**l values
if (allocated(zil)) deallocate(zil)
allocate(zil(0:input%groundstate%lmaxapw))
do l=0, input%groundstate%lmaxapw
  zil(l)=zi**l
end do

!------------------------------------!
!     index to atoms and species     !
!------------------------------------!
! check if the system is an isolated molecule
if (input%structure%molecule) then
  input%structure%primcell=.false.
  input%structure%tshift=.false.
end if
! find primitive cell if required
if (input%structure%primcell) call findprim
natmmax=0
ias=0
do is=1, nspecies
  do ia=1, natoms(is)
    ias=ias+1
    idxas(ia, is)=ias
  end do
! maximum number of atoms over all species
  natmmax=max(natmmax, natoms(is))
end do
! total number of atoms
natmtot=ias

!------------------------!
!     spin variables     !
!------------------------!
if (isspinspiral()) then
  select case(task)
  case(2, 3, 15, 51, 52, 53, 61, 62, 63, 120, 121)
    write(*, *)
    write(*, '("Error(init0): spin-spirals do not work with task ", I4)') task
    write(*, *)
    stop
  end select
  if (input%groundstate%xctypenumber.lt.0) then
    write(*, *)
    write(*, '("Error(init0): spin-spirals do not work with the OEP method")')
    write(*, *)
    stop
  end if
end if
! spin-orbit coupling or fixed spin moment implies spin-polarised calculation

! number of spinor components and maximum allowed occupancy
if (associated(input%groundstate%spin)) then
  nspinor=2
  occmax=1.d0
else
  nspinor=1
  occmax=2.d0
end if
! number of spin-dependent first-variational functions per state
if (isspinspiral()) then
  nspnfv=2
else
  nspnfv=1
end if
! spin-polarised calculations require second-variational eigenvectors
if (associated(input%groundstate%spin)) input%groundstate%tevecsv=.true.
! Hartree-Fock/RDMFT requires second-variational eigenvectors
if ((task.eq.5).or.(task.eq.6).or.(task.eq.300)) input%groundstate%tevecsv=.true.
! get exchange-correlation functional data
call getxcdata(input%groundstate%xctypenumber, xcdescr, xcspin, xcgrad)
if ((associated(input%groundstate%spin)).and.(xcspin.eq.0)) then
  write(*, *)
  write(*, '("Error(init0): requested spin-polarised run with &
   &spin-unpolarised")')
  write(*, '(" exchange-correlation functional")')
  write(*, *)
  stop
end if
! check for collinearity in the z-direction and set the dimension of the
! magnetisation and exchange-correlation vector fields
if (associated(input%groundstate%spin)) then
  ndmag=1
  if&
    &((abs(input%groundstate%spin%bfieldc(1)).gt.input%structure%epslat).or.(abs(input%groundstate%spin%bfieldc(2)).gt.&
    &input%structure%epslat)) ndmag = 3
  do is=1, nspecies
    do ia=1, natoms(is)
      if&
    &((abs(input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(1)).gt.input%structure%epslat).or.(abs(inp&
    &ut%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(2)).gt.input%structure%epslat)) &
       ndmag = 3
    end do
  end do
! source-free fields and spin-spirals are non-collinear in general
  if ((nosource).or.(isspinspiral())) ndmag=3
! spin-orbit coupling is non-collinear in general
  if (isspinorb()) ndmag=3
else
  ndmag=0
end if
! set the non-collinear flag
if (ndmag.eq.3) then
  ncmag=.true.
else
  ncmag=.false.
end if
if ((ncmag).and.(xcgrad.gt.0)) then
  write(*, *)
  write(*, '("Warning(init0): GGA inconsistent with non-collinear magnetism")')
end if
! set fixed spin moment effective field to zero
bfsmc(:)=0.d0
! set muffin-tin FSM fields to zero
bfsmcmt(:, :, :)=0.d0

!-------------------------------------!
!     lattice and symmetry set up     !
!-------------------------------------!
! generate the reciprocal lattice vectors and unit cell volume
call reciplat
! compute the inverse of the lattice vector matrix
call r3minv(input%structure%crystal%basevect, ainv)
! compute the inverse of the reciprocal vector matrix
call r3minv(bvec, binv)
do is=1, nspecies
  do ia=1, natoms(is)
! map atomic lattice coordinates to [0,1) if not in molecule mode
    if (.not.input%structure%molecule) call r3frac(input%structure%epslat, &
    &input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), iv)
! determine atomic Cartesian coordinates
    call r3mv(input%structure%crystal%basevect, &
    &input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), atposc(:, ia, is))
! lattice coordinates of the muffin-tin magnetic fields
    call r3mv(ainv, input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:), bflmt(:, ia, is))
  end do
end do
! lattice coordinates of the global magnetic field
if(associated(input%groundstate%spin))then
	tv3=input%groundstate%spin%bfieldc
else
	tv3=0
endif
call r3mv(ainv, tv3, bfieldl)
! Cartesian coordinates of the spin-spiral vector
if(associated(input%groundstate%spin))then
	tv3=input%groundstate%spin%vqlss
else
	tv3=0
endif
call r3mv(bvec, tv3, vqcss)
! find Bravais lattice symmetries
call findsymlat
! use only the identity if required
if (input%groundstate%nosym) nsymlat=1
! find the crystal symmetries and shift atomic positions if required
call findsymcrys
! find the site symmetries
call findsymsite
#ifdef XS
! determine inverse symmery elements
call findsymi(input%structure%epslat, maxsymcrys, nsymcrys, symlat, lsplsymc, vtlsymc, isymlat, &
     scimap)
! generate symmetrization array for rank 2 tensors
call gensymt2(maxsymcrys, nsymcrys, symlatc, lsplsymc, symt2)
! calculate advanced information on symmetry group
call setupsym
#endif
! automatically determine the muffin-tin radii if required
if (input%structure%autormt) call autoradmt
! check for overlapping muffin-tins
call checkmt

!-----------------------!
!     radial meshes     !
!-----------------------!
nrmtmax=1
nrcmtmax=1
js=1
do is=1, nspecies
! make the muffin-tin mesh commensurate with lradstp
  nrmt(is)=nrmt(is)-mod(nrmt(is)-1, input%groundstate%lradstep)
  nrmtmax=max(nrmtmax, nrmt(is))
! number of coarse radial mesh points
  nrcmt(is)=(nrmt(is)-1)/input%groundstate%lradstep+1
  nrcmtmax=max(nrcmtmax, nrcmt(is))
! smallest muffin-tin radius
  if (rmt(is).lt.rmt(js)) js=is
end do
if ((input%groundstate%isgkmax.lt.1).or.(input%groundstate%isgkmax.gt.nspecies)) input%groundstate%isgkmax = js
! set up atomic and muffin-tin radial meshes
call genrmesh

!--------------------------------------!
!     charges and number of states     !
!--------------------------------------!
chgzn=0.d0
chgcr=0.d0
chgval=0.d0
spnstmax=0
do is=1, nspecies
! nuclear charge
  chgzn=chgzn+spzn(is)*dble(natoms(is))
! find the maximum number of atomic states
  spnstmax=max(spnstmax, spnst(is))
! compute the electronic charge for each species, as well as the total core and
! valence charge
  spze(is)=0.d0
  do ist=1, spnst(is)
    spze(is)=spze(is)+spocc(ist, is)
    if (spcore(ist, is)) then
      chgcr=chgcr+dble(natoms(is))*spocc(ist, is)
    else
      chgval=chgval+dble(natoms(is))*spocc(ist, is)
    end if
  end do
end do
! add excess charge
chgval=chgval+input%groundstate%chgexs
! total charge
chgtot=chgcr+chgval
if (chgtot.lt.1.d-8) then
  write(*, *)
  write(*, '("Error(init0): zero total charge")')
  write(*, *)
  stop
end if
! effective Wigner radius
rwigner=(3.d0/(fourpi*(chgtot/omega)))**(1.d0/3.d0)

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
allocate(sfacg(ngvec, natmtot))
! generate structure factors for G-vectors
call gensfacgp(ngvec, vgc, ngvec, sfacg)
! generate the characteristic function
call gencfun

!-------------------------!
!     atoms and cores     !
!-------------------------!
#ifdef XS
if (init0symonly) goto 10
#endif
! solve the Kohn-Sham-Dirac equations for all atoms
call allatoms
! allocate core state eigenvalue array and set to default
if (allocated(evalcr)) deallocate(evalcr)
allocate(evalcr(spnstmax, natmtot))
do is=1, nspecies
  do ia=1, natoms(is)
    ias=idxas(ia, is)
    do ist=1, spnst(is)
      evalcr(ist, ias)=speval(ist, is)
    end do
  end do
end do
! allocate core state radial wavefunction array
if (allocated(rwfcr)) deallocate(rwfcr)
allocate(rwfcr(spnrmax, 2, spnstmax, natmtot))
! allocate core state charge density array
if (allocated(rhocr)) deallocate(rhocr)
allocate(rhocr(spnrmax, natmtot))
#ifdef XS
10 continue
#endif

!---------------------------------------!
!     charge density and potentials     !
!---------------------------------------!
! allocate charge density arrays
if (allocated(rhomt)) deallocate(rhomt)
allocate(rhomt(lmmaxvr, nrmtmax, natmtot))
if (allocated(rhoir)) deallocate(rhoir)
allocate(rhoir(ngrtot))
! allocate magnetisation arrays
if (allocated(magmt)) deallocate(magmt)
if (allocated(magir)) deallocate(magir)
if (associated(input%groundstate%spin)) then
  allocate(magmt(lmmaxvr, nrmtmax, natmtot, ndmag))
  allocate(magir(ngrtot, ndmag))
end if
! Coulomb potential
if (allocated(vclmt)) deallocate(vclmt)
allocate(vclmt(lmmaxvr, nrmtmax, natmtot))
if (allocated(vclir)) deallocate(vclir)
allocate(vclir(ngrtot))
! exchange-correlation potential
if (allocated(vxcmt)) deallocate(vxcmt)
allocate(vxcmt(lmmaxvr, nrmtmax, natmtot))
if (allocated(vxcir)) deallocate(vxcir)
allocate(vxcir(ngrtot))
! exchange-correlation magnetic field
if (allocated(bxcmt)) deallocate(bxcmt)
if (allocated(bxcir)) deallocate(bxcir)
if (associated(input%groundstate%spin)) then
  allocate(bxcmt(lmmaxvr, nrmtmax, natmtot, ndmag))
  allocate(bxcir(ngrtot, ndmag))
end if
! exchange energy density
if (allocated(exmt)) deallocate(exmt)
allocate(exmt(lmmaxvr, nrmtmax, natmtot))
if (allocated(exir)) deallocate(exir)
allocate(exir(ngrtot))
! correlation energy density
if (allocated(ecmt)) deallocate(ecmt)
allocate(ecmt(lmmaxvr, nrmtmax, natmtot))
if (allocated(ecir)) deallocate(ecir)
allocate(ecir(ngrtot))
! effective potential
if (allocated(veffmt)) deallocate(veffmt)
allocate(veffmt(lmmaxvr, nrmtmax, natmtot))
if (allocated(veffir)) deallocate(veffir)
allocate(veffir(ngrtot))
if (allocated(veffig)) deallocate(veffig)
allocate(veffig(ngvec))
! allocate muffin-tin charge and moment arrays
if (allocated(chgmt)) deallocate(chgmt)
allocate(chgmt(natmtot))
if (allocated(mommt)) deallocate(mommt)
allocate(mommt(3, natmtot))

!--------------------------------------------!
!     forces and structural optimisation     !
!--------------------------------------------!
if (allocated(forcehf)) deallocate(forcehf)
allocate(forcehf(3, natmtot))
if (allocated(forcecr)) deallocate(forcecr)
allocate(forcecr(3, natmtot))
if (allocated(forceibs)) deallocate(forceibs)
allocate(forceibs(3, natmtot))
if (allocated(forcetot)) deallocate(forcetot)
allocate(forcetot(3, natmtot))
if (allocated(forcetp)) deallocate(forcetp)
allocate(forcetp(3, natmtot))
if (allocated(tauatm)) deallocate(tauatm)
allocate(tauatm(natmtot))
! initialise the previous force
forcetp(:, :)=0.d0
! initial step sizes
if(associated(input%structureoptimization))then
tauatm(:)=input%structureoptimization%tau0atm
else
tauatm(:)=0
endif

!-------------------------!
!     LDA+U variables     !
!-------------------------!
if ((ldapu.ne.0).or.(task.eq.17)) then
! LDA+U requires second-variational eigenvectors
  input%groundstate%tevecsv=.true.
! density matrices
  if (allocated(dmatlu)) deallocate(dmatlu)
  allocate(dmatlu(lmmaxlu, lmmaxlu, nspinor, nspinor, natmtot))
! potential matrix elements
  if (allocated(vmatlu)) deallocate(vmatlu)
  allocate(vmatlu(lmmaxlu, lmmaxlu, nspinor, nspinor, natmtot))
! zero the potential
  vmatlu(:, :, :, :, :)=0.d0
! energy for each atom
  if (allocated(engyalu)) deallocate(engyalu)
  allocate(engyalu(natmtot))
! interpolation constants (alpha)
  if (allocated(alphalu)) deallocate(alphalu)
  allocate(alphalu(natmtot))
end if

!-----------------------!
!     miscellaneous     !
!-----------------------!
! determine the nuclear-nuclear energy
call energynn
! get smearing function data
call getsdata(input%groundstate%stypenumber, sdescr)
! generate the spherical harmonic transform (SHT) matrices
call genshtmat

! allocate 1D plotting arrays
if (allocated(dvp1d)) deallocate(dvp1d)
allocate(dvp1d(nvp1d))
if (allocated(vplp1d)) deallocate(vplp1d)
allocate(vplp1d(3, npp1d))
if (allocated(dpp1d)) deallocate(dpp1d)
allocate(dpp1d(npp1d))
! zero self-consistent loop number
iscl=0
tlast=.false.

call timesec(ts1)
timeinit=timeinit+ts1-ts0

return
end subroutine
!EOC
