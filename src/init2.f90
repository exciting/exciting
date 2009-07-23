


! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine init2
use modmain
use modinput
#ifdef XS
use modxs
#endif
implicit none
! local variables
integer::is, ia, ist, ic, m
real(8)::ts0, ts1
real(8)::boxl(3, 4)
#ifdef XS
real(8) :: v(3),t1
integer :: iq,iv(3)
#endif

call timesec(ts0)

!---------------------!
!     q-point set     !
!---------------------!
! check if the system is an isolated molecule
if (input%structure%molecule) ngridq(:)=1
! OEP, Hartree-Fock or RDMFT
if ((input%groundstate%xctypenumber.lt.0).or.(task.eq.5).or.(task.eq.6).or.(task.eq.300)) then
  ngridq(:)=input%groundstate%ngkgrid(:)
  input%phonons%reduceq=.false.
end if
#ifdef XS
if (task.le.300) then
#endif
! allocate the q-point arrays
if (allocated(ivq)) deallocate(ivq)
allocate(ivq(3,ngridq(1)*ngridq(2)*ngridq(3)))
if (allocated(vql)) deallocate(vql)
allocate(vql(3,ngridq(1)*ngridq(2)*ngridq(3)))
if (allocated(vqc)) deallocate(vqc)
allocate(vqc(3,ngridq(1)*ngridq(2)*ngridq(3)))
if (allocated(wqpt)) deallocate(wqpt)
allocate(wqpt(ngridq(1)*ngridq(2)*ngridq(3)))
if (allocated(iqmap)) deallocate(iqmap)
allocate(iqmap(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
! setup the q-point box (offset should always be zero)
boxl(:,:)=0.d0
boxl(1,2)=1.d0; boxl(2,3)=0.d0; boxl(3,4)=0.d0
! generate the q-point set, note that the vectors vql and vqc are mapped to the
! first Brillouin zone
call genppts(input%phonons%reduceq, .true., ngridq, boxl, nqpt, iqmap, ivq, vql, vqc, wqpt)
#ifdef XS
end if
#endif
#ifdef XS
! Q-/q-point set should have no offset
! setup the q-point box (offset should always be zero)
boxl(:,:)=0.d0
boxl(1,2)=1.d0; boxl(2,3)=1.d0; boxl(3,4)=1.d0
boxl(:,2)=boxl(:,2)+boxl(:,1)
boxl(:,3)=boxl(:,3)+boxl(:,1)
boxl(:,4)=boxl(:,4)+boxl(:,1)
! assign momentum transfer Q-points set to q-point set
if (((task.ge.301).and.(task.le.399))) then
   nqpt=size(input%xs%qpointset%qpoint,1)
   if (allocated(vqlmt)) deallocate(vqlmt)
   allocate(vqlmt(3,nqpt))
   if (allocated(ivgmt)) deallocate(ivgmt)
   allocate(ivgmt(3,nqpt))
   if (allocated(vql)) deallocate(vql)
   allocate(vql(3,nqpt))
   if (allocated(vqc)) deallocate(vqc)
   allocate(vqc(3,nqpt))
   do iq=1,nqpt
      v(:)=input%xs%qpointset%qpoint(:,iq)
      iv(:)=0
      ! map Q-point to reciprocal unit cell
      if (input%xs%tddft%mdfqtype.eq.1) call r3frac(input%structure%epslat, v, iv)
      vqlmt(:,iq)=v(:)
      ivgmt(:,iq)=iv(:)
      vql(:,iq)=vqlmt(:,iq)
      vqc(:,iq)=vql(1,iq)*bvec(:,1)+vql(2,iq)*bvec(:,2)+ &
	   vql(3,iq)*bvec(:,3)
   end do
else
   ! determine only integer-part of Q-points
   if (allocated(ivgmt)) deallocate(ivgmt)
   allocate(ivgmt(3,size(input%xs%qpointset%qpoint,1)))
   do iq=1,size(input%xs%qpointset%qpoint,1)
      v(:)=input%xs%qpointset%qpoint(:,iq)
      iv(:)=0
      ! map Q-point to reciprocal unit cell
      if (input%xs%tddft%mdfqtype.eq.1) call r3frac(input%structure%epslat, v, iv)
      ivgmt(:,iq)=iv(:)
   end do
end if
! generate q-point set from grid
if ((task.ge.400).and.(task.le.439)) then
   if (allocated(ivq)) deallocate(ivq)
   allocate(ivq(3,ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(vql)) deallocate(vql)
   allocate(vql(3,ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(vqc)) deallocate(vqc)
   allocate(vqc(3,ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(wqpt)) deallocate(wqpt)
   allocate(wqpt(ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(iqmap)) deallocate(iqmap)
   allocate(iqmap(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
   ! generate reduced q-point set
   call genppts(input%xs%reduceq, input%xs%BSE%fbzq, ngridq, boxl, nqpt, iqmap, ivq, vql, vqc, wqpt)
   nqptr=nqpt
end if
if ((task.eq.440).or.(task.eq.441).or.(task.eq.445).or.(task.eq.450).or.(task.eq.451).or. &
     (task.eq.499).or.(task.eq.700)) then
   if (allocated(ivqr)) deallocate(ivqr)
   allocate(ivqr(3,ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(vqlr)) deallocate(vqlr)
   allocate(vqlr(3,ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(vqcr)) deallocate(vqcr)
   allocate(vqcr(3,ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(wqptr)) deallocate(wqptr)
   allocate(wqptr(ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(iqmapr)) deallocate(iqmapr)
   allocate(iqmapr(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
   ! generate reduced q-point set
   call genppts(input%xs%reduceq, input%xs%BSE%fbzq, ngridq, boxl, nqptr, iqmapr, ivqr, vqlr, vqcr, &
      wqptr)
   if (allocated(ivq)) deallocate(ivq)
   allocate(ivq(3,ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(vql)) deallocate(vql)
   allocate(vql(3,ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(vqc)) deallocate(vqc)
   allocate(vqc(3,ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(wqpt)) deallocate(wqpt)
   allocate(wqpt(ngridq(1)*ngridq(2)*ngridq(3)))
   if (allocated(iqmap)) deallocate(iqmap)
   allocate(iqmap(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
   ! generate non-reduced q-point set
   call genppts(.false., input%xs%BSE%fbzq, ngridq, boxl, nqpt, iqmap, ivq, vql, vqc, wqpt)
end if
! find (little/small) group of q
if (allocated(nsymcrysq)) deallocate(nsymcrysq)
allocate(nsymcrysq(nqpt))
if (allocated(scqmap)) deallocate(scqmap)
allocate(scqmap(nsymcrys,nqpt))
if (allocated(ivscwrapq)) deallocate(ivscwrapq)
allocate(ivscwrapq(3,nsymcrys,nqpt))
!do iq=1,nqpt
!   call findgroupq(fbzq,vql(1,iq),epslat,bvec,symlat,nsymcrys,lsplsymc,&
!	nsymcrysq(iq),scqmap(1,iq),ivscwrapq(1,1,iq))
!end do

!-----------------------!
!     k+q-point set	!
!-----------------------!
if (allocated(qvkloff)) deallocate(qvkloff)
allocate(qvkloff(3,0:nqpt))
if (allocated(ikmapikq)) deallocate(ikmapikq)
allocate(ikmapikq(nkpt,nqpt))
qvkloff(:, 0)=input%groundstate%vkloff(:)
do iq=1,nqpt
   ! offset for k+q-point set derived from q-point
   call genqvkloff(vql(1,iq),qvkloff(1,iq))
   ! map from k-point index to k+q point index for same k
   call findkmapkq(vql(1,iq),qvkloff(1,iq),ikmapikq(1,iq))
end do

!---------------------!
!     G+q-point set   !
!---------------------!
! checking
if (input%xs%gqmax.ge.gkmax) then
   write( * , '(a, 2g18.10)') 'Warning(init2/xs): input%xs%gqmax >= gkmax: ', input%xs%gqmax, &
	gkmax
end if
! maximum number of G+q vectors for all q
call getngqmax

! allocate the G+q-vector arrays
if (allocated(ngq)) deallocate(ngq)
allocate(ngq(nqpt))
if (allocated(igqig)) deallocate(igqig)
allocate(igqig(ngqmax,nqpt))
if (allocated(vgql)) deallocate(vgql)
allocate(vgql(3,ngqmax,nqpt))
if (allocated(vgqc)) deallocate(vgqc)
allocate(vgqc(3,ngqmax,nqpt))
if (allocated(gqc)) deallocate(gqc)
allocate(gqc(ngqmax,nqpt))
if (allocated(tpgqc)) deallocate(tpgqc)
allocate(tpgqc(2,ngqmax,nqpt))
if (allocated(sfacgq)) deallocate(sfacgq)
allocate(sfacgq(ngqmax,natmtot,nqpt))
if (allocated(ylmgq)) deallocate(ylmgq)
allocate(ylmgq(lmmaxapw,ngqmax,nqpt))
if (allocated(ivgigq)) deallocate(ivgigq)
allocate(ivgigq(intgqv(1,1):intgqv(1,2),intgqv(2,1):intgqv(2,2), &
     intgqv(3,1):intgqv(3,2),nqpt))
do iq=1,nqpt
   ! generate G+q vectors
   call gengqvec(iq,vql(1,iq),vqc(1,iq),ngq(iq),igqig(1,iq), &
	vgql(1,1,iq),vgqc(1,1,iq),gqc(1,iq),tpgqc(1,1,iq))
   ! generate structure factors for G-vectors
   call gensfacgp(ngq(iq),vgqc(1,1,iq),ngqmax,sfacgq(1,1,iq))
   ! spherical harmonics for G+q-vectors
   call genylmgq(iq, input%groundstate%lmaxvr)
end do
t1=maxval(gqc(1,:))
!if ((task.ge.400).and.(task.le.499).and.(gqmax.lt.t1)) then
!   write(*,*)
!   write(*,'("Error(init2): gqmax too small:",g18.10)') gqmax
!   write(*,'(" smalles possible value:",g18.10)') t1
!   write(*,*)
!   call terminate
!end if

!---------------------------!
!     Coulomb potential     !
!---------------------------!
if (allocated(sptclg)) deallocate(sptclg)
allocate(sptclg(ngqmax,nqpt))
do iq=1,nqpt
   ! set up Coulomb potential square root
   call genptclg('nocutoff',ngqmax,ngq(iq),vgqc(:,:,iq),gqc(:,iq),sptclg(:,iq))
end do

!------------------------!
!     radial functions   !
!------------------------!
! read density and potentials from file (STATE.OUT) exclusively
isreadstate0=.true.
call readstate
isreadstate0=.false.
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
#endif

!-----------------------------------------------!
!     OEP, Hartree-Fock and RDMFT variables     !
!-----------------------------------------------!
if ((input%groundstate%xctypenumber.lt.0).or.(task.eq.5).or.(task.eq.6).or.(task.eq.300)) then
! determine the 1/q^2 integral weights if required
  call genwiq2
! output the 1/q^2 integrals to WIQ2.OUT
  call writewiq2
end if
if (input%groundstate%xctypenumber.lt.0) then
! initialise OEP residual magnitude
  resoep=1.d0
! find maximum core states over all species
  ncrmax=0
  do is=1,nspecies
    do ia=1,natoms(is)
      ic=0
      do ist=1,spnst(is)
        if (spcore(ist,is)) then
          do m=-spk(ist,is),spk(ist,is)-1
            ic=ic+1
          end do
        end if
      end do
      ncrmax=max(ncrmax,ic)
    end do
  end do
! allocate and zero the complex exchange potential and field
  if (allocated(zvxmt)) deallocate(zvxmt)
  allocate(zvxmt(lmmaxvr,nrcmtmax,natmtot))
  zvxmt(:,:,:)=0.d0
  if (allocated(zvxir)) deallocate(zvxir)
  allocate(zvxir(ngrtot))
  zvxir(:)=0.d0
  if (associated(input%groundstate%spin)) then
    if (allocated(zbxmt)) deallocate(zbxmt)
    allocate(zbxmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
    zbxmt(:,:,:,:)=0.d0
    if (allocated(zbxir)) deallocate(zbxir)
    allocate(zbxir(ngrtot,ndmag))
    zbxir(:,:)=0.d0
  end if
end if
if ((task.eq.5).or.(task.eq.6).or.(task.eq.300)) then
! allocate the kinetic matrix elements for Hartree-Fock/RDMFT
  if (allocated(kinmatc)) deallocate(kinmatc)
  allocate(kinmatc(nstsv,nstsv,nkpt))
end if
if (task.eq.300) then
  if (allocated(vclmat)) deallocate(vclmat)
  allocate(vclmat(nstsv,nstsv,nkpt))
  if (allocated(dkdc)) deallocate(dkdc)
  allocate(dkdc(nstsv,nstsv,nkpt))
  if (allocated(vnlrdm)) deallocate(vnlrdm)
  allocate(vnlrdm(nstsv,nkpt,nstsv,nkptnr))
end if

call timesec(ts1)
timeinit=timeinit+ts1-ts0

return
end subroutine
