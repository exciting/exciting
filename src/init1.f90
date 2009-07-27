


! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: init1
! !INTERFACE:


subroutine init1
! !USES:
use modinput
use modmain
#ifdef TETRA
use modtetra
#endif
#ifdef XS
use modxs
#endif
! !DESCRIPTION:
!   Generates the $k$-point set and then allocates and initialises global
!   variables which depend on the $k$-point set.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer::ik, is, ia, ias, io, ilo
integer::i1, i2, i3, ispn, iv(3)
integer::l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
real(8)::vl(3), vc(3), boxl(3, 4)
real(8)::ts0, ts1
! external functions
complex(8) gauntyry
external gauntyry

call timesec(ts0)

!---------------------!
!     k-point set     !
!---------------------!
! check if the system is an isolated molecule
if (input%structure%molecule) then
  input%groundstate%ngkgrid(:)=1
  input%groundstate%vkloff(:)=0.d0
  input%groundstate%autokpt=.false.
end if
! setup the default k-point box
boxl(:, 1)=input%groundstate%vkloff(:)/dble(input%groundstate%ngkgrid(:))
boxl(:, 2)=boxl(:, 1); boxl(:, 3)=boxl(:, 1); boxl(:, 4)=boxl(:, 1)
boxl(1, 2)=boxl(1, 2)+1.d0
boxl(2, 3)=boxl(2, 3)+1.d0
boxl(3, 4)=boxl(3, 4)+1.d0
! k-point set and box for Fermi surface plots
if ((task.eq.100).or.(task.eq.101)) then
  input%groundstate%ngkgrid(:)=np3d(:)
  boxl(:, :)=vclp3d(:, :)
end if
if ((task.eq.20).or.(task.eq.21)) then
! for band structure plots generate k-points along a line
  nvp1d=size(input%properties%bandstructure%plot1d%path%pointarray)
  npp1d=input%properties%bandstructure%plot1d%path%steps
  if (allocated(dvp1d)) deallocate(dvp1d)
  allocate(dvp1d(nvp1d))
  if (allocated(vplp1d)) deallocate(vplp1d)
  allocate(vplp1d(3, npp1d))
  if (allocated(dpp1d)) deallocate(dpp1d)
  allocate(dpp1d(npp1d))
  call connect(bvec, input%properties%bandstructure%plot1d, nvp1d,  npp1d, vplp1d, dvp1d, dpp1d)

  nkpt=input%properties%bandstructure%plot1d%path%steps
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3, nkpt))
  if (allocated(vkc)) deallocate(vkc)
  allocate(vkc(3, nkpt))
  do ik=1, nkpt
    vkl(:, ik)=vplp1d(:, ik)
    call r3mv(bvec, vkl(:, ik), vkc(:, ik))
  end do
else if (task.eq.25) then
! effective mass calculation
  nkpt=(2*input%properties%masstensor%ndspem+1)**3
  if (allocated(ivk)) deallocate(ivk)
  allocate(ivk(3, nkpt))
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3, nkpt))
  if (allocated(vkc)) deallocate(vkc)
  allocate(vkc(3, nkpt))
! map vector to [0,1)
  call r3frac(input%structure%epslat, input%properties%masstensor%vklem, iv)
  ik=0
  do i3=-input%properties%masstensor%ndspem, input%properties%masstensor%ndspem
    do i2=-input%properties%masstensor%ndspem, input%properties%masstensor%ndspem
      do i1=-input%properties%masstensor%ndspem, input%properties%masstensor%ndspem
	ik=ik+1
	ivk(1, ik)=i1; ivk(2, ik)=i2; ivk(3, ik)=i3
	vc(1)=dble(i1); vc(2)=dble(i2); vc(3)=dble(i3)
	vc(:)=vc(:)*input%properties%masstensor%deltaem
	call r3mv(binv, vc, vl)
	vkl(:, ik)=input%properties%masstensor%vklem(:)+vl(:)
	call r3mv(bvec, vkl(:, ik), vkc(:, ik))
      end do
    end do
  end do
else
! determine the k-point grid automatically from radkpt if required
  if (input%groundstate%autokpt) then
    input%groundstate%ngkgrid(:) = int(input%groundstate%radkpt/sqrt(input%structure%crystal%basevect(1, :) ** 2 +&
    &input%structure%crystal%basevect(2, :) ** 2 + input%structure%crystal%basevect(3, :) ** 2)) + 1
  end if
! allocate the reduced k-point set arrays
  if (allocated(ivk)) deallocate(ivk)
  allocate(ivk(3, input%groundstate%ngkgrid(1) * input%groundstate%ngkgrid(2) * input%groundstate%ngkgrid(3)))
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3, input%groundstate%ngkgrid(1) * input%groundstate%ngkgrid(2) * input%groundstate%ngkgrid(3)))
  if (allocated(vkc)) deallocate(vkc)
  allocate(vkc(3, input%groundstate%ngkgrid(1) * input%groundstate%ngkgrid(2) * input%groundstate%ngkgrid(3)))
  if (allocated(wkpt)) deallocate(wkpt)
  allocate(wkpt(input%groundstate%ngkgrid(1) * input%groundstate%ngkgrid(2) * input%groundstate%ngkgrid(3)))
  if (allocated(ikmap)) deallocate(ikmap)
  allocate(ikmap(0:input%groundstate%ngkgrid(1) - 1, 0:input%groundstate%ngkgrid(2) - 1, &
    &0:input%groundstate%ngkgrid(3)&
    &- 1))
! generate the reduced k-point set
  call genppts(input%groundstate%reducek, .false., input%groundstate%ngkgrid, boxl, nkpt, ikmap, ivk, vkl, vkc, wkpt)
! allocate the non-reduced k-point set arrays
  nkptnr=input%groundstate%ngkgrid(1)*input%groundstate%ngkgrid(2)*input%groundstate%ngkgrid(3)
  if (allocated(ivknr)) deallocate(ivknr)
  allocate(ivknr(3, nkptnr))
  if (allocated(vklnr)) deallocate(vklnr)
  allocate(vklnr(3, nkptnr))
  if (allocated(vkcnr)) deallocate(vkcnr)
  allocate(vkcnr(3, nkptnr))
  if (allocated(wkptnr)) deallocate(wkptnr)
  allocate(wkptnr(nkptnr))
  if (allocated(ikmapnr)) deallocate(ikmapnr)
  allocate(ikmapnr(0:input%groundstate%ngkgrid(1) - 1, 0:input%groundstate%ngkgrid(2) - 1, &
    &0:input%groundstate%ngkgrid(3) - 1))
! generate the non-reduced k-point set
  call genppts(.false., .false., input%groundstate%ngkgrid, boxl, nkptnr, ikmapnr, ivknr, vklnr, vkcnr, &
   wkptnr)
#ifdef TETRA
  ! call to module routine
   if(associated(input%xs)) then
  if(associated(input%xs%tetra)) then
 	 if (input%xs%tetra%tetraocc.or.tetraopt.or.input%xs%tetra%tetradf)then
 	  call genkpts_tet(filext, input%structure%epslat, bvec, &
       maxsymcrys, nsymcrys, lsplsymc, symlat, input%groundstate%reducek, input%groundstate%ngkgrid, &
       input%groundstate%vkloff, nkpt, ikmap,  vkl, wkpt)
     endif
  endif
  endif
#endif
#ifdef XS
  ! determine inverse symmery elements
  call findsymi(input%structure%epslat, maxsymcrys, nsymcrys, symlat, lsplsymc, vtlsymc, isymlat, &
       scimap)
#endif
end if

!---------------------!
!     G+k vectors     !
!---------------------!
! determine gkmax
if ((input%groundstate%isgkmax.ge.1).and.(input%groundstate%isgkmax.le.nspecies)) then
  gkmax=input%groundstate%rgkmax/rmt(input%groundstate%isgkmax)
else
  gkmax=input%groundstate%rgkmax/2.d0
end if
if (2.d0*gkmax.gt.input%groundstate%gmaxvr+input%structure%epslat) then
  write(*, *)
  write(*, '("Error(init1): 2*gkmax > gmaxvr  ", 2G18.10)') 2.d0*gkmax, input%groundstate%gmaxvr
  write(*, *)
  stop
end if
! find the maximum number of G+k-vectors
call getngkmax
! allocate the G+k-vector arrays
if (allocated(ngk)) deallocate(ngk)
allocate(ngk(nspnfv, nkpt))
if (allocated(igkig)) deallocate(igkig)
allocate(igkig(ngkmax, nspnfv, nkpt))
if (allocated(vgkl)) deallocate(vgkl)
allocate(vgkl(3, ngkmax, nspnfv, nkpt))
if (allocated(vgkc)) deallocate(vgkc)
allocate(vgkc(3, ngkmax, nspnfv, nkpt))
if (allocated(gkc)) deallocate(gkc)
allocate(gkc(ngkmax, nspnfv, nkpt))
if (allocated(tpgkc)) deallocate(tpgkc)
allocate(tpgkc(2, ngkmax, nspnfv, nkpt))
if (allocated(sfacgk)) deallocate(sfacgk)
allocate(sfacgk(ngkmax, natmtot, nspnfv, nkpt))
do ik=1, nkpt
  do ispn=1, nspnfv
    if (isspinspiral()) then
! spin-spiral case
      if (ispn.eq.1) then
	vl(:)=vkl(:, ik)+0.5d0*input%groundstate%spin%vqlss(:)
	vc(:)=vkc(:, ik)+0.5d0*vqcss(:)
      else
	vl(:)=vkl(:, ik)-0.5d0*input%groundstate%spin%vqlss(:)
	vc(:)=vkc(:, ik)-0.5d0*vqcss(:)
      end if
    else
      vl(:)=vkl(:, ik)
      vc(:)=vkc(:, ik)
    end if
! generate the G+k-vectors
    call gengpvec(vl, vc, ngk(ispn, ik), igkig(:, ispn, ik), vgkl(:, :, ispn, ik), &
     vgkc(:, :, ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, ik))
! generate structure factors for G+k-vectors
    call gensfacgp(ngk(ispn, ik), vgkc(:, :, ispn, ik), ngkmax, sfacgk(:, :, ispn, ik))
  end do
end do

#ifdef XS
if (init1norealloc) goto 10
#endif
!---------------------------------!
!     APWs and local-orbitals     !
!---------------------------------!
! allocate linearisation energy arrays
if (allocated(apwe)) deallocate(apwe)
allocate(apwe(maxapword, 0:input%groundstate%lmaxapw, natmtot))
if (allocated(lorbe)) deallocate(lorbe)
allocate(lorbe(maxlorbord, maxlorb, natmtot))
nlomax=0
lolmax=0
apwordmax=0
do is=1, nspecies
! find the maximum APW order
  do l1=0, input%groundstate%lmaxapw
    apwordmax=max(apwordmax, apword(l1, is))
  end do
! set the APW linearisation energies to the default
  do ia=1, natoms(is)
    ias=idxas(ia, is)
    do l1=0, input%groundstate%lmaxapw
      do io=1, apword(l1, is)
	apwe(io, l1, ias)=apwe0(io, l1, is)
      end do
    end do
  end do
! find the maximum number of local-orbitals
  nlomax=max(nlomax, nlorb(is))
! set the local-orbital linearisation energies to the default
  do ia=1, natoms(is)
    ias=idxas(ia, is)
    do ilo=1, nlorb(is)
      lolmax=max(lolmax, lorbl(ilo, is))
      do io=1, lorbord(ilo, is)
	lorbe(io, ilo, ias)=lorbe0(io, ilo, is)
      end do
    end do
  end do
end do
lolmmax=(lolmax+1)**2
! generate the local-orbital index
call genidxlo
! allocate radial function arrays
if (allocated(apwfr)) deallocate(apwfr)
allocate(apwfr(nrmtmax, 2, apwordmax, 0:input%groundstate%lmaxapw, natmtot))
if (allocated(apwdfr)) deallocate(apwdfr)
allocate(apwdfr(apwordmax, 0:input%groundstate%lmaxapw, natmtot))
if (allocated(lofr)) deallocate(lofr)
allocate(lofr(nrmtmax, 2, nlomax, natmtot))
#ifdef XS
10 continue
#endif

!------------------------------------!
!     secular equation variables     !
!------------------------------------!
! number of first-variational states
nstfv=int(chgval/2.d0)+input%groundstate%nempty+1
! overlap and Hamiltonian matrix sizes
if (allocated(nmat)) deallocate(nmat)
allocate(nmat(nspnfv, nkpt))
if (allocated(npmat)) deallocate(npmat)
allocate(npmat(nspnfv, nkpt))
nmatmax=0
do ik=1, nkpt
  do ispn=1, nspnfv
    nmat(ispn, ik)=ngk(ispn, ik)+nlotot
    nmatmax=max(nmatmax, nmat(ispn, ik))
! packed matrix sizes
    npmat(ispn, ik)=(nmat(ispn, ik)*(nmat(ispn, ik)+1))/2
! the number of first-variational states should not exceed the matrix size
    nstfv=min(nstfv, nmat(ispn, ik))
  end do
end do
! number of second-variational states
nstsv=nstfv*nspinor
#ifdef XS
if (init1norealloc) goto 20
#endif
! allocate second-variational arrays
if (allocated(evalsv)) deallocate(evalsv)
allocate(evalsv(nstsv, nkpt))
if (allocated(occsv)) deallocate(occsv)
allocate(occsv(nstsv, nkpt))
occsv(:, :)=0.d0
#ifdef XS
if (allocated(occsv0)) deallocate(occsv0)
allocate(occsv0(nstsv, nkpt))
occsv(:, :)=0.d0
if (allocated(isto0)) deallocate(isto0)
allocate(isto0(nkpt))
isto0(:)=0.d0
if (allocated(isto)) deallocate(isto)
allocate(isto(nkpt))
isto(:)=0.d0
if (allocated(istu0)) deallocate(istu0)
allocate(istu0(nkpt))
istu0(:)=0.d0
if (allocated(istu)) deallocate(istu)
allocate(istu(nkpt))
istu(:)=0.d0
#endif
! allocate overlap and Hamiltonian integral arrays
if (allocated(oalo)) deallocate(oalo)
allocate(oalo(apwordmax, nlomax, natmtot))
if (allocated(ololo)) deallocate(ololo)
allocate(ololo(nlomax, nlomax, natmtot))
if (allocated(haa)) deallocate(haa)
allocate(haa(apwordmax, 0:input%groundstate%lmaxmat, apwordmax, 0:input%groundstate%lmaxapw, lmmaxvr, natmtot))
if (allocated(hloa)) deallocate(hloa)
allocate(hloa(nlomax, apwordmax, 0:input%groundstate%lmaxmat, lmmaxvr, natmtot))
if (allocated(hlolo)) deallocate(hlolo)
allocate(hlolo(nlomax, nlomax, lmmaxvr, natmtot))
! allocate and generate complex Gaunt coefficient array
if (allocated(gntyry)) deallocate(gntyry)
allocate(gntyry(lmmaxmat, lmmaxvr, lmmaxapw))
do l1=0, input%groundstate%lmaxmat
  do m1=-l1, l1
    lm1=idxlm(l1, m1)
    do l2=0, input%groundstate%lmaxvr
      do m2=-l2, l2
	lm2=idxlm(l2, m2)
	do l3=0, input%groundstate%lmaxapw
	  do m3=-l3, l3
	    lm3=idxlm(l3, m3)
	    gntyry(lm1, lm2, lm3)=gauntyry(l1, l2, l3, m1, m2, m3)
	  end do
	end do
      end do
    end do
  end do
end do
#ifdef XS
20 continue
#endif

call timesec(ts1)
timeinit=timeinit+ts1-ts0

return
end subroutine
!EOC
