
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init2
use modmain
implicit none
! local variables
integer is,ia,ist,ic,m
real(8) ts0,ts1
real(8) boxl(3,4)

call timesec(ts0)

!---------------------!
!     q-point set     !
!---------------------!
! check if the system is an isolated molecule
if (molecule) ngridq(:)=1
! OEP, Hartree-Fock or RDMFT
if ((xctype.lt.0).or.(task.eq.5).or.(task.eq.6).or.(task.eq.300)) then
  ngridq(:)=ngridk(:)
  reduceq=.false.
end if
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
call genppts(reduceq,.true.,ngridq,boxl,nqpt,iqmap,ivq,vql,vqc,wqpt)

!-----------------------------------------------!
!     OEP, Hartree-Fock and RDMFT variables     !
!-----------------------------------------------!
if ((xctype.lt.0).or.(task.eq.5).or.(task.eq.6).or.(task.eq.300)) then
! determine the 1/q^2 integral weights if required
  call genwiq2
! output the 1/q^2 integrals to WIQ2.OUT
  call writewiq2
end if
if (xctype.lt.0) then
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
  if (spinpol) then
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

