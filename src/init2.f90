
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init2
use modmain
implicit none
! local variables
integer is,ia,ist,ic,m
real(8) cpu0,cpu1
real(8) vqloff(3)

call cpu_time(cpu0)

!---------------------!
!     q-point set     !
!---------------------!
! check if the system is an isolated molecule
if (molecule) ngridq(:)=1
! OEP, Hartree-Fock or RDMFT
if ((xctype.lt.0).or.(task.eq.5).or.(task.eq.300)) then
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
! the q-point offset should always be zero
vqloff(:)=0.d0
! generate the q-point set
call genppts(reduceq,ngridq,vqloff,nqpt,iqmap,ivq,vql,vqc,wqpt)

!-----------------------------------------------!
!     OEP, Hartree-Fock and RDMFT variables     !
!-----------------------------------------------!
if ((xctype.lt.0).or.(task.eq.5).or.(task.eq.300)) then
! determine the 1/q^2 integral weights if required
  call genwiq2
! output the 1/q^2 integrals to WIQ2.OUT
  call writewiq2
end if
if (xctype.lt.0) then
! initialise OEP residue magnitude
  resoep=1.d0
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
! maximum core states over all species
      ncrmax=max(ncrmax,ic)
    end do
  end do
end if
if ((task.eq.5).or.(task.eq.300)) then
! allocate the kinetic matrix elements for Hartree-Fock/RDMFT
  if (allocated(kinmatc)) deallocate(kinmatc)
  allocate(kinmatc(nstsv,nstsv,nkpt))
end if
if (task.eq.300) then
  if (allocated(vclmat)) deallocate(vclmat)
  allocate(vclmat(nstsv,nstsv,nkpt))
  if (allocated(dkdc)) deallocate(dkdc)
  allocate(dkdc(nstsv,nstsv,nkpt))
  if (allocated(vnlmatr)) deallocate(vnlmatr)
  allocate(vnlmatr(nstsv,nkpt,nstsv,nkptnr))
  if (allocated(vnlmat)) deallocate(vnlmat)
  allocate(vnlmat(nstsv,nstsv,nkpt,nstsv,nkptnr))
end if

call cpu_time(cpu1)
timeinit=timeinit+cpu1-cpu0

return
end subroutine

