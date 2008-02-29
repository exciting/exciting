
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rdmminc
! minimises the total energy w.r.t. evecsv using steepest descent
use modmain
implicit none
integer it,ik,idm
real(8) sum,sp,ds
! parameter to check energy convergence
real(8), parameter :: eps=1.d-10
! allocatable arrays
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
! allocate arrays
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
sp=0.d0
ds=0.d0
open(61,file='RDMC_ENERGY.OUT',action='WRITE',form='FORMATTED')
if (spinpol) then
  open(62,file='RDMC_MOMENT.OUT',action='WRITE',form='FORMATTED')
end if
write(*,*)
! begin iteration loop
do it=1,maxitc
  write(*,'("Info(rdmminc): iteration ",I4," of ",I4)') it,maxitc
! vary evecsv and orthogonalise it
  call rdmvaryc(sum)
! zero the density
  rhomt(:,:,:)=0.d0
  rhoir(:)=0.d0
! zero the magnetisation
  if (spinpol) then
    magmt(:,:,:,:)=0.d0
    magir(:,:)=0.d0
  end if
  do ik=1,nkpt
! get the eigenvectors and values from file
    call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
    call getevecsv(vkl(1,ik),evecsv)
! calculate the density
    call rhovalk(ik,evecfv,evecsv)
  end do
! symmetrise the density
  call symrf(lradstp,rhomt,rhoir)
! convert the muffin-tin density from coarse to a fine grid
  call rfmtctof(rhomt)
  if (spinpol) then
! symmetrise the magnetisation
    call symrvf(lradstp,magmt,magir)
! convert the magnetisation from a coarse to a fine radial mesh
    do idm=1,ndmag
      call rfmtctof(magmt(1,1,1,idm))
    end do
  end if
! add core density to the valence density
  call addrhocr
! calculate the charges
  call charge
! calculate the magnetic moment
  if (spinpol) then
    call moment
    write(62,'(I6,3G18.10)') it,momtot(1:ndmag)
    call flushifc(62)
  end if
! normalise the density
  call rhonorm
! calculate the Coulomb potential
  call potcoul
! calculate Coulomb matrix elements
  call genvmat(vclmt,vclir,vclmat)
! calculate derivative of kinetic energy w.r.t. evecsv
  call rdmdkdc
! calculate the energy
  call rdmenergy
! check for convergence of derivative of energy w.r.t. evecsv
  if (it.gt.1) then
    ds=sp-sum
    sp=sum
  end if
! write energy and convergence factor to a file
  write(61,'(I6,2G18.10)') it,engytot,ds
  call flushifc(61)
! end iteration loop
end do
close(61)
close(62)
deallocate(evecfv,evecsv)
return
end subroutine

