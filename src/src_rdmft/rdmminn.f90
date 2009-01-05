
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rdmminn
! minimise the total energy w.r.t. occupation numbers
use modmain
implicit none
! allocatable arrays
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
integer ik,it,idm
real(8) ep,de
! parameter to check energy convergence
real(8), parameter :: eps=1.d-8
! allocate arrays
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
open(61,file='RDMN_ENERGY.OUT',action='WRITE',form='FORMATTED')
if (spinpol) then
  open(62,file='RDMN_MOMENT.OUT',action='WRITE',form='FORMATTED')
end if
! calculate the non-local matrix elements (i-jj-i)
if ((rdmxctype.ne.0).and.(maxitc.lt.1)) then
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ik=1,nkpt
!$OMP CRITICAL
    write(*,'("Info(rdmminn): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
    call rdmvnln(ik)
  end do
!$OMP END DO
!$OMP END PARALLEL
end if
ep=0.d0
! begin iteration loop
do it=1,maxitn
  write(*,'("Info(rdmminn): iteration ",I4," of ",I4)') it,maxitn
! vary the occupation numbers
  call rdmvaryn
! zero the density
  rhomt(:,:,:)=0.d0
  rhoir(:)=0.d0
! zero the magnetisation
  if (spinpol) then
    magmt(:,:,:,:)=0.d0
    magir(:,:)=0.d0
  end if
! compute the charge density and magnetisation with the new occupancies
  do ik=1,nkpt
! get the eigenvectors from file
    call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
    call getevecsv(vkl(:,ik),evecsv)
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
      call rfmtctof(magmt(:,:,:,idm))
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
! calculate Coulomb potential matrix elements (RDM states)
  call genvmat(vclmt,vclir,vclmat)
! calculate the energy
  call rdmenergy
! check for convergence
  de=ep-engytot
  if (it.gt.1) then
    if (abs(de).lt.eps) goto 10
  end if
  ep=engytot
! write energy and convergence factor to a file
  write(61,'(I6,2G18.10)') it,engytot,de
  call flushifc(61)
! end iteration loop
end do
10 continue
close(61)
if (spinpol) close(62)
deallocate(evecfv,evecsv)
return
end subroutine

