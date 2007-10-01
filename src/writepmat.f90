
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepmat
! !INTERFACE:
subroutine writepmat
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT.OUT}.
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik,recl
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: pmat(:,:,:)
! initialise universal variables
call init0
call init1
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
! allocate the momentum matrix elements array
allocate(pmat(3,nstsv,nstsv))
! read in the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! find the record length
inquire(iolength=recl) pmat
open(50,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 status='REPLACE',recl=recl)
do ik=1,nkpt
! get the eigenvectors from file
  call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
! find the matching coefficients
  call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
! calculate the momentum matrix elements
  call genpmat(ngk(ik,1),igkig(1,ik,1),vgkc(1,1,ik,1),apwalm,evecfv,evecsv,pmat)
! write the matrix elements to direct-access file
  write(50,rec=ik) pmat
end do
close(50)
write(*,*)
write(*,'("Info(writepmat):")')
write(*,'(" momentum matrix elements written to file PMAT.OUT")')
write(*,*)
deallocate(apwalm,evecfv,evecsv,pmat)
end subroutine
!EOC

