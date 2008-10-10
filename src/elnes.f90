
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine elnes
use modmain
implicit none
! local variables
integer ik,ist,jst
integer n,nsk(3),iw
real(8) wd,dw,w,t1
real(8) vecqc(3),qc
! allocatable arrays
real(8), allocatable :: e(:,:,:)
real(8), allocatable :: f(:,:,:)
real(8), allocatable :: eps2(:)
complex(8), allocatable :: emat(:,:)
! initialise universal variables
call init0
call init1
! allocate local arrays
allocate(e(nstsv,nstsv,nkpt))
allocate(f(nstsv,nstsv,nkpt))
allocate(eps2(nwdos))
! allocate the matrix elements array for < i,k+G+q | exp(iq.r) | j,k >
allocate(emat(nstsv,nstsv))
! read in the density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! loop over k-points
do ik=1,nkpt
! get the second-variational eigenvalues and occupancies from file
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
! compute < i,k+G+q | exp(iq.r) | j,k > matrix elements
  call genexpiqr(ik,emat)
  do ist=1,nstsv
    do jst=1,nstsv
      e(ist,jst,ik)=evalsv(jst,ik)-evalsv(ist,ik)
      t1=dble(emat(ist,jst))**2+aimag(emat(ist,jst))**2
      f(ist,jst,ik)=t1*occsv(ist,ik)*(occmax-occsv(jst,ik))
    end do
  end do
end do
! number of subdivisions used for interpolation
nsk(:)=max(ngrdos/ngridk(:),1)
n=nstsv*nstsv
! integrate over the Brillouin zone
call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,n,n,e,f,eps2)
! q-vector in Cartesian coordinates
call r3mv(bvec,vecql,vecqc)
qc=sqrt(vecqc(1)**2+vecqc(2)**2+vecqc(3)**2)
t1=occmax/omega
if (qc.gt.epslat) t1=t1/qc**2
eps2(:)=t1*eps2(:)
open(50,file='ELNES.OUT',action='WRITE',form='FORMATTED')
wd=wdos(2)-wdos(1)
dw=wd/dble(nwdos)
do iw=1,nwdos
  w=dw*dble(iw-1)+wdos(1)
  write(50,'(2G18.10)') w,eps2(iw)
end do
close(50)
write(*,*)
write(*,'("Info(elnes):")')
write(*,'(" ELNES intensity distribution written to ELNES.OUT")')
write(*,*)
deallocate(e,f,eps2,emat)
return
end subroutine

