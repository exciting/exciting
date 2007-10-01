
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: seceqnfv
! !INTERFACE:
subroutine seceqnfv(nmatp,ngp,igpig,vgpc,apwalm,evalfv,evecfv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nmatp  : order of overlap and Hamiltonian matrices (in,integer)
!   ngp    : number of G+k-vectors for augmented plane waves (in,integer)
!   igpig  : index from G+k-vectors to G-vectors (in,integer(ngkmax))
!   vgpc   : G+k-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evalfv : first-variational eigenvalues (out,real(nstfv))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
! !DESCRIPTION:
!   Solves the secular equation,
!   $$ (H-\epsilon O)\Phi=0, $$
!   for the all the first-variational states of the input $k$-point.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nmatp
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(out) :: evalfv(nstfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv)
! local variables
integer is,ia,i,m,np,info
real(8) vl,vu,abstol
real(8) cpu0,cpu1
! external functions
real(8) dlamch
external dlamch
! allocatable arrays
integer, allocatable :: iwork(:)
integer, allocatable :: ifail(:)
real(8), allocatable :: w(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: v(:)
complex(8), allocatable :: h(:)
complex(8), allocatable :: o(:)
complex(8), allocatable :: work(:)
np=(nmatp*(nmatp+1))/2
allocate(iwork(5*nmatp))
allocate(ifail(nmatp))
allocate(w(nmatp))
allocate(rwork(7*nmatp))
allocate(v(1))
allocate(h(np))
allocate(o(np))
allocate(work(2*nmatp))
!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!
call cpu_time(cpu0)
! set the matrices to zero
h(:)=0.d0
o(:)=0.d0
! muffin-tin contributions
do is=1,nspecies
  do ia=1,natoms(is)
    call hmlaa(.false.,is,ia,ngp,apwalm,v,h)
    call hmlalo(.false.,is,ia,ngp,apwalm,v,h)
    call hmllolo(.false.,is,ia,ngp,v,h)
    call olpaa(.false.,is,ia,ngp,apwalm,v,o)
    call olpalo(.false.,is,ia,ngp,apwalm,v,o)
    call olplolo(.false.,is,ia,ngp,v,o)
  end do
end do
! interstitial contributions
call hmlistl(.false.,ngp,igpig,vgpc,v,h)
call olpistl(.false.,ngp,igpig,v,o)
call cpu_time(cpu1)
!$OMP CRITICAL
timemat=timemat+cpu1-cpu0
!$OMP END CRITICAL
!------------------------------------!
!     solve the secular equation     !
!------------------------------------!
call cpu_time(cpu0)
vl=0.d0
vu=0.d0
abstol=2.d0*dlamch('S')
! LAPACK 3.0 call
call zhpgvx(1,'V','I','U',nmatp,h,o,vl,vu,1,nstfv,abstol,m,w,evecfv,nmatmax, &
 work,rwork,iwork,ifail,info)
evalfv(1:nstfv)=w(1:nstfv)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(seceqnfv): diagonalisation failed")')
  write(*,'(" ZHPGVX returned INFO = ",I8)') info
  if (info.gt.nmatp) then
    i=info-nmatp
    write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
    write(*,'("  is not positive definite")')
    write(*,'(" Order of overlap matrix : ",I8)') nmatp
    write(*,*)
  end if
  stop
end if
call cpu_time(cpu1)
!$OMP CRITICAL
timefv=timefv+cpu1-cpu0
!$OMP END CRITICAL
deallocate(iwork,ifail,w,rwork,v,h,o,work)
return
end subroutine
!EOC

