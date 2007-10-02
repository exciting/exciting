
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: seceqnfv
! !INTERFACE:
subroutine seceqnfv(ik,ispn,apwalm,evalfv,evecfv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point number (in,integer)
!   ispn   : first-variational spin index (in,integer)
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
integer, intent(in) :: ik
integer, intent(in) :: ispn
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(out) :: evalfv(nstfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv)
! local variables
integer is,ia,i,m,n,np,info
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
if ((ik.lt.1).or.(ik.gt.nkpt)) then
  write(*,*)
  write(*,'("Error(seceqnfv): k-point out of range : ",I8)') ik
  write(*,*)
  stop
end if
n=nmat(ik,ispn)
np=npmat(ik,ispn)
allocate(iwork(5*n))
allocate(ifail(n))
allocate(w(n))
allocate(rwork(7*n))
allocate(v(1))
allocate(h(np))
allocate(o(np))
allocate(work(2*n))
!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!
call hamiltonandoverlapsetup(np,ik,ispn,apwalm,h,o)
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
call zhpgvx(1,'V','I','U',n,h,o,vl,vu,1,nstfv,abstol,m,w,evecfv,nmatmax,work, &
 rwork,iwork,ifail,info)
evalfv(1:nstfv)=w(1:nstfv)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(seceqnfv): diagonalisation failed for k-point ",I8)') ik
  write(*,'(" ZHPGVX returned INFO = ",I8)') info
  if (info.gt.n) then
    i=info-n
    write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
    write(*,'("  is not positive definite")')
    write(*,'(" Order of overlap matrix : ",I8)') n
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

