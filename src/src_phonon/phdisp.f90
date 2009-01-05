
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdisp
use modmain
implicit none
! local variables
integer iq,i,n,iv
real(8) wmin,wmax
! allocatable arrays
real(8), allocatable :: wp(:,:)
complex(8), allocatable :: ev(:,:)
complex(8), allocatable :: dynq(:,:,:)
complex(8), allocatable :: dynp(:,:)
complex(8), allocatable :: dynr(:,:,:)
! initialise universal variables
call init0
call init2
n=3*natmtot
allocate(wp(n,npp1d))
allocate(ev(n,n))
allocate(dynq(n,n,nqpt))
allocate(dynp(n,n))
allocate(dynr(n,n,ngridq(1)*ngridq(2)*ngridq(3)))
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! Fourier transform the dynamical matrices to real-space
call dynqtor(dynq,dynr)
! generate a set of q-point vectors along a path in the Brillouin zone
call connect(bvec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
wmin=0.d0
wmax=0.d0
! compute the phonon frequencies along the path
do iq=1,npp1d
! compute the dynamical matrix at this particular q-point
  call dynrtoq(vplp1d(:,iq),dynr,dynp)
! find the phonon frequencies and eigenvectors
  call dyndiag(dynp,wp(:,iq),ev)
  wmin=min(wmin,wp(1,iq))
  wmax=max(wmax,wp(n,iq))
end do
wmax=wmax+(wmax-wmin)*0.5d0
wmin=wmin-(wmax-wmin)*0.5d0
! output the vertex location lines
open(50,file='PHDLINES.OUT',action='WRITE',form='FORMATTED')
do iv=1,nvp1d
  write(50,'(2G18.10)') dvp1d(iv),wmin
  write(50,'(2G18.10)') dvp1d(iv),wmax
  write(50,'("     ")')
end do
close(50)
! output the phonon dispersion
open(50,file='PHDISP.OUT',action='WRITE',form='FORMATTED')
do i=1,n
  do iq=1,npp1d
    write(50,'(2G18.10)') dpp1d(iq),wp(i,iq)
  end do
  write(50,'("     ")')
end do
close(50)
write(*,*)
write(*,'("Info(phdisp):")')
write(*,'(" phonon dispersion written to PHDISP.OUT")')
write(*,'(" vertex location lines written to PHDLINES.OUT")')
write(*,*)
deallocate(wp,ev,dynq,dynp,dynr)
return
end subroutine

