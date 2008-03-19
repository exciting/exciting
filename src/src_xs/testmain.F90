
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine testmain
  use modmain
  use m_ftfun
  implicit none




  complex(8), allocatable :: fmt(:,:,:)
  complex(8), allocatable :: fir(:)
  complex(8),allocatable :: gft(:)
  integer :: ir

  call init0

  allocate(fmt(lmmaxvr,nrmtmax,natmtot))
  allocate(fir(ngrtot),gft(ngvec))

  fmt(:,:,:)=0.d0
  fmt(1,:,:)=1.d0/y00
  fir(:)=1.d0
  call ftfun(ngvec,tir=.true.,tmt=.true.,gir=fir,gmt=fmt,ftg=gft)

  do ir=1,ngvec
     write(2000,'(3g18.10)') gft(ir), abs(gft(ir))
  end do

  write(*,*) 'normalization:',sum(abs(gft)**2)

  deallocate(fir,fmt,gft)


end subroutine

