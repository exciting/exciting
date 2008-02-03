
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine pmatgather
  use modmain
  use modxs
  use modmpi
  use m_filedel
  use m_getpmat
  use m_putpmat
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'pmatgather'
  complex(8), allocatable :: pm(:,:,:)
  integer :: ik,iproc

  ! allocate matrix elements array
  allocate(pm(3,nstsv,nstsv))

  ! file extension for q-point
  do iproc=0,procs-1
     call genfilname(basename='PMAT_XS',procs=procs,rank=iproc,&
          filnam=fnpmat_t)
     kpari=firstofset(iproc,nkpt)
     kparf=lastofset(iproc,nkpt)
     do ik=kpari,kparf
        ! momentum matrix elements
        call getpmat(ik,vkl,.false.,trim(fnpmat_t),pm)
        call putpmat(ik,.true.,trim(fnpmat),pm)
     end do
  end do

  ! delete partial files
  do iproc=0,procs-1
     call genfilname(basename='PMAT_XS',procs=procs,rank=iproc,&
          filnam=fnpmat_t)
     call filedel(trim(fnpmat_t))
  end do

  write(unitout,'(a,i8)') 'Info('//thisnam//'): Momentum matrix &
       &elements gathered.'

  ! deallocate
  deallocate(pm)

end subroutine pmatgather
