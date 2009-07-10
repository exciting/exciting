
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dfgather
  use modmain
  use modxs
  use modmpi
  use m_filedel
  use m_getx0
  use m_putx0
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='dfgather'
  integer :: n,iq,iw,iproc
  complex(8), allocatable :: chi0(:,:),chi0wg(:,:,:),chi0hd(:,:)
  logical :: tq0
  logical, external :: tqgamma
  ! loop over q-points
  do iq=1,nqpt
     tq0=tqgamma(iq)
     ! calculate k+q and G+k+q related variables
     call init1offs(qvkloff(1,iq))
     ! size of local field effects
     n=ngq(iq)
     ! allocate
     allocate(chi0(n,n),chi0wg(n,2,3),chi0hd(3,3))
     ! file extension for q-point
     do iproc=0,procs-1
        call genfilname(basename='X0',bzsampl=bzsampl,acont=acont,&
             nar=.not.aresdf,tord=torddf,markfxcbse=tfxcbse,iqmt=iq,procs=procs,rank=iproc,filnam=fnchi0_t)
        wpari=firstofset(iproc,nwdf)
        wparf=lastofset(iproc,nwdf)
        do iw=wpari,wparf
           ! exponential factor matrix elements
           call getx0(tq0,iq,iw-wpari+1,trim(fnchi0_t),'',chi0,&
                chi0wg,chi0hd)
           call putx0(tq0,iq,iw,trim(fnchi0),'',chi0,chi0wg,chi0hd)
        end do
     end do
     do iproc=0,procs-1
        call genfilname(basename='X0',iqmt=iq,procs=procs,rank=iproc,&
             filnam=fnchi0_t)
     end do
     deallocate(chi0,chi0wg,chi0hd)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Kohn Sham response &
          &function gathered for q-point:',iq
  end do
end subroutine dfgather
