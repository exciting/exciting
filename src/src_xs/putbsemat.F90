
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putbsemat(fname,zmat,ikkp,iknr,jknr,iq,iqr,n1,n2,n3,n4)
  use modmain
  use modmpi
  use m_getunit
  implicit none
  ! arguments
  character(*), intent(in) :: fname
  complex(8), intent(in) :: zmat(n1,n2,n3,n4)
  integer, intent(in) :: ikkp,iknr,jknr,iq,iqr,n1,n2,n3,n4
  ! local variables
  integer :: recl,un,iknrr,jknrr,ikkpr,iq_r,iqr_r
#ifdef MPI
  integer :: nkkp,iproc,tag,status(MPI_STATUS_SIZE)
#endif
  ikkpr=ikkp
  iknrr=iknr
  jknrr=jknr
  iq_r=0
  iqr_r=0
  call getunit(un)
  inquire(iolength=recl) ikkp,iknr,jknr,iq,iqr,n1,n2,n3,n4,zmat
#ifdef MPI
  tag=77
  nkkp=(nkptnr*(nkptnr+1))/2
  if (rank.ne.0) call mpi_send(zmat,size(zmat),MPI_DOUBLE_COMPLEX,0,tag, &
       MPI_COMM_WORLD,ierr)
  if (rank.eq.0) then
     do iproc=0,lastproc(ikkp,nkkp)
        ikkpr=firstofset(iproc,nkkp)-1+ikkp
        call kkpmap(ikkpr,nkptnr,iknrr,jknrr)
        
        if (iproc.ne.0) then
           ! receive data from slaves
           call mpi_recv(zmat,size(zmat),MPI_DOUBLE_COMPLEX,iproc,tag, &
                MPI_COMM_WORLD,status,ierr)
        end if
#endif
        ! only the master is performing file I/O
        open(unit=un,file=trim(fname),form='unformatted',action='write', &
             access='direct',recl=recl)
        write(un,rec=ikkpr) ikkpr,iknrr,jknrr,iq_r,iqr_r,n1,n2,n3,n4,zmat
        close(un)
#ifdef MPI
     end do
  end if
#endif
end subroutine putbsemat
