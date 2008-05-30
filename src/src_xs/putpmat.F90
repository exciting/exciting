
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_putpmat
  implicit none
contains

  subroutine putpmat(ik,tarec,filnam,pm)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: ik
    ! true if absolut record position is ik
    logical, intent(in) :: tarec
    character(*), intent(in) :: filnam
    complex(8), intent(in) :: pm(:,:,:)
    ! local variables
    integer :: un, recl, ikr
#ifdef MPI
    integer :: n,n0,iproc,tag,dest,source,status(MPI_STATUS_SIZE)
#endif
    ! record position for k-point
    ikr=ik
    ! record position is not absolute k-point index
    if (.not.tarec) call getridx(procs,nkpt,ik,ikr)
    ! I/O record length
    inquire(iolength=recl) nstsv,nkpt,vkl(:,ik),pm
    call getunit(un)
#ifdef MPI
    tag=77
    if (rank.ne.0) then
       dest=0
       call mpi_send(pm,size(pm),MPI_DOUBLE_COMPLEX,dest,tag, &
            MPI_COMM_WORLD,ierr)
    end if
    if (rank.eq.0) then
       n0=nofset(0,nkpt)
       do iproc=0,procs-1
          ikr=firstofset(iproc,nkpt)-1+ik
          n=nofset(iproc,nkpt)
          if ((ik.eq.n0).and.(n.lt.n0)) cycle
          if (iproc.ne.0) then
             ! receive data from slaves
             source=iproc
             call mpi_recv(pm,size(pm),MPI_DOUBLE_COMPLEX,source,tag, &
                  MPI_COMM_WORLD,status,ierr)
          end if
#endif
          ! only master is performing I/O
          open(unit=un,file=trim(filnam),form='unformatted',action='write', &
               access='direct',recl=recl)
          write(un,rec=ikr) nstsv,nkpt,vkl(:,ikr),pm
          close(un)
#ifdef MPI
       end do
    end if
#endif

  end subroutine putpmat

end module m_putpmat
