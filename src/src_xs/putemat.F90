
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_putemat
  implicit none
contains

  subroutine putemat(iq,ik,tarec,filnam,x1,x2)
    use modmain
    use modmpi
    use modxs
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    logical :: tarec
    character(*), intent(in) :: filnam
    complex(8), intent(in) :: x1(:,:,:)
    complex(8), optional, intent(in) :: x2(:,:,:)
    ! local variables
    integer :: un,recl,ikr
    logical :: tarec_
#ifdef MPI
    integer :: iproc,tag1,tag2,status(MPI_STATUS_SIZE)
#endif
    !TODO: use "tarec"; access subset of bands!!! UNFINISHED
    tarec_=tarec
    ikr=ik
    call getunit(un)
    if (present(x2)) then
#ifdef MPI
       tag1=77
       tag2=78
       if (rank.ne.0) call mpi_send(x1,size(x1),MPI_DOUBLE_COMPLEX,0,tag1, &
            MPI_COMM_WORLD,ierr)
       if (rank.ne.0) call mpi_send(x2,size(x2),MPI_DOUBLE_COMPLEX,0,tag2, &
            MPI_COMM_WORLD,ierr)
       if (rank.eq.0) then
          do iproc=0,lastproc(ik,nkpt)
             ikr=firstofset(iproc,nkpt)-1+ik
             if (iproc.ne.0) then
                ! receive data from slaves
                call mpi_recv(x1,size(x1),MPI_DOUBLE_COMPLEX,iproc,tag1, &
                     MPI_COMM_WORLD,status,ierr)
                call mpi_recv(x2,size(x2),MPI_DOUBLE_COMPLEX,iproc,tag2, &
                     MPI_COMM_WORLD,status,ierr)
             end if
#endif
             ! I/O record length
             inquire(iolength=recl) vql(:,iq),vkl(:,ikr),nstsv,ngq(iq), &
                  nst1,nst2,nst3,nst4,x1,x2
             open(unit=un,file=trim(filnam),form='unformatted', &
                  action='write',access='direct',recl=recl)
             write(un,rec=ikr) vql(:,iq),vkl(:,ikr),nstsv,ngq(iq), &
                  nst1,nst2,nst3,nst4,x1,x2
#ifdef MPI
          end do
       end if
#endif
    else
#ifdef MPI
       tag1=77
       if (rank.ne.0) call mpi_send(x1,size(x1),MPI_DOUBLE_COMPLEX,0,tag1, &
            MPI_COMM_WORLD,ierr)
       if (rank.eq.0) then
          do iproc=0,lastproc(ik,nkpt)
             ikr=firstofset(iproc,nkpt)-1+ik
             if (iproc.ne.0) then
                ! receive data from slaves
                call mpi_recv(x1,size(x1),MPI_DOUBLE_COMPLEX,iproc,tag1, &
                     MPI_COMM_WORLD,status,ierr)
             end if
#endif
             ! I/O record length
             inquire(iolength=recl) vql(:,iq),vkl(:,ikr),nstsv,ngq(iq), &
                  nst1,nst2,x1
             open(unit=un,file=trim(filnam),form='unformatted', &
                  action='write',access='direct',recl=recl)
             write(un,rec=ikr) vql(:,iq),vkl(:,ikr),nstsv,ngq(iq), &
                  nst1,nst2,x1
#ifdef MPI
          end do
       end if
#endif
    end if
    close(un)
  end subroutine putemat

end module m_putemat
