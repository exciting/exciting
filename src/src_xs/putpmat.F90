
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_putpmat
  implicit none
contains

  subroutine putpmat(ik,tarec,filnam,pm)
    use modmain
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
    integer :: un,recl,ikr
    logical :: tarect
#ifdef MPI
    integer :: iproc,tag,status(MPI_STATUS_SIZE)
#endif
    ! TODO: use "tarec"
    tarect=tarec
    ikr=ik
    inquire(iolength=recl) vkl(:,ik),nstsv,pm
    call getunit(un)
#ifdef MPI
    tag=77
    if (rank.ne.0) call mpi_send(pm,size(pm),MPI_DOUBLE_COMPLEX,0,tag, &
         MPI_COMM_WORLD,ierr)
    if (rank.eq.0) then
       do iproc=0,lastproc(ik,nkpt)
          ikr=firstofset(iproc,nkpt)-1+ik
          if (iproc.ne.0) then
             ! receive data from slaves
             call mpi_recv(pm,size(pm),MPI_DOUBLE_COMPLEX,iproc,tag, &
                  MPI_COMM_WORLD,status,ierr)
          end if
#endif
          ! only master is performing I/O
          open(unit=un,file=trim(filnam),form='unformatted',action='write', &
               access='direct',recl=recl)
          write(un,rec=ikr) vkl(:,ikr),nstsv,pm
          close(un)
#ifdef MPI
       end do
    end if
#endif
  end subroutine putpmat

end module m_putpmat
