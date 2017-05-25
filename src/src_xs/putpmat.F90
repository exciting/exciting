! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module m_putpmat

  implicit none

  contains

    !BOP
    ! !ROUTINE: putpmat
    ! !INTERFACE:
    subroutine putpmat(ik, filnam, pm, tag)
    ! !USES:
      use modmain
      use modmpi
      use m_getunit
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: ik
    ! character(*) :: filnam
    ! integer(4) :: tag
    ! IN/OUT:
    ! complex(8) :: pm(:,:,:) 
    !
    ! !DESCRIPTION:
    !   The routine collects the momentum matrix elements
    !   for each k-point form the {\tt MPI} processes and
    !   writes them to a direct access file.
    !   Note: The content of pm is destroyed on exit.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. 2016 (Aurich)
    !   Removed parts which literately did nothing. (Aurich)
    !EOP
    !BOC

      implicit none

      ! arguments
      integer, intent(in) :: ik
      character(*), intent(in) :: filnam
      complex(8), intent(inout) :: pm(:, :, :)
      integer, intent(in), optional :: tag

      integer :: un, reclen, ikr

#ifdef MPI
      integer :: iproc, mpitag, stat(mpi_status_size)
#endif

      ikr = ik
      inquire(iolength=reclen) vkl(:, ik), nstsv, pm
      call getunit(un)

#ifdef MPI
      if(present(tag)) then
        mpitag = tag
      else
        mpitag = 77
      end if

      if(rank .ne. 0) then 
        call mpi_send(pm, size(pm),&
          & mpi_double_complex, 0, mpitag, mpi_comm_world, ierr)
      end if

      if(rank .eq. 0) then

        ! For each call form rank 0 there are lastproc(ik, nkpt) 
        ! sends from the other ranks to rank 0.
        do iproc = 0, lastproc(ik, nkpt)

          ! Calculate ik form sender
          ikr = firstofset(iproc, nkpt) - 1 + ik

          if(iproc .ne. 0) then
            ! receive data from slaves
            call mpi_recv(pm, size(pm), mpi_double_complex,&
              & iproc, mpitag, mpi_comm_world, stat, ierr)
          end if
#endif
          ! only master is performing i/o
          open(unit=un, file=trim(filnam), form='unformatted',&
            & action='write', access='direct', recl=reclen)
          write(un, rec=ikr) vkl(:, ikr), nstsv, pm
          close(un) 

#ifdef MPI
        end do

      end if
#endif

    end subroutine putpmat
    !EOC

end module m_putpmat
