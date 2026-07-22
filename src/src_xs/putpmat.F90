! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module m_putpmat
  use mod_large_io, only: inquire_large, open_direct_unformatted_large
  use precision, only: i32, long_int, dp

  implicit none

  contains

    !BOP
    ! !ROUTINE: putpmat
    ! !INTERFACE:
    subroutine putpmat(ik, filnam, pm, tag)
    ! !USES:
      use modmain
      use modmpi
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
      integer(i32), intent(in) :: ik
      character(*), intent(in) :: filnam
      complex(dp), intent(inout) :: pm(:, :, :)
      integer(i32), intent(in), optional :: tag

      integer(i32) :: un, ikr
      integer(long_int) :: reclen

#ifdef MPI
      integer(i32) :: iproc, mpitag, stat(mpi_status_size)
#endif

      ikr = ik
      call inquire_large( reclen, vkl(:, ik), [nstsv], pm )

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
          call open_direct_unformatted_large( un, trim( filnam ), "write", reclen, "unknown" )
          write(un, rec=ikr) vkl(:, ikr), nstsv, pm
          close(un) 

#ifdef MPI
        end do

      end if
#endif

    end subroutine putpmat
    !EOC

end module m_putpmat
