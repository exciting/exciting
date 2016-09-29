! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module m_putpmat

  implicit none

  contains

    subroutine putpmat(ik, tarec, filnam, pm, tag)
      use modmain
      use modmpi
      use m_getunit

      implicit none

      ! arguments
      integer, intent(in) :: ik
      integer, intent(in), optional :: tag

      ! true if absolut record position is ik
      logical, intent(in) :: tarec
      character(*), intent(in) :: filnam
      complex(8), intent(inout) :: pm(:, :, :)

      ! local variables
      integer :: un, reclen, ikr
      logical :: tarect

#ifdef MPI
      integer :: iproc, mpitag, stat(mpi_status_size)
#endif

      tarect = tarec
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
        call mpi_send(pm, size(pm), mpi_double_complex, 0, mpitag, mpi_comm_world, ierr)
      end if

      if(rank .eq. 0) then

        do iproc = 0, lastproc(ik, nkpt)

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

end module m_putpmat
