! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: putbsemat
! !INTERFACE:
subroutine putbsemat(fname, tag, zmat, ikkp, iknr, jknr, iq, iqr, n1, n2, n3, n4)
! !USES:
  use mod_kpoint, only: nkptnr
  use modmpi
  use m_getunit
! !DESCRIPTION:
!   The routine writes complex 4d-array to a direct access file and
!   is intended for the use in be BSE part of the code. It is used
!   to write the screened coulomb interaction {\tt SCCLI.OUT} and 
!   the exchange interaction {\tt EXCLI.OUT} to file.
!
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC

  implicit none

  ! Arguments
  character(*), intent(in) :: fname
  complex(8), intent(in) :: zmat(n1, n2, n3, n4)
  integer, intent(in) :: ikkp, iknr, jknr, iq, iqr, n1, n2, n3, n4
  integer, intent(in) :: tag

  ! Local variables
  integer :: reclen, un, iknrr, jknrr, ikkpr, iq_r, iqr_r
#ifdef MPI
  integer :: nkkp, iproc, stat(mpi_status_size)
#endif

  ikkpr = ikkp
  iknrr = iknr
  jknrr = jknr
  iq_r = 0
  iqr_r = 0

  call getunit(un)
  inquire(iolength=reclen) ikkp, iknr, jknr, iq, iqr, n1, n2, n3, n4, zmat

#ifdef MPI
  nkkp = (nkptnr*(nkptnr+1)) / 2

  if(rank .ne. 0) then
    call mpi_send(zmat, size(zmat), mpi_double_complex, 0, tag, mpi_comm_world, ierr)
  end if

  if(rank .eq. 0) then
    do iproc = 0, lastproc(ikkp, nkkp)
      ikkpr = firstofset(iproc, nkkp) - 1 + ikkp
      call kkpmap(ikkpr, nkptnr, iknrr, jknrr)

      if(iproc .ne. 0) then
        ! Receive data from slaves
        call mpi_recv(zmat, size(zmat), mpi_double_complex, iproc, tag,&
          & mpi_comm_world, stat, ierr)
      end if
#endif
      ! Only the master is performing file i/o
      open(unit=un, file=trim(fname), form='unformatted', action='write',&
        & access='direct', recl=reclen)
      write(un, rec=ikkpr) ikkpr, iknrr, jknrr, iq_r, iqr_r, n1, n2, n3, n4, zmat
      close(un)
#ifdef MPI
    end do
  end if
#endif
end subroutine putbsemat
!EOC
