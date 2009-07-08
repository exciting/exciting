

! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_putemat
  implicit none
contains


subroutine putemat(iq, ik, tarec, filnam, l1, h1, l2, h2, x12, l3, h3, l4, h4, x34)
    use modmain
    use modmpi
    use modxs
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq, ik
    logical :: tarec
    character(*), intent(in) :: filnam
    integer, intent(in) :: l1, h1, l2, h2
    complex(8), intent(in) :: x12(:, :, :)
    integer, optional, intent(in) :: l3, h3, l4, h4
    complex(8), optional, intent(in) :: x34(:, :, :)
    ! local variables
    integer :: un, recl, ikr
    logical :: tarec_
#ifdef MPI
    integer :: iproc, tag1, tag2, status(MPI_STATUS_SIZE)
#endif
    ! check if all optional variables are present if any is present
    if ((present(l3).or.present(h3).or.present(l4).or.present(h4).or. &
	 present(x34)).and.(.not.( &
	 present(l3).and.present(h3).and.present(l4).and.present(h4).and. &
	 present(x34)))) then
       write(*, *)
       write(*, '("Error(putemat): optional parameters not complete - check &
	    &routine")')
       write(*, *)
       call terminate
    end if
    !TODO: use "tarec"
    tarec_=tarec
    ikr=ik
    call getunit(un)
    if (present(x34)) then
#ifdef MPI
       tag1=77
       tag2=78
       if (rank.ne.0) call mpi_send(x12, size(x12), MPI_DOUBLE_COMPLEX, 0, tag1, &
	    MPI_COMM_WORLD, ierr)
       if (rank.ne.0) call mpi_send(x34, size(x34), MPI_DOUBLE_COMPLEX, 0, tag2, &
	    MPI_COMM_WORLD, ierr)
       if (rank.eq.0) then
	  do iproc=0, lastproc(ik, nkpt)
	     ikr=firstofset(iproc, nkpt)-1+ik
	     if (iproc.ne.0) then
                ! receive data from slaves
		call mpi_recv(x12, size(x12), MPI_DOUBLE_COMPLEX, iproc, tag1, &
		     MPI_COMM_WORLD, status, ierr)
		call mpi_recv(x34, size(x34), MPI_DOUBLE_COMPLEX, iproc, tag2, &
		     MPI_COMM_WORLD, status, ierr)
	     end if
#endif
             ! I/O record length
	     inquire(iolength = recl) vql(:, iq), vkl(:, ikr), nstsv, ngq(iq), &
		  l1, h1, l2, h2, l3, h3, l4, h4, x12, x34
	     open(unit = un, file = trim(filnam), form = 'unformatted', &
		  action = 'write', access = 'direct', recl = recl)
	     write(un, rec = ikr) vql(:, iq), vkl(:, ikr), nstsv, ngq(iq), &
		  l1, h1, l2, h2, l3, h3, l4, h4, x12, x34
#ifdef MPI
	  end do
       end if
#endif
    else
#ifdef MPI
       tag1=77
       if (rank.ne.0) call mpi_send(x12, size(x12), MPI_DOUBLE_COMPLEX, 0, tag1, &
	    MPI_COMM_WORLD, ierr)
       if (rank.eq.0) then
	  do iproc=0, lastproc(ik, nkpt)
	     ikr=firstofset(iproc, nkpt)-1+ik
	     if (iproc.ne.0) then
                ! receive data from slaves
		call mpi_recv(x12, size(x12), MPI_DOUBLE_COMPLEX, iproc, tag1, &
		     MPI_COMM_WORLD, status, ierr)
	     end if
#endif
             ! I/O record length
	     inquire(iolength = recl) vql(:, iq), vkl(:, ikr), nstsv, ngq(iq), &
		  l1, h1, l2, h2, x12
	     open(unit = un, file = trim(filnam), form = 'unformatted', &
		  action = 'write', access = 'direct', recl = recl)
	     write(un, rec = ikr) vql(:, iq), vkl(:, ikr), nstsv, ngq(iq), &
		  l1, h1, l2, h2, x12
#ifdef MPI
	  end do
       end if
#endif
    end if
    close(un)
  end subroutine putemat

end module m_putemat
