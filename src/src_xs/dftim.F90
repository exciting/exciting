
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dftim

contains

  subroutine dftim(iq,ik,filnam,t1,t2,t3,t4,leta)
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    integer, intent(in),optional :: leta(4)
    character(*), intent(in) :: filnam
    real(8), intent(in) :: t1,t2,t3,t4
    ! local variables
    integer :: un
    call getunit(un)
    open(un,file=trim(filnam),action='write',form='formatted', &
         position='append')
    write(un,*)
    write(un,'("Timings (CPU seconds) for q-point/k-point: ",2i6)')iq,ik
    write(un,'("  reading matrix elements          : ",f14.4)') t1
    write(un,'("  generating oscillators           : ",f14.4)') t2
    write(un,'("  updating response function       : ",f14.4)') t3
    write(un,'("  total                            : ",f14.4)') t4
    if (present(leta)) then
       write(un,'(a,i6,a,i6,a,i6,a,i6,a)') '  ETA :', &
            leta(1), ' d ', &
            leta(2), ' h ', &
            leta(3), ' m ', &
            leta(4), ' s '
    end if
    close(un)
  end subroutine dftim

end module m_dftim
