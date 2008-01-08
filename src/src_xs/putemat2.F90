
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_putemat2
  implicit none
contains

  subroutine putemat2(iq,ik,filnam,x1,x2)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    character(*), intent(in) :: filnam
    complex(8), intent(in) :: x1(:,:,:)
    complex(8), optional, intent(in) :: x2(:,:,:)
    ! local variables
    integer :: un,recl
    call getunit(un)
    if (present(x2)) then
       ! I/O record length
       inquire(iolength=recl) nstval,nstcon,nkpt,ngq(iq),vql(:,iq),vkl(:,ik), &
            x1,x2
       open(unit=un,file=trim(filnam),form='unformatted', &
            action='write',access='direct',recl=recl)
       write(un,rec=ik) nstval,nstcon,nkpt,ngq(iq),vql(:,iq),vkl(:,ik),x1,x2
    else
       ! I/O record length
       inquire(iolength=recl) nstval,nstcon,nkpt,ngq(iq),vql(:,iq),vkl(:,ik), &
            x1
       open(unit=un,file=trim(filnam),form='unformatted', &
            action='write',access='direct',recl=recl)
       write(un,rec=ik) nstval,nstcon,nkpt,ngq(iq),vql(:,iq),vkl(:,ik),x1
    end if
    close(un)
  end subroutine putemat2

end module m_putemat2
