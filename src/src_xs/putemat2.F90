
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
    complex(8), intent(in) :: x1(:,:,:),x2(:,:,:)
    ! local variables
    integer :: un,recl
    ! I/O record length
    inquire(iolength=recl) nstval,nstcon,nkpt,ngq(iq),vql(:,iq),vkl(:,ik),x1,x2
    call getunit(un)
    open(unit=un,file=trim(filnam),form='unformatted', &
         action='write',access='direct',recl=recl)
    write(un,rec=ik) nstval,nstcon,nkpt,ngq(iq),vql(:,iq),vkl(:,ik),x1,x2
    close(un)

write(*,*) 'putemat2:abs(x1),abs(x2)',sum(abs(x1)),sum(abs(x2))
write(*,*) 'putemat2:shape(x1),shape(x2)',shape(x1),shape(x2)


  end subroutine putemat2

end module m_putemat2
