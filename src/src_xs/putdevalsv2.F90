
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_putdevalsv2
  implicit none
contains

  subroutine putdevalsv2(iq,ik,tarec,filnam,e1,o1,e2,o2)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    logical :: tarec
    character(*) :: filnam
    real(8), intent(in) :: e1(:,:),o1(:,:)
    real(8), optional, intent(in) :: e2(:,:),o2(:,:)
    ! local variables
    integer :: un, ikr, recl
    ! record position for k-point
    ikr=ik
    if (.not.tarec) call getridx(procs,nkpt,ik,ikr)
    call getunit(un)
    if (present(e2).and.(.not.present(o2))) then
       write(*,*)
       write(*,'("Error(putdevalsv2): either type 1 or type 2 arguments &
            required for eigenvalue and occupancy differences - refer to &
            developers manual")')
       write(*,*)
    end if
    if (present(e2)) then
       ! I/O record length
       inquire(iolength=recl) nstval,nstcon,nkpt,ngq(iq),vql(:,iq), &
            vkl(:,ik),e1,o1,e2,o2
       open(unit=un,file=trim(filnam),form='unformatted', &
            action='write',access='direct',recl=recl)
       write(un,rec=ikr) nstval,nstcon,nkpt,ngq(iq),vql(:,iq),vkl(:,ik), &
            e1,o1,e2,o2
    else
       ! I/O record length
       inquire(iolength=recl) nstval,nstcon,nkpt,ngq(iq),vql(:,iq), &
            vkl(:,ik),e1,o1
       open(unit=un,file=trim(filnam),form='unformatted', &
            action='write',access='direct',recl=recl)
       write(un,rec=ikr) nstval,nstcon,nkpt,ngq(iq),vql(:,iq),vkl(:,ik), &
            e1,o1
    end if
    close(un)
  end subroutine putdevalsv2

end module m_putdevalsv2
