
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_putemat
  implicit none
contains

  subroutine putemat(iq,ik,tarec,filnam,xou,xuo)
    use modmain
    use modxs
    use modmpi
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    logical :: tarec
    character(*), intent(in) :: filnam
    complex(8), intent(in) :: xou(:,:,:), xuo(:,:,:)
    ! local variables
    integer :: un, recl, ikr

    ! record position for k-point
    ikr=ik
    if (.not.tarec) call getridx(procs,nkpt,ik,ikr)
    ! I/O record length
    inquire(iolength=recl) nstval, nstcon, nkpt, ngq(iq), vql(:,iq), &
         vkl(:,ik), xou, xuo
    call getunit(un)
    open(unit=un,file=trim(filnam),form='unformatted', &
         action='write',access='direct',recl=recl)
    write(un,rec=ikr) nstval, nstcon, nkpt, ngq(iq), vql(:,iq), vkl(:,ik), &
         xou, xuo
    close(un)

  end subroutine putemat

end module m_putemat
