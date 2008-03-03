
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine atkp
  use modmain
  implicit none
  ! local variables
  real(8), allocatable :: evalsv_(:)
  real(8) vkl_(3)
  integer recl,nstsv_,ios,ik

  call init0
  call init1

  ! find number of bands
  inquire(iolength=recl) vkl_,nstsv_
  open(70,file=trim(scrpath)//'EVALSV'//trim(filext),action='READ', &
       form='UNFORMATTED',access='DIRECT',recl=recl)
  read(70,rec=1,iostat=ios) vkl_,nstsv_
  if (ios /= 0) then
     write(*,*) 'Cannot read in EVALSV.OUT'
     stop
  end if
  close(70)

  ! warn if different number of bands in file and input
  if (nstsv_ /= nstsv) then
     write(*,'(a,2i9)') 'Warning: different number of bands (current/&
          &EVALSV.OUT):',nstsv,nstsv_
  end if

  ! find current k-point
  allocate(evalsv_(nstsv_))
  inquire(iolength=recl) vkl_,nstsv_,evalsv_
  open(70,file=trim(scrpath)//'EVALSV'//trim(filext),action='READ', &
       form='UNFORMATTED',access='DIRECT',recl=recl)
  do ik=1,nkpt
     read(70,rec=ik,iostat=ios) vkl_,nstsv_,evalsv_
     if (any(vkl(:,ik) /= vkl_)) exit
  end do
  close(70)
  deallocate(evalsv_)
  write(*,'(a,i9)') 'EVALSV.OUT: at k-point',ik-1

end subroutine atkp

