
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

program atkp
  implicit none
  ! local variables
  character(256) :: fname
  real(8), allocatable :: evalsv_(:)
  real(8) vkl_(3)
  integer recl,nstsv_,ios,ik

  fname='EVALSV.OUT'
  if (iargc()>1) then
     write(*,*)
     write(*,'("Usage: exciting [INPUT FILE]")')
     write(*,*)
     stop
  end if
  if (iargc().eq.1) then
     call getarg(1,fname)
  end if

  ! check for file
  inquire(file=trim(fname),exist=ios)
  if (.not.ios) then
    write(*,*)
    write(*,'(a)') 'No such file: '//trim(fname)
    write(*,*)
    stop
  end if
  write(*,'(a)') 'Analyzing file: '//trim(fname)

  ! find number of bands
  inquire(iolength=recl) vkl_,nstsv_
  open(70,file=trim(fname),action='READ', &
       form='UNFORMATTED',access='DIRECT',recl=recl)
  read(70,rec=1,iostat=ios) vkl_,nstsv_
  if (ios /= 0) then
     write(*,*)
     write(*,'(a)') 'Cannot read in '//trim(fname)
     write(*,*)
    stop
  end if
  close(70)

  ! report number of bands
  write(*,'(a,i8)') 'Number of bands: ',nstsv_
  write(*,'(a)') 'ik, vkl below'

  ! report k-points
  allocate(evalsv_(nstsv_))
  inquire(iolength=recl) vkl_,nstsv_,evalsv_
  open(70,file=trim(fname),action='READ', &
       form='UNFORMATTED',access='DIRECT',recl=recl)
  do ik=1,1000000000
     read(70,rec=ik,iostat=ios) vkl_,nstsv_,evalsv_
     if (ios /= 0) exit
     write(*,'(i8,3g18.10)') ik,vkl_
  end do
  
  close(70)
  deallocate(evalsv_)

end program atkp

