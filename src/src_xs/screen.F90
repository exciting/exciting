
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine screen
  use modmain
  use modmpi
  use modxs
  use m_filedel
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='screen'
  integer :: taskt,iq

  ! map variables for screening
  call initscr

  ! initialize universal variables
  call init0

  ! initialize k-point set
  call init1

  ! initialize q-point set
  call init2

  ! calculate eigenvectors, -values and occupancies for basic k-mesh
  taskt=task; task=1
  isreadstate0=.true.
  call genfilname(setfilext=.true.,dotext='_SCR.OUT')
  call gndstate
  task=taskt
  if (rank == 0) then
     ! safely remove unnecessary files
     call filedel('EQATOMS'//trim(filext))
     call filedel('EVALCORE'//trim(filext))
     call filedel('FERMIDOS'//trim(filext))
     call filedel('GEOMETRY'//trim(filext))
     call filedel('LATTICE'//trim(filext))
     call filedel('IADIST'//trim(filext))
     call filedel('LINENGY'//trim(filext))
     call filedel('SYMCRYS'//trim(filext))
     call filedel('SYMLAT'//trim(filext))
     call filedel('SYMSITE'//trim(filext))
     call filedel('TOTENERGY'//trim(filext))
     call filedel('EVALFV'//trim(filext))
     call filedel('RMSDVEFF'//trim(filext))
  end if

  ! *** TEST ***
  do iq=1,nqpt
     call init1
     call updateq(iq)
     write(*,'(a,i6,3f12.3,3x,3f12.3)') 'TEST: iq/vql/vqlcu',iq,vql(:,iq),vqlcu
  end do

  write(unitout,'(a)') "Info("//trim(thisnam)//"): Screening finished"

end subroutine screen
