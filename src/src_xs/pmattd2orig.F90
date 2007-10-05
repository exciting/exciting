
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine pmattd2orig
  use modmain
  use modtddft
  use modpar
  use m_getpmat
  use m_getunit
  implicit none
  ! local variables
  character(*),parameter :: thisnam='pmattd2orig'
  complex(8), allocatable :: pm(:,:,:)
  integer :: un,ik,recl

  if (rank.eq.1) then
     call init0
     call init1
     call init2td
     allocate(pm(3,nstsv,nstsv))

     inquire(iolength=recl) pm
     call getunit(un)
     open(un,file='PMAT.OUT',form='unformatted',action='write',&
          status='replace',access='direct',recl=recl)

     do ik=1,nkpt
        call getpmat(ik,vkl,.true.,'PMAT_TD.OUT',pm)
        write(un,rec=ik) pm
     end do

     close(un)
     deallocate(pm)
  end if

  call getunit(un)
  call barrier(rank=rank,nproc=nproc,un=un,async=0,string='.barrier')

  write(unitout,'(a)') "Info("//trim(thisnam)//"): conversion of PMAT &
       &finished"

end subroutine pmattd2orig
