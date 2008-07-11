
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine pmatxs2orig
  use modmain
  use modxs
  use modmpi
  use m_getunit
  use m_getpmat
  implicit none
  ! local variables
  character(*),parameter :: thisnam='pmatxs2orig'
  complex(8), allocatable :: pm(:,:,:)
  integer :: un,ik,recl
  if (rank == 0) then
     call init0
     call init1
     call init2xs
     allocate(pm(3,nstsv,nstsv))
     inquire(iolength=recl) pm
     call getunit(un)
     open(un,file='PMAT.OUT',form='unformatted',action='write',&
          status='replace',access='direct',recl=recl)
     do ik=1,nkpt
        call getpmat(ik,vkl,.true.,'PMAT_XS.OUT',pm)
        write(un,rec=ik) pm
     end do
     close(un)
     deallocate(pm)
  end if
  call barrier
  write(unitout,'(a)') "Info("//trim(thisnam)//"): conversion of PMAT &
       &to original format finished"
end subroutine pmatxs2orig
