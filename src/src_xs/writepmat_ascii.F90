
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writepmat_ascii
  use modmain
  use m_getunit
  use m_getpmat
  implicit none
  complex(8), allocatable :: pmat(:,:,:)
  character(16) :: f1,f2,f
  integer :: un,ik,ist1,ist2

  call init0
  call init1
  call init2xs

  allocate(pmat(3,nstsv,nstsv))
  call getunit(un)
  open(un,file='PMAT_TD_ASC.OUT',action='write')

  do ik=1,nkpt
     ! read momentum matrix elements if required
     call getpmat(ik,vkl,.true.,'PMAT_TD.OUT',pmat)
     do ist1=1,nstsv
        f1='v'
        if (ist1.gt.(nstsv-nempty-1)) f1='c'
        do ist2=1,nstsv
           f2='v'
           if (ist2.gt.(nstsv-nempty-1)) f2='c'
           f='  '//trim(f1)//'-'//trim(f2)//'  '
           write(un,'(3i8,a,3g18.10)') ik,ist1,ist2,f,abs(pmat(:,ist1,ist2))
        end do
     end do
  end do
  
  close(un)
  deallocate(pmat)

end subroutine writepmat_ascii
