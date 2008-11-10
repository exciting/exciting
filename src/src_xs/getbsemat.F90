
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getbsemat(fname,ikkp,n1,n2,zmat)
  use m_getunit
  implicit none
  ! arguments
  character(*), intent(in) :: fname
  integer, intent(in) :: ikkp,n1,n2
  complex(8), intent(out) :: zmat(n1,n2,n1,n2)
  ! local variables
  integer :: un,recl,ikkp_,iknr_,jknr_,iq_,iqr_,n1_,n2_,n3_,n4_
  complex(8), allocatable :: zm(:,:,:,:)
  call getunit(un)
  ! get sizes of stored matrix
  inquire(iolength=recl) ikkp_,iknr_,jknr_,iq_,iqr_,n1_,n2_,n3_,n4_
  open(unit=un,file=trim(fname),form='unformatted',action='read', &
       access='direct',recl=recl)
  read(un,rec=1) ikkp_,iknr_,jknr_,iq_,iqr_,n1_,n2_,n3_,n4_
  close(un)
  ! check if requested size can be retrieved
  if ((n1.gt.n1_).or.(n2.gt.n2_)) then
     write(*,*)
     write(*,'("Error(getbsemat): requested matrix size out of range")')
     write(*,'(" requested size : ",2i8)') n1,n2
     write(*,'(" stored size    : ",2i8)') n1_,n2_
     write(*,*)
     call terminate
  end if
  ! read matrix
  allocate(zm(n1_,n2_,n3_,n4_))
  inquire(iolength=recl) ikkp_,iknr_,jknr_,iq_,iqr_,n1_,n2_,n3_,n4_,zm
  open(unit=un,file=trim(fname),form='unformatted',action='read', &
       access='direct',recl=recl)
  read(un,rec=ikkp) ikkp_,iknr_,jknr_,iq_,iqr_,n1_,n2_,n3_,n4_,zm
  close(un)
  ! check kkp-index
  if (ikkp.ne.ikkp_) then
     write(*,*)
     write(*,'("Error(getbsemat): inconsistent (k,kp)-index")')
     write(*,'(" requested                    : ",i8)') ikkp
     write(*,'(" stored at requested position : ",i8)') ikkp_
     write(*,*)
     call terminate
  end if
  ! cut matrix
  zmat(:,:,:,:)=zm(n1_-n1+1:,:n2,n1_-n1+1:,:n2)
  deallocate(zm)
end subroutine getbsemat
