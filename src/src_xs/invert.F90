
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module invert
  implicit none
contains
  
  subroutine zinvert_lapack(m,mi)
    ! arguments
    complex(8), intent(in) :: m(:,:)
    complex(8), intent(out) :: mi(:,:)
    ! local variables
    character(*), parameter :: thisnam='zinvert_lapack'
    complex(8), allocatable :: zwork(:)
    integer, allocatable :: ipiv(:)
    integer :: lwork,info,sh(2),n
    sh=shape(m)
    n=sh(1)
    allocate(ipiv(n))
    lwork=2*n
    allocate(zwork(lwork))
    mi(:,:)=m(:,:)
    call zgetrf(n,n,mi,n,ipiv,info)
    if (info.ne.0) then
       write(*,*)
       write(*,'("Error(",a,"): zgetrf returned non-zero info : ",I8)') &
            thisnam,info
       write(*,*)
       call terminate
    end if
    call zgetri(n,mi,n,ipiv,zwork,lwork,info)
    if (info.ne.0) then
       write(*,*)
       write(*,'("Error(",a,"): zgetri returned non-zero info : ",I8)') &
            thisnam,info
       write(*,*)
       call terminate
    end if
    deallocate(ipiv,zwork)
  end subroutine zinvert_lapack

end module invert
