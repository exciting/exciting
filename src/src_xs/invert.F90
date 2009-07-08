

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module invert
  implicit none
contains


subroutine zinvert_lapack(m, mi)
    implicit none
    ! arguments
    complex(8), intent(in) :: m(:, :)
    complex(8), intent(out) :: mi(:, :)
    ! local variables
    character(*), parameter :: thisnam='zinvert_lapack'
    complex(8), allocatable :: zwork(:)
    integer, allocatable :: ipiv(:)
    integer :: lwork, info, sh(2), n
    sh=shape(m)
    n=sh(1)
    allocate(ipiv(n))
    lwork=2*n
    allocate(zwork(lwork))
    mi(:, :)=m(:, :)
    call zgetrf(n, n, mi, n, ipiv, info)
    if (info.ne.0) then
       write(*, *)
       write(*, '("Error(", a, "): zgetrf returned non-zero info : ", I8)') &
	    thisnam, info
       write(*, *)
       call terminate
    end if
    call zgetri(n, mi, n, ipiv, zwork, lwork, info)
    if (info.ne.0) then
       write(*, *)
       write(*, '("Error(", a, "): zgetri returned non-zero info : ", I8)') &
	    thisnam, info
       write(*, *)
       call terminate
    end if
    deallocate(ipiv, zwork)
  end subroutine zinvert_lapack


subroutine zinvert_hermitian(flag, m, mi)
    implicit none
    ! arguments
    integer, intent(in) :: flag
    complex(8), intent(in) :: m(:, :)
    complex(8), intent(out) :: mi(:, :)
    ! local variables
    character(*), parameter :: thisnam='zinvert_hermitian'
    character(1) :: uplo
    integer :: info, n, j, sh(2)
    complex(8), allocatable :: tm(:, :)
    ! we do not check if both arguments have same shapes and are square matrices
    sh=shape(m)
    n=sh(1)
    allocate(tm(sh(1), sh(2)))
    tm(:, :)=m(:, :)
    info=0
    select case(flag)
    case(0)
       ! invert full matrix (matrix is allowed to be not strictly Hermitian)
       call zinvert_lapack(tm, mi)
    case(1)
       ! Hermitian average matrix
       tm=0.5d0*(tm+conjg(transpose(tm)))
       uplo='u'
    case(2)
       ! assume Hermitian and use upper triangle for inversion
       uplo='u'
    case(3)
       ! assume Hermitian and use lower triangle for inversion
       uplo='l'
    case default
       write(*, *)
       write(*, '("Error(", a, "): not a valid flag:", i6)') trim(thisnam), flag
       write(*, *)
       call terminate
    end select
    select case(flag)
    case(1, 2, 3)
       ! set up unity matrix for zposv
       mi(:, :)=(0.d0, 0.d0)
       forall(j=1:n)
	  mi(j, j)=(1.d0, 0.d0)
       end forall
       ! invert using upper/lower triangle of matrix
       call zposv(uplo, n, n, tm, n, mi, n, info)
    end select
    if (info.ne.0) then
       write(*, *)
       write(*, '("Error(", a, "): zposv returned non-zero info : ", I8)') &
	    trim(thisnam), info
       write(*, *)
       call terminate
    end if
    deallocate(tm)
  end subroutine zinvert_hermitian

end module invert
