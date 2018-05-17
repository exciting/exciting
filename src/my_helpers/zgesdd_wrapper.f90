! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine zgesdd_wrapper( mat, xsize, ysize, sval, lsvec, rsvec)
    implicit none
    ! arguments
    integer, intent( in) :: xsize, ysize
    complex (8), intent( in) :: mat( xsize, ysize)
    complex (8), intent( out) :: lsvec( xsize, xsize), rsvec( ysize, ysize)
    real(8), intent( out) :: sval( min( xsize, ysize))
    ! local variables
    integer :: lwork, info, mn, mx
    complex(8) :: cpy( xsize, ysize)
    ! allocatable arrays
    complex(8), allocatable :: work(:)
    real(8), allocatable :: rwork(:)
    integer, allocatable :: iwork(:)
    ! check input
    if( size( shape( mat), 1) .ne. 2) then
      write(*,*) "The Dimension of mat is not equal to 2. Please insert a 2-dimensional matrix."
      call terminate
    end if
    ! LAPACK 3.0 call
    cpy = mat
    lwork= - 1
    mn = min( xsize, ysize)
    mx = max( xsize, ysize)
    allocate( work(1))
    allocate( rwork( max( 5*mn*mn + 5*mn, 2*mx*mn + 2*mn*mn + mn)))
    allocate( iwork( 8*mn))
    call zgesdd( 'A', xsize, ysize, &
         cpy, xsize, &
         sval, &
         lsvec, xsize, &
         rsvec, ysize, &
         work, -1, rwork, iwork, info)
    lwork = nint( dble( work(1)))
    if( allocated( work)) deallocate( work)

    allocate( work( lwork))

    call zgesdd( 'A', xsize, ysize, &
         cpy, xsize, &
         sval, &
         lsvec, xsize, &
         rsvec, ysize, &
         work, lwork, rwork, iwork, info)

    if( info .ne. 0) then
      write(*,*)
      write(*, '("Error( zgesdd_wrapper): zgesdd returned non-zero info: ", i6)') info
      write (*,*)
      call terminate
    end if
    deallocate( work, rwork, iwork)
end subroutine zgesdd_wrapper
