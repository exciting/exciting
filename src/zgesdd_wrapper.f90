!
!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine zgesdd_wrapper( mat, xsize, ysize, sval, lsvec, rsvec)
      Implicit None
  ! arguments
      Integer, Intent (In) :: xsize, ysize
      Complex (8), Intent (In) :: mat( xsize, ysize)
      Complex (8), Intent (Out) :: lsvec( xsize, xsize), rsvec( ysize, ysize)
      real(8), intent( out) :: sval( min( xsize, ysize))
  ! local variables
      Integer :: lwork, info
      complex(8) :: cpy( xsize, ysize)
  ! allocatable arrays
      Complex (8), Allocatable :: work (:)
      real(8), allocatable :: rwork(:)
      integer, allocatable :: iwork(:)
  ! check input
    IF( SIZE( SHAPE( mat), 1).NE.2) THEN
      write(*,*) "The Dimension of mat is not equal to 2. Please insert a 2-dimensional matrix."
      stop
    END IF
  ! LAPACK 3.0 call
     cpy(:,:) = mat(:,:)
     lwork=-1
     Allocate (work(1), &
               rwork( max( 5*min( xsize, ysize)**2 + 5*min( xsize, ysize), 2*max( xsize, ysize)*min( xsize, ysize) + 2*min( xsize, ysize)**2 + min( xsize, ysize))), &
               iwork( 8*min( xsize, ysize)))
     call zgesdd( 'A', xsize, ysize, &
          cpy, xsize, &
          sval, &
          lsvec, xsize, &
          rsvec, ysize, &
          work, -1, rwork, iwork, info)
     lwork = nint( abs( work(1)))
     deallocate( work)

!     write(*,*) lrwork,liwork,lwork
      Allocate( work( lwork))

      call zgesdd( 'A', xsize, ysize, &
           cpy, xsize, &
           sval, &
           lsvec, xsize, &
           rsvec, ysize, &
           work, lwork, rwork, iwork, info)

      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error( zgesdd_wrapper): zgesdd returned non-zero info:&
        &", i6)') info
         Write (*,*)
         Call terminate
      End If
  ! deallocate Hamiltonian array
      Deallocate( work, rwork, iwork)

End Subroutine zgesdd_wrapper
