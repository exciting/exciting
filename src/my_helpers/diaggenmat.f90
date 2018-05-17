!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine diaggenmat ( matsize, mat, eval, evec)
      Implicit None
  ! arguments
      Integer, Intent (In) :: matsize
      Complex (8), Intent (In) :: mat( matsize, matsize)
      Complex (8), Intent (Out) :: eval( matsize)
      Complex (8), Intent (Out) :: evec( matsize, matsize)
  ! local variables
      Integer :: lwork, info
      Complex (8) :: vl( matsize)
  ! allocatable arrays
      Complex (8), Allocatable :: work (:), rwork(:)
  ! check input
    IF( SIZE( SHAPE( mat), 1).NE.2) THEN
      write(*,*) "The Dimension of mat is not equal to 2. Please insert a 2-dimensional matrix."
      stop
    END IF
    IF( SIZE( mat, 2).NE.matsize) THEN
      write(*,*) "The Matrix mat is not a square matrix. Please insert a sqare matrix."
      stop
    END IF
  ! LAPACK 3.0 call
     lwork=-1
     Allocate (work(1), rwork( 2*matsize))
     Call zgeev ('N', 'V', matsize, mat, matsize, eval, vl, matsize, evec, matsize, work, lwork, rwork, info)
     lwork=int( work(1))
     deallocate( work)

!     write(*,*) lrwork,liwork,lwork
      Allocate( work( lwork))

     Call zgeev ('N', 'V', matsize, mat, matsize, eval, vl, matsize, evec, matsize, work, lwork, rwork, info)

      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error( diaggenmat): zgeev returned non-zero info:&
        &", i6)') info
         Write (*,*)
         Call terminate
      End If
  ! deallocate Hamiltonian array
      Deallocate( work)

End Subroutine diaggenmat

