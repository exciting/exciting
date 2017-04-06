!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine diaghermat ( matsize, mat, eval, evec)
      Implicit None
  ! arguments
      Integer, Intent (In) :: matsize
      Complex (8), Intent (In) :: mat( matsize, matsize)
      Real (8), Intent (Out) :: eval( matsize)
      Complex (8), Intent (Out) :: evec( matsize, matsize)
  ! local variables
      Real (8) :: vl, vu, abstol
      Integer :: il, iu, neval, lwork, info, lrwork, liwork
      complex(8) :: matcpy( matsize, matsize)
  ! allocatable arrays
      Complex (8), Allocatable :: work (:)
      Real (8), Allocatable :: rwork (:)
      Integer, Allocatable :: iwork (:), ifail (:),isuppz(:)
  ! external functions
      Real (8), External :: dlamch
  ! check input
    IF( SIZE( SHAPE( mat), 1).NE.2) THEN
      write(*,*) "The Dimension of mat is not equal to 2. Please insert a 2-dimensional matrix."
      stop
    END IF
    IF( SIZE( mat, 2).NE.matsize) THEN
      write(*,*) "The Matrix mat is not a square matrix. Please insert a sqare matrix."
      stop
    END IF
    matcpy = mat
  ! smallest and largest eigenvalue indices
      il = 1
      iu = matsize
  ! tolerance parameter
      abstol = 2.d0 * dlamch ('S')
  ! LAPACK 3.0 call
     lrwork=-1 !hamsiz
     liwork=-1 !5*hamsiz
     lwork=-1
     iu=matsize
     Allocate (work(1), rwork(1), iwork(1), isuppz(1))
     Call zheevr ('V', 'A', 'U', matsize, matcpy, matsize, vl, vu, il, iu, &
     & abstol, neval, eval, evec, matsize, isuppz, work, lwork, rwork, lrwork, iwork, liwork, &
     & info)
     lrwork=int(rwork(1))
     liwork=int(iwork(1))
     lwork=int(work(1))
     deallocate(work, rwork, iwork, isuppz)

!     write(*,*) lrwork,liwork,lwork
      Allocate (work(lwork), rwork(lrwork), iwork(liwork))
      allocate(isuppz(matsize*2))

      Call zheevr ('V', 'A', 'U', matsize, matcpy, matsize, vl, vu, il, iu, &
     & abstol, neval, eval, evec, matsize, isuppz, work, lwork, rwork, lrwork, iwork, liwork, &
     & info)
      deallocate(isuppz)

      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("ERROR( diaghermat): zheevr returned non-zero info:&
        &", i6)') info
         Write (*,*)
         Call terminate
      End If
  ! deallocate Hamiltonian array
      Deallocate (work, rwork, iwork)

End Subroutine diaghermat

