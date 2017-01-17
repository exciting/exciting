!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine diagsymmat ( matsize, mat, eval, evec)
    implicit none

    integer, intent (in) :: matsize
    real (8), intent (in) :: mat( matsize, matsize)
    real (8), intent (out) :: eval( matsize)
    real (8), intent (out) :: evec( matsize, matsize)

    integer :: lwork, info, i, j
    real(8) :: matcpy( matsize, matsize), diag( matsize), offd( matsize-1), tau( matsize-1)
    real(8) :: v( matsize, 1), h( matsize, matsize)
    real(8), allocatable :: work (:)
  
    ! check input
    if( size( shape( mat), 1).ne.2) then
      write(*,*) "Error (diagsymmat): The Dimension of mat is not equal to 2. Please insert a 2-dimensional matrix."
      stop
    end if
    if( size( mat, 2).ne.matsize) then
      write(*,*) "Error (diagsymmat): The Matrix mat is not a square matrix. Please insert a sqare matrix."
      stop
    end if
    
    matcpy = mat

    lwork = -1
    allocate( work(1))
    call dsytrd( 'U', matsize, matcpy, matsize, diag, offd, tau, work, lwork, info)
    lwork = int( work(1))
    deallocate( work)
    allocate( work( lwork))
    call dsytrd( 'U', matsize, matcpy, matsize, diag, offd, tau, work, lwork, info)

    evec = 0.d0
    do i = 1, matsize
      evec( i, i) = 1.d0
    end do
    do i = matsize-1, 1, -1
      v = 0.d0
      v( i, 1) = 1.d0
      if( i .gt. 1) v( 1:i-1, 1) = matcpy( 1:i-1, i+1)
      h = -tau( i)*matmul( v, transpose( v))
      do j = 1, matsize
        h( j, j) = h( j, j) + 1.d0
      end do
      evec = matmul( evec, h)
    end do

    deallocate( work)
    allocate( work( max( 1, 2*matsize-2)))
    eval = diag
    call dsteqr( 'V', matsize, eval, offd, evec, matsize, work, info)

    if( info .ne. 0) then
       write(*,*)
       write(*, '("Error( diagsymmat): dsteqr returned non-zero info:", i6)') info
       write (*,*)
       call terminate
    end if
    deallocate( work)

End Subroutine diagsymmat

