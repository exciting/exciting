!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine diagrealmat ( matsize, mat, eval, evec)
    implicit none

    integer, intent (in) :: matsize
    real (8), intent (in) :: mat( matsize, matsize)
    complex (8), intent (out) :: eval( matsize)
    complex (8), intent (out) :: evec( matsize, matsize)

    integer :: lwork, info, i, j
    real(8) :: matcpy( matsize, matsize)
    real(8) :: vl( matsize, matsize), vr( matsize, matsize), wr( matsize), wi( matsize)
    real(8), allocatable :: work (:)
  
    ! check input
    if( size( shape( mat), 1).ne.2) then
      write(*,*) "Error (diagrealmat): The Dimension of mat is not equal to 2. Please insert a 2-dimensional matrix."
      stop
    end if
    if( size( mat, 2).ne.matsize) then
      write(*,*) "Error (diagrealmat): The Matrix mat is not a square matrix. Please insert a sqare matrix."
      stop
    end if
    
    matcpy = mat

    lwork = -1
    allocate( work(1))
    call dgeev( 'N', 'V', matsize, matcpy, matsize, wr, wi, vl, matsize, vr, matsize, work, lwork, info)
    lwork = int( work(1))
    deallocate( work)
    allocate( work( lwork))
    call dgeev( 'N', 'V', matsize, matcpy, matsize, wr, wi, vl, matsize, vr, matsize, work, lwork, info)

    evec = cmplx( 0, 0, 8)
    eval = cmplx( 0, 0, 8)
    do i = 1, matsize
      if( (abs( eval( i)) .lt. 1.d-6) .and. (wi( i) .lt. 1.d-6)) then
        eval( i) = cmplx( wr( i), 0, 8)
        evec( :, i) = vr( :, i)
      else
        if( abs( eval( i)) .lt. 1.d-6) then
          eval( i) = cmplx( wr( i), wi( i), 8)
          eval( i+1) = conjg( eval( i))
        end if
        if( norm2( abs( evec( :, i))) .lt. 1.d-6) then
          evec( :, i) = cmplx( vr( :, i), vr( :, i+1), 8)
          evec( :, i+1) = conjg( evec( :, i))
        end if
      end if
    end do

    if( info .ne. 0) then
       write(*,*)
       write(*, '("Error( diagrealmat): dgeev returned non-zero info:", i6)') info
       write (*,*)
       call terminate
    end if
    deallocate( work)

End Subroutine diagrealmat

