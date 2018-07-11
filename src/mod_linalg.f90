module m_linalg
  use mod_constants
  implicit none
  real(8), external :: dlamch

  contains

    !***********************************
    !          DIAGONALIZATION
    !***********************************
    ! diagonalize complex general matrix
    subroutine zgediag( mat, eval, evec)
      complex(8), intent( in) :: mat(:,:)
      complex(8), intent( out) :: eval(:), evec(:,:)

      integer :: m, info, lwork
      complex(8), allocatable :: cpy(:,:), work(:), vl(:)
      real(8), allocatable :: rwork(:)

      m = size( mat, 1)

      if( m .le. 0) then
        write(*,'("Error (zgediag): Invalid matrix dimension.")')
        stop
      end if
      if( size( mat, 2) .ne. m) then
        write(*,'("Error (zgediag): The matrix must be squared.")')
        stop
      end if
      if( size( eval, 1) .ne. m) then
        write(*,'("Error (zgediag): The eigenvalue array must have the same dimension as the matrix.")')
        stop
      end if
      if( (size( evec, 1) .ne. m) .or. (size( evec, 2) .ne. m)) then
        write(*,'("Error (zgediag): The eigenvector array must have the same dimension as the matrix.")')
        stop
      end if

      allocate( cpy( m, m))
      allocate( vl( m))
      allocate( work(1), rwork( 2*m))

      cpy = mat

      call zgeev( 'n', 'v', m, cpy, m, eval, vl, m, evec, m, work, -1, rwork, info)
      lwork = nint( dble( work( 1)))
      deallocate( work)
      allocate( work( lwork))
      call zgeev( 'n', 'v', m, cpy, m, eval, vl, m, evec, m, work, lwork, rwork, info)

      if( info .ne. 0) then
        write(*,'("Error (zgediag): Diagonalization failed. ZGEEV returned info ",i4)') info
        stop
      end if

      deallocate( cpy, work, rwork, vl)
      return
    end subroutine zgediag

    ! diagonalize complex hermitian matrix
    subroutine zhediag( mat, eval, evec)
      complex(8), intent( in) :: mat(:,:)
      real(8), intent( out) :: eval(:)
      complex(8), intent( out) :: evec(:,:)

      integer :: m, n, info, lrwork, liwork, lwork
      real(8) :: vl, vu, eps
      integer, allocatable :: iwork(:), isuppz(:)
      complex(8), allocatable :: cpy(:,:), work(:)
      real(8), allocatable :: rwork(:)

      eps = 2.d0*dlamch( 's')
      m = size( mat, 1)

      if( m .le. 0) then
        write(*,'("Error (zhediag): Invalid matrix dimension (",i6,")")'), m
        stop
      end if
      if( size( mat, 2) .ne. m) then
        write(*,'("Error (zhediag): The matrix must be squared.")')
        stop
      end if
      if( size( eval, 1) .ne. m) then
        write(*,'("Error (zhediag): The eigenvalue array must have the same dimension as the matrix.")')
        stop
      end if
      if( (size( evec, 1) .ne. m) .or. (size( evec, 2) .ne. m)) then
        write(*,'("Error (zhediag): The eigenvector array must have the same dimension as the matrix.")')
        stop
      end if

      allocate( cpy( m, m))
      allocate( work(1), rwork(1), iwork(1), isuppz( 2*m))

      cpy = mat

      call zheevr( 'v', 'a', 'u', m, cpy, m, vl, vu, 1, m, eps, n, eval, evec, m, isuppz, work, -1, rwork, -1, iwork, -1, info)
      lrwork = nint( rwork( 1))
      liwork = iwork( 1)
      lwork = nint( dble( work( 1)))
      deallocate( work, rwork, iwork)
      allocate( work( lwork), rwork( lrwork), iwork( liwork))
      call zheevr( 'v', 'a', 'u', m, cpy, m, vl, vu, 1, m, eps, n, eval, evec, m, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)

      if( info .ne. 0) then
        write(*,'("Error (zhediag): Diagonalization failed. ZHEEVR returned info ",i4)') info
        stop
      end if

      deallocate( cpy, work, rwork, iwork, isuppz)
      return
    end subroutine zhediag

    ! diagonalize real symmetric matrix
    subroutine rsydiag( mat, eval, evec)
      real(8), intent( in) :: mat(:,:)
      real(8), intent( out) :: eval(:), evec(:,:)

      integer :: m, i, j, info, lwork
      integer, allocatable :: iwork(:), isuppz(:)
      real(8), allocatable :: cpy(:,:), work(:), diag(:), offd(:), tau(:), v(:,:), h(:,:)

      m = size( mat, 1)

      if( m .le. 0) then
        write(*,'("Error (zhediag): Invalid matrix dimension.")')
        stop
      end if
      if( size( mat, 2) .ne. m) then
        write(*,'("Error (zhediag): The matrix must be squared.")')
        stop
      end if
      if( size( eval, 1) .ne. m) then
        write(*,'("Error (zhediag): The eigenvalue array must have the same dimension as the matrix.")')
        stop
      end if
      if( (size( evec, 1) .ne. m) .or. (size( evec, 2) .ne. m)) then
        write(*,'("Error (zhediag): The eigenvector array must have the same dimension as the matrix.")')
        stop
      end if

      allocate( cpy( m, m))
      allocate( work(1), diag( m), offd( m-1), tau( m-1), v( m, 1), h( m, m))

      cpy = mat

      call dsytrd( 'u', m, cpy, m, diag, offd, tau, work, -1, info)
      lwork = int( work(1))
      deallocate( work)
      allocate( work( lwork))
      call dsytrd( 'u', m, cpy, m, diag, offd, tau, work, lwork, info)

      if( info .ne. 0) then
        write(*,'("Error (rsydiag): Tridiagonalization failed. DSYTRD returned info ",i4)') info
        stop
      end if

      evec = 0.d0
      do i = 1, m
        evec( i, i) = 1.d0
      end do
      do i = m-1, 1, -1
        v = 0.d0
        v( i, 1) = 1.d0
        if( i .gt. 1) v( 1:i-1, 1) = cpy( 1:i-1, i+1)
        h = -tau( i)*matmul( v, transpose( v))
        do j = 1, m
          h( j, j) = h( j, j) + 1.d0
        end do
        evec = matmul( evec, h)
      end do

      deallocate( work)
      allocate( work( max( 1, 2*m-2)))
      eval = diag
      call dsteqr( 'v', m, eval, offd, evec, m, work, info)

      if( info .ne. 0) then
        write(*, '("Error( rsydiag): Diagonalisation failed. DSTEQR returned info", i4)') info
        stop
      end if

      deallocate( cpy, work, offd, tau, v, h)
      return
    end subroutine rsydiag

    !***********************************
    !    SINGULAR VALUE DECOMPOSITION
    !***********************************
    ! SVD of complex matrix
    subroutine zsvd( mat, sval, lsvec, rsvec)
      complex(8), intent( in) :: mat(:,:)
      real(8), intent( out) :: sval(:)
      complex(8), intent( out) :: lsvec(:,:), rsvec(:,:)

      integer :: m, n, l, k, info, lwork
      integer, allocatable :: iwork(:)
      real(8), allocatable :: rwork(:)
      complex(8), allocatable :: cpy(:,:), work(:)
    
      m = size( mat, 1)
      n = size( mat, 2)
      l = min( m, n)
      k = max( m, n)
    
      if( (m .le. 0) .or. (n .le. 0)) then
        write(*,'("Error (zsvd): Invalid matrix dimension.")')
        stop
      end if
      if( size( sval, 1) .ne. l) then
        write(*,'("Error (zsvd): The singular value array must have the dimension of the smaller side of the matrix.")')
        stop
      end if
      if( (size( lsvec, 1) .ne. m) .or. (size( lsvec, 2) .ne. m)) then
        write(*,'("Error (zsvd): The left singular vector array must have the dimension of the first side of the matrix.")')
        stop
      end if
      if( (size( rsvec, 1) .ne. n) .or. (size( rsvec, 2) .ne. n)) then
        write(*,'("Error (zsvd): The right singular vector array must have the dimension of the second side of the matrix.")')
        stop
      end if

      allocate( cpy( m, n))
      allocate( iwork( 8*l))
      allocate( work(1))
      allocate( rwork( max( 5*l*(l+1), 2*l*(l+k) + l)))

      cpy = mat
      call zgesdd( 'a', m, n, cpy, m, sval, lsvec, m, rsvec, n, work, -1, rwork, iwork, info)
      lwork = nint( dble( work(1)))
      deallocate( work)
      allocate( work( lwork))
      call zgesdd( 'a', m, n, cpy, m, sval, lsvec, m, rsvec, n, work, lwork, rwork, iwork, info)
      if( info .ne. 0) then
        write(*,'("Error (zsvd): SVD failed. ZGESDD returned info ",i4)') info
        stop
      end if

      deallocate( cpy, iwork, rwork, work)
      return
    end subroutine zsvd

    ! SVD of real matrix
    subroutine rsvd( mat, sval, lsvec, rsvec)
      real(8), intent( in) :: mat(:,:)
      real(8), intent( out) :: sval(:)
      real(8), intent( out) :: lsvec(:,:), rsvec(:,:)

      integer :: m, n, l, k, info, lwork
      integer, allocatable :: iwork(:)
      real(8), allocatable :: cpy(:,:), work(:)
    
      m = size( mat, 1)
      n = size( mat, 2)
      l = min( m, n)
      k = max( m, n)
    
      if( (m .le. 0) .or. (n .le. 0)) then
        write(*,'("Error (zsvd): Invalid matrix dimension.")')
        stop
      end if
      if( size( sval, 1) .ne. l) then
        write(*,'("Error (zsvd): The singular value array must have the dimension of the smaller side of the matrix.")')
        stop
      end if
      if( (size( lsvec, 1) .ne. m) .or. (size( lsvec, 2) .ne. m)) then
        write(*,'("Error (zsvd): The left singular vector array must have the dimension of the first side of the matrix.")')
        stop
      end if
      if( (size( rsvec, 1) .ne. n) .or. (size( rsvec, 2) .ne. n)) then
        write(*,'("Error (zsvd): The right singular vector array must have the dimension of the second side of the matrix.")')
        stop
      end if

      allocate( cpy( m, n))
      allocate( iwork( 8*l))
      allocate( work(1))

      cpy = mat
      call dgesdd( 'a', m, n, cpy, m, sval, lsvec, m, rsvec, n, work, -1, iwork, info)
      lwork = nint( work(1))
      deallocate( work)
      allocate( work( lwork))
      call dgesdd( 'a', m, n, cpy, m, sval, lsvec, m, rsvec, n, work, lwork, iwork, info)
      if( info .ne. 0) then
        write(*,'("Error (rsvd): SVD failed. DGESDD returned info ",i4)') info
        stop
      end if

      deallocate( cpy, iwork, work)
      return
    end subroutine rsvd

    !***********************************
    !           PSEUDO INVERSE
    !***********************************
    ! compute pseudo inverse of real matrix
    subroutine zpinv( mat, pinv)
      complex(8), intent( in) :: mat(:,:)
      complex(8), intent( out) :: pinv(:,:)
    
      integer :: m, n, l, i
      real(8) :: eps
      real(8), allocatable :: sval(:)
      complex(8), allocatable :: lsvec(:,:), rsvec(:,:)
      
      m = size( mat, 1)
      n = size( mat, 2)
      l = min( m, n)
    
      if( (m .le. 0) .or. (n .le. 0)) then
        write(*,'("Error (rpinv): Invalid matrix dimension.")')
        stop
      end if
      if( (size( pinv, 1) .ne. n) .or. (size( pinv, 2) .ne. m)) then
        write(*,'("Error (rpinv): Matrix dimensions do not match. pinv has to have the dimension of mat^T.")')
        stop
      end if
    
      allocate( sval( l))
      allocate( lsvec( m, m), rsvec( n, n))
    
      call zsvd( mat, sval, lsvec, rsvec)
      eps = dlamch( 'e')*dble( max( m, n))*maxval( sval)
    
      if( m .ge. n) then
        do i = 1, l
          if( sval( i) .gt. eps) then
            lsvec( :, i) = lsvec( :, i)/cmplx( sval( i), 0.d0, 8)
          else
            lsvec( :, i) = 0.d0
          end if
        end do
      else
        do i = 1, l
          if( sval( i) .gt. eps) then
            rsvec( i, :) = rsvec( i, :)/cmplx( sval( i), 0.d0, 8)
          else
            rsvec( i, :) = 0.d0
          end if
        end do
      end if
      call zgemm( 'c', 'c', n, m, l, zone, rsvec, n, lsvec, m, zzero, pinv, n)
    
      deallocate( sval, lsvec, rsvec)
      return
    end subroutine zpinv

    ! compute pseudo inverse of real matrix
    subroutine rpinv( mat, pinv)
      real(8), intent( in) :: mat(:,:)
      real(8), intent( out) :: pinv(:,:)
    
      integer :: m, n, l, i
      real(8) :: eps
      real(8), allocatable :: sval(:), lsvec(:,:), rsvec(:,:)
      
      m = size( mat, 1)
      n = size( mat, 2)
      l = min( m, n)
    
      if( (m .le. 0) .or. (n .le. 0)) then
        write(*,'("Error (rpinv): Invalid matrix dimension.")')
        stop
      end if
      if( (size( pinv, 1) .ne. n) .or. (size( pinv, 2) .ne. m)) then
        write(*,'("Error (rpinv): Matrix dimensions do not match. pinv has to have the dimension of mat^T.")')
        stop
      end if
    
      allocate( sval( l))
      allocate( lsvec( m, m), rsvec( n, n))
    
      call rsvd( mat, sval, lsvec, rsvec)
      eps = dlamch( 'e')*dble( max( m, n))*maxval( sval)
    
      if( m .ge. n) then
        do i = 1, l
          if( sval( i) .gt. eps) then
            lsvec( :, i) = lsvec( :, i)/sval( i)
          else
            lsvec( :, i) = 0.d0
          end if
        end do
      else
        do i = 1, l
          if( sval( i) .gt. eps) then
            rsvec( i, :) = rsvec( i, :)/sval( i)
          else
            rsvec( i, :) = 0.d0
          end if
        end do
      end if
      call dgemm( 't', 't', n, m, l, 1.d0, rsvec, n, lsvec, m, 0.d0, pinv, n)
    
      deallocate( sval, lsvec, rsvec)
      return
    end subroutine rpinv
end module m_linalg
