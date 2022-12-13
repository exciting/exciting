! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! !REVISION HISTORY:
!   Created July 2020 (SeTi)
module m_linalg
  use constants, only: zzero, zone
  use xlapack, only: svd_divide_conquer

  implicit none
  real(8), external :: dlamch

  contains

    !***********************************
    !          DIAGONALIZATION
    !***********************************
    ! diagonalize complex general matrix
    subroutine zgediag( mat, eval, levec, revec)
      complex(8), intent( in) :: mat(:,:)
      complex(8), intent( out) :: eval(:)
      complex(8), optional, target, intent( out) :: levec(:,:), revec(:,:)

      integer :: m, info, lwork
      character :: lv, rv
      complex(8), pointer :: lv_ptr, rv_ptr
      complex(8), allocatable, target :: cpy(:,:)
      complex(8), allocatable :: work(:)
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
      lv = 'n'
      if( present( levec)) then
        lv = 'v'
        lv_ptr => levec(1,1)
        if( (size( levec, 1) .ne. m) .or. (size( levec, 2) .ne. m)) then
          write(*,'("Error (zgediag): The eigenvector array must have the same dimension as the matrix.")')
          stop
        end if
      end if
      rv = 'n'
      if( present( revec)) then
        rv = 'v'
        rv_ptr => revec(1,1)
        if( (size( revec, 1) .ne. m) .or. (size( revec, 2) .ne. m)) then
          write(*,'("Error (zgediag): The eigenvector array must have the same dimension as the matrix.")')
          stop
        end if
      end if

      allocate( cpy( m, m))
      allocate( work(1), rwork( 2*m))

      cpy = mat
      if( lv == 'n') lv_ptr => cpy(1,1)
      if( rv == 'n') rv_ptr => cpy(1,1)

      call zgeev( lv, rv, m, cpy, m, eval, lv_ptr, m, rv_ptr, m, work, -1, rwork, info)
      lwork = nint( dble( work( 1)))
      deallocate( work)
      allocate( work( lwork))
      call zgeev( lv, rv, m, cpy, m, eval, lv_ptr, m, rv_ptr, m, work, lwork, rwork, info)

      if( info .ne. 0) then
        write(*,'("Error (zgediag): Diagonalization failed. ZGEEV returned info ",i4)') info
        stop
      end if

      deallocate( cpy, work, rwork)
      return
    end subroutine zgediag

    ! diagonalize complex hermitian matrix
    subroutine zhediag( mat, eval, evec)
      complex(8), intent( in) :: mat(:,:)
      real(8), intent( out) :: eval(:)
      complex(8), optional, intent( out) :: evec(:,:)

      integer :: m, n, info, lwork
      character :: rv
      complex(8), allocatable :: cpy(:,:), work(:)
      real(8), allocatable :: rwork(:)

      m = size( mat, 1)

      if( m .le. 0) then
        write(*,'("Error (zhediag): Invalid matrix dimension (",i6,")")') m
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
      rv = 'n'
      if( present( evec)) then
        rv = 'v'
        if( (size( evec, 1) .ne. m) .or. (size( evec, 2) .ne. m)) then
          write(*,'("Error (zhediag): The eigenvector array must have the same dimension as the matrix.")')
          stop
        end if
      end if

      allocate( cpy( m, m))
      allocate( work(1), rwork( 3*m-2))

      cpy = mat

      call zheev( rv, 'u', m, cpy, m, eval, work, -1, rwork, info)
      lwork = nint( dble( work( 1)))
      deallocate( work)
      allocate( work( lwork))
      call zheev( rv, 'u', m, cpy, m, eval, work, lwork, rwork, info)
      if( rv == 'v') evec = cpy

      if( info .ne. 0) then
        write(*,'("Error (zhediag): Diagonalization failed. ZHEEV returned info ",i4)') info
        stop
      end if

      deallocate( cpy, work, rwork)
      return
    end subroutine zhediag

    ! find unique gauge of eigenvectors
    ! First, degenerate eigenvectors \({\bf U}_{:,m:n}\) are rotated using a
    ! \((n-m+1)\)-dimensional unitary matrix \({\bf Q}\) that is obtained from
    ! a QR decomposition \({\bf U}_{:,m:n}^\dagger = {\bf Q}^\dagger \cdot {\bf R}^\dagger\), 
    ! and the degenerate eigenvectors \({\bf U}_{:,m:n}\) are replaced by 
    ! \({\bf R} = {\bf U}_{:,m:n} \cdot {\bf Q}\).
    ! In a second step, all eigenvectors are multiplied with a constant phase factor
    ! such that the first occurence of their largest component is real and positive.
    subroutine zhegauge( eval, evec, eps)
      use m_plotmat
      real(8), intent( in)           :: eval(:)
      complex(8), intent( inout)     :: evec(:,:)
      real(8), optional, intent( in) :: eps

      integer :: m, n, i, j, k, l
      real(8) :: epsil, t1
      complex(8) :: z1
      complex(8), allocatable :: aux1(:,:), aux2(:,:)

      m = size( evec, 1)
      n = size( evec, 2)

      if( size( eval) < n) then
        write(*,'("Error (zhegauge): Less eigenvalues than eigenvectors.")')
        stop
      end if

      epsil = 1.d-10
      if( present( eps)) epsil = eps

      i = 1; j = 1
      do while( j <= n)
        if( j < n) then
          if( abs( eval(j+1) - eval(j)) < epsil) then
            j = j + 1
            cycle
          end if
        end if
        k = j - i + 1
        ! states i to j are degenerate
        if( k > 1) then
          ! get rotated eigenvector from QR decomposition
          allocate( aux1(k,m), aux2(k,m))
          aux1 = conjg( transpose( evec(:,i:j)))
          where( abs( aux1) < 1.d-12) aux1 = zzero
          call zqr( aux1, r=aux2)
          evec(:,i:j) = conjg( transpose( aux2))
          deallocate( aux1, aux2)
        end if
        ! rotate eigenvectors such that the first occurance of the largest component is real and positive
        do l = i, j
          t1 = maxval( abs( evec(:,l)))
          do k = 1, m
            if( abs( evec(k,l)) + epsil > t1) exit
          end do
          z1 = conjg( evec(k,l))/abs( evec(k,l))
          evec(:,l) = evec(:,l)*z1
        end do
        ! increment counters
        i = j + 1
        j = i
      end do
    end subroutine zhegauge

    ! diagonalize real symmetric matrix
    subroutine rsydiag( mat, eval, evec)
      real(8), intent( in) :: mat(:,:)
      real(8), intent( out) :: eval(:)
      real(8), optional, intent( out) :: evec(:,:)

      integer :: m, i, j, info, lwork
      logical :: rv
      integer, allocatable :: iwork(:), isuppz(:)
      real(8), allocatable :: cpy(:,:), work(:), diag(:), offd(:), tau(:), v(:,:), h(:,:)

      m = size( mat, 1)

      if( m .le. 0) then
        write(*,'("Error (rsydiag): Invalid matrix dimension.")')
        stop
      end if
      if( size( mat, 2) .ne. m) then
        write(*,'("Error (rsydiag): The matrix must be squared.")')
        stop
      end if
      if( size( eval, 1) .ne. m) then
        write(*,'("Error (rsydiag): The eigenvalue array must have the same dimension as the matrix.")')
        stop
      end if
      rv = .false.
      if( present( evec)) then
        rv = .true.
        if( (size( evec, 1) .ne. m) .or. (size( evec, 2) .ne. m)) then
          write(*,'("Error (rsydiag): The eigenvector array must have the same dimension as the matrix.")')
          stop
        end if
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

      eval = diag
      deallocate( work)
      allocate( work( max( 1, 2*m-2)))
      if( rv) then
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

        call dsteqr( 'v', m, eval, offd, evec, m, work, info)
      else
        call dsteqr( 'n', m, eval, offd, cpy, m, work, info)
      end if

      if( info .ne. 0) then
        write(*, '("Error( rsydiag): Diagonalisation failed. DSTEQR returned info", i4)') info
        stop
      end if

      deallocate( cpy, work, offd, tau, v, h)
      return
    end subroutine rsydiag

    !***********************************
    !   GENERALIZED DIAGONALIZATION
    !***********************************
    ! diagonalize complex general matrices
    subroutine zgegdiag( mat1, mat2, alpha, beta, levec, revec)
      complex(8), intent( in) :: mat1(:,:), mat2(:,:)
      complex(8), intent( out) :: alpha(:), beta(:)
      complex(8), optional, intent( out) :: levec(:,:), revec(:,:)

      integer :: m, info, lwork
      character :: lv, rv
      complex(8), allocatable :: cpy1(:,:), cpy2(:,:), work(:)
      real(8), allocatable :: rwork(:)

      m = size( mat1, 1)

      if( m .le. 0) then
        write(*,'("Error (zgegdiag): Invalid matrix dimension.")')
        stop
      end if
      if( size( mat1, 2) .ne. m) then
        write(*,'("Error (zgegdiag): The matrix must be squared.")')
        stop
      end if
      if( (size( mat2, 1) .ne. m) .or. (size( mat2, 2) .ne. m)) then
        write(*,'("Error (zgegdiag): Both matrices must have equal shapes.")')
        stop
      end if
      if( (size( alpha, 1) .ne. m) .or. (size( beta, 1) .ne. m)) then
        write(*,'("Error (zgegdiag): The eigenvalue arrays must have the same dimension as the matrices.")')
        stop
      end if
      lv = 'n'
      if( present( levec)) then
        lv = 'v'
        if( (size( levec, 1) .ne. m) .or. (size( levec, 2) .ne. m)) then
          write(*,'("Error (zgegdiag): The eigenvector array must have the same dimension as the matrices.")')
          stop
        end if
      end if
      rv = 'n'
      if( present( revec)) then
        rv = 'v'
        if( (size( revec, 1) .ne. m) .or. (size( revec, 2) .ne. m)) then
          write(*,'("Error (zgegdiag): The eigenvector array must have the same dimension as the matrices.")')
          stop
        end if
      end if

      allocate( cpy1( m, m), cpy2( m, m))
      allocate( work(1), rwork( 8*m))

      cpy1 = mat1
      cpy2 = mat2

      call zggev( lv, rv, m, cpy1, m, cpy2, m, alpha, beta, levec, m, revec, m, work, -1, rwork, info)
      lwork = work(1)
      deallocate( work)
      allocate( work( lwork))
      call zggev( lv, rv, m, cpy1, m, cpy2, m, alpha, beta, levec, m, revec, m, work, lwork, rwork, info)
        
      if( info .ne. 0) then
        write(*, '("Error( zgegdiag): Diagonalisation failed. ZGGEV returned info", i4)') info
        stop
      end if

      deallocate( cpy1, cpy2, work, rwork)
      return
    end subroutine zgegdiag

    ! diagonalize complex hermitian matrices
    ! mat2 must be positive definite
    subroutine zhegdiag( mat1, mat2, eval, evec, irange, erange)
      complex(8), intent( in) :: mat1(:,:), mat2(:,:)
      real(8), intent( out) :: eval(:)
      complex(8), optional, intent( out) :: evec(:,:)
      integer, optional, intent( in) :: irange(2)
      real(8), optional, intent( in) :: erange(2)

      integer :: m, n, i, j, il, iu, info, lwork
      character :: rv, rr
      real(8) :: phase, vl, vu
      integer, allocatable :: iwork(:), ifail(:)
      complex(8), allocatable :: cpy1(:,:), cpy2(:,:), work(:)
      real(8), allocatable :: rwork(:)

      m = size( mat1, 1)
      n = m

      if( m .le. 0) then
        write(*,'("Error (zhegdiag): Invalid matrix dimension.")')
        stop
      end if
      if( size( mat1, 2) .ne. m) then
        write(*,'("Error (zhegdiag): The matrix must be squared.")')
        stop
      end if
      if( (size( mat2, 1) .ne. m) .or. (size( mat2, 2) .ne. m)) then
        write(*,'("Error (zhegdiag): Both matrices must have equal shapes.")')
        stop
      end if
      rr = 'a'
      if( present( erange)) then
        vl = min( erange(1), erange(2))
        vu = max( erange(1), erange(2))
        rr = 'v'
      end if
      if( present( irange)) then
        il = max( 1, min( irange(1), irange(2)))
        iu = min( m, max( irange(1), irange(2)))
        n = iu - il + 1
        rr = 'i'
      end if
      if( size( eval, 1) .ne. n) then
        write(*,'("Error (zhegdiag): The eigenvalue array has incorrect dimensions.")')
        stop
      end if
      rv = 'n'
      if( present( evec)) then
        rv = 'v'
        if( (size( evec, 1) .ne. m) .or. (size( evec, 2) .ne. n)) then
          write(*,'("Error (zhegdiag): The eigenvector array has incorrect dimensions.")')
          stop
        end if
      end if

      allocate( cpy1, source=mat1)
      allocate( cpy2, source=mat2)
      allocate( work(1), rwork(7*m), iwork(5*m), ifail(m))

      call zhegvx( 1, rv, rr, 'u', m, cpy1, m, cpy2, m, vl, vu, il, iu, 2*DLAMCH('S'), i, eval, cpy1, m, work, -1, rwork, iwork, ifail, info)
      lwork = nint( dble( work(1)))
      deallocate( work)
      allocate( work( lwork))
      if( rv == 'v') then
        call zhegvx( 1, rv, rr, 'u', m, cpy1, m, cpy2, m, vl, vu, il, iu, 2*DLAMCH('S'), i, eval, evec, m, work, lwork, rwork, iwork, ifail, info)
      else
        call zhegvx( 1, rv, rr, 'u', m, cpy1, m, cpy2, m, vl, vu, il, iu, 2*DLAMCH('S'), i, eval, cpy1, m, work, lwork, rwork, iwork, ifail, info)
      end if

      if( info .ne. 0) then
        write(*, '("Error( zhegdiag): Diagonalisation failed. ZHEGVX returned info", i4)') info
        stop
      end if

      deallocate( cpy1, cpy2, work, rwork, iwork, ifail)
      return
    end subroutine zhegdiag

    !***********************************
    !           PSEUDO INVERSE
    !***********************************
    ! compute pseudo inverse of complex matrix
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
    
      call svd_divide_conquer( mat, sval, lsvec, rsvec)
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
    
      call svd_divide_conquer( mat, sval, lsvec, rsvec)
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

    !***********************************
    !            DETERMINANT
    !***********************************
    ! Determinant of complex matrix
    subroutine zdet( mat, det)
      complex(8), intent( in)  :: mat(:,:)
      complex(8), intent( out) :: det

      integer :: m, i, info
      integer, allocatable :: ipiv(:), jpiv(:)
      complex(8), allocatable :: cpy(:,:)

      m = size( mat, 1)

      if( m .le. 0) then
        write(*,'("Error (zdet): Invalid matrix dimension (",i6,")")') m
        stop
      end if
      if( size( mat, 2) .ne. m) then
        write(*,'("Error (zdet): The matrix must be squared.")')
        stop
      end if

      allocate( cpy, source=mat)
      allocate( ipiv( m), jpiv( m))

      call zgetc2( m, cpy, m, ipiv, jpiv, info)
      if( info .ne. 0) then
        write(*,'("Error (zdet): Determinant failed. ZGETC2 returned info ",i4)') info
        stop
      end if
      det = zone
      do i = 1, m
        det = det*cpy( i, i)
      end do
      do i = 1, m
        if( i .ne. ipiv( i)) det = -det
        if( i .ne. jpiv( i)) det = -det
      end do

      deallocate( ipiv, jpiv)
      return
    end subroutine zdet

    ! Determinant of real matrix
    subroutine rdet( mat, det)
      real(8), intent( in)  :: mat(:,:)
      real(8), intent( out) :: det

      integer :: m, i, info
      integer, allocatable :: ipiv(:), jpiv(:)
      real(8), allocatable :: cpy(:,:)

      m = size( mat, 1)

      if( m .le. 0) then
        write(*,'("Error (rdet): Invalid matrix dimension (",i6,")")') m
        stop
      end if
      if( size( mat, 2) .ne. m) then
        write(*,'("Error (rdet): The matrix must be squared.")')
        stop
      end if

      allocate( cpy, source=mat)
      allocate( ipiv( m), jpiv( m))

      call dgetc2( m, cpy, m, ipiv, jpiv, info)
      if( info .ne. 0) then
        write(*,'("Error (rdet): Determinant failed. DGETC2 returned info ",i4)') info
        stop
      end if
      det = 1.d0
      do i = 1, m
        det = det*cpy( i, i)
      end do
      do i = 1, m
        if( i .ne. ipiv( i)) det = -det
        if( i .ne. jpiv( i)) det = -det
      end do

      deallocate( ipiv, jpiv)
      return
    end subroutine rdet

    !***********************************
    !   LINEAR LEAST SQUARES PROBLEM
    !***********************************
    ! complex least squares problem
    subroutine zlsp( mat, rhs, sol)
      complex(8), intent( in)  :: mat(:,:), rhs(:,:)
      complex(8), intent( out) :: sol(:,:)

      integer :: m, n, k, info, rank, lwork
      integer, allocatable :: iwork(:)
      real(8), allocatable :: rwork(:)
      complex(8), allocatable :: cpy1(:,:), cpy2(:,:), s(:), work(:)

      m = size( mat, 1)
      n = size( mat, 2)
      k = size( rhs, 2)

      if( m .le. 0) then
        write(*,'("Error (zlsp): Invalid matrix dimension (m=",i6,")")') m
        stop
      end if
      if( n .le. 0) then
        write(*,'("Error (zlsp): Invalid matrix dimension (n=",i6,")")') n
        stop
      end if
      if( k .le. 0) then
        write(*,'("Error (zlsp): Invalid matrix dimension (k=",i6,")")') k
        stop
      end if
      if( size( rhs, 1) .ne. m) then
        write(*,'("Error (zlsp): The right hand side matrix must have as many rows as the coefficient matrix.")')
        stop
      end if
      if( size( sol, 1) .ne. n) then
        write(*,'("Error (zlsp): The solution matrix must have as many rows as the coefficient matrix has columns.")')
        stop
      end if
      if( size( sol, 2) .ne. k) then
        write(*,'("Error (zlsp): The solution matrix must have as many columns as the right hand side matrix.")')
        stop
      end if

      allocate( cpy1( m, n), cpy2( m, k), s( min( m, n)))
      cpy1 = mat
      cpy2 = rhs

      allocate( work( 1), iwork( 1), rwork( 1))
      call zgelsd( m, n, k, cpy1, m, cpy2, m, s, -1.d0, rank, work, -1, rwork, iwork, info)
      lwork = work(1); deallocate( work); allocate( work( lwork))
      info = iwork(1); deallocate( iwork); allocate( iwork( info))
      info = nint( rwork(1)); deallocate( rwork); allocate( rwork( info))
      call zgelsd( m, n, k, cpy1, m, cpy2, m, s, -1.d0, rank, work, lwork, rwork, iwork, info)
      if( info .ne. 0) then
        write(*,'("Error (zlsp): Least squares problem failed. ZGELSD returned info ",i4)') info
        stop
      end if
      sol = cpy2( 1:n, :)
      deallocate( cpy1, cpy2, s, work, iwork, rwork)
      return
    end subroutine zlsp

    ! real least squares problem
    subroutine rlsp( mat, rhs, sol)
      real(8), intent( in)  :: mat(:,:), rhs(:,:)
      real(8), intent( out) :: sol(:,:)

      integer :: m, n, k, info, rank, lwork
      integer, allocatable :: iwork(:)
      real(8), allocatable :: cpy1(:,:), cpy2(:,:), s(:), work(:)

      m = size( mat, 1)
      n = size( mat, 2)
      k = size( rhs, 2)

      if( m .le. 0) then
        write(*,'("Error (rlsp): Invalid matrix dimension (m=",i6,")")') m
        stop
      end if
      if( n .le. 0) then
        write(*,'("Error (rlsp): Invalid matrix dimension (n=",i6,")")') n
        stop
      end if
      if( k .le. 0) then
        write(*,'("Error (rlsp): Invalid matrix dimension (k=",i6,")")') k
        stop
      end if
      if( size( rhs, 1) .ne. m) then
        write(*,'("Error (rlsp): The right hand side matrix must have as many rows as the coefficient matrix.")')
        stop
      end if
      if( size( sol, 1) .ne. n) then
        write(*,'("Error (rlsp): The solution matrix must have as many rows as the coefficient matrix has columns.")')
        stop
      end if
      if( size( sol, 2) .ne. k) then
        write(*,'("Error (rlsp): The solution matrix must have as many columns as the right hand side matrix.")')
        stop
      end if

      allocate( cpy1( m, n), cpy2( m, k), s( min( m, n)))
      cpy1 = mat
      cpy2 = rhs

      allocate( work( 1), iwork( 1))
      call dgelsd( m, n, k, cpy1, m, cpy2, m, s, -1.d0, rank, work, -1, iwork, info)
      lwork = work(1)
      info = iwork(1)
      deallocate( work, iwork)
      allocate( work( lwork), iwork( info))
      call dgelsd( m, n, k, cpy1, m, cpy2, m, s, -1.d0, rank, work, lwork, iwork, info)
      if( info .ne. 0) then
        write(*,'("Error (rlsp): Least squares problem failed. DGELSD returned info ",i4)') info
        stop
      end if
      sol = cpy2( 1:n, :)
      deallocate( cpy1, cpy2, s, work, iwork)
      return
    end subroutine rlsp

    !***********************************
    !         QR FACTORIZATION
    !***********************************
    ! complex QR factorization
    subroutine zqr( mat, q, r)
      complex(8), intent( in) :: mat(:,:)
      complex(8), optional, intent( out) :: q(:,:)
      complex(8), optional, intent( out) :: r(:,:)

      integer :: m, n, lwork, info, i, j, k, l
      complex(8) :: v1, v2
      logical :: doq, dor
      complex(8), allocatable :: cpy(:,:), tau(:), work(:)

      doq = .false.
      if( present( q)) doq = .true.
      dor = .false.
      if( present( r)) dor = .true.

      m = size( mat, 1)
      n = size( mat, 2)
      k = min( m, n)

      if( m .le. 0) then
        write(*,'("Error (zqr): Invalid matrix dimension (m=",i6,")")') m
        stop
      end if
      if( n .le. 0) then
        write(*,'("Error (zqr): Invalid matrix dimension (n=",i6,")")') n
        stop
      end if
      if( doq) then
        if( size( q, 1) .ne. m) then
          write(*,'("Error (zqr): The matrix Q must have as many rows as the input matrix has.")') n
          stop
        end if
        if( size( q, 2) .ne. n) then
          write(*,'("Error (zqr): The matrix Q must have as many columns as the input matrix has.")') n
          stop
        end if
      end if
      if( dor) then
        if( size( r, 1) .ne. m) then
          write(*,'("Error (zqr): The matrix R must have as many rows as the input matrix has.")') n
          stop
        end if
        if( size( r, 2) .ne. n) then
          write(*,'("Error (zqr): The matrix R must have as many columns as the input matrix has.")') n
          stop
        end if
      end if

      allocate( cpy, source=mat)

      allocate( tau(k), work(1))
      call zgeqrfp( m, n, cpy, m, tau, work, -1, info)
      lwork = nint( dble( work(1)))
      deallocate( work)
      allocate( work( lwork))
      call zgeqrfp( m, n, cpy, m, tau, work, lwork, info)
      if( info .ne. 0) then
        write(*,'("Error (zqr): QR factorization failed. ZGEQRFP returned info ",i4)') info
        stop
      end if

      if( dor) then
        r = zzero
        do i = 1, k
          do j = i, n
            r(i,j) = cpy(i,j)
          end do
        end do
      end if

      if( doq) then
        call zungqr( m, k, k, cpy, m, tau, work, -1, info)
        lwork = nint( dble( work(1)))
        deallocate( work)
        allocate( work( lwork))
        call zungqr( m, k, k, cpy, m, tau, work, lwork, info)
        if( info .ne. 0) then
          write(*,'("Error (zqr): QR factorization failed. ZUNGQR returned info ",i4)') info
          stop
        end if
        q = cpy
      end if
      deallocate( cpy, tau, work)
      return
    end subroutine zqr

    ! real QR factorization
    subroutine rqr( mat, q, r)
      real(8), intent( in)  :: mat(:,:)
      real(8), optional, intent( out) :: q(:,:)
      real(8), optional, intent( out) :: r(:,:)

      integer :: m, n, lwork, info, i, j, k, l
      real(8) :: v1, v2
      logical :: doq, dor
      real(8), allocatable :: cpy(:,:), tau(:), work(:)

      doq = .false.
      if( present( q)) doq = .true.
      dor = .false.
      if( present( r)) dor = .true.

      m = size( mat, 1)
      n = size( mat, 2)
      k = min( m, n)

      if( m .le. 0) then
        write(*,'("Error (rqr): Invalid matrix dimension (m=",i6,")")') m
        stop
      end if
      if( n .le. 0) then
        write(*,'("Error (rqr): Invalid matrix dimension (n=",i6,")")') n
        stop
      end if
      if( doq) then
        if( size( q, 1) .ne. m) then
          write(*,'("Error (rqr): The matrix Q must have as many rows as the input matrix has.")') n
          stop
        end if
        if( size( q, 2) .ne. n) then
          write(*,'("Error (rqr): The matrix Q must have as many columns as the input matrix has.")') n
          stop
        end if
      end if
      if( dor) then
        if( size( r, 1) .ne. m) then
          write(*,'("Error (rqr): The matrix R must have as many rows as the input matrix has.")') n
          stop
        end if
        if( size( r, 2) .ne. n) then
          write(*,'("Error (rqr): The matrix R must have as many columns as the input matrix has.")') n
          stop
        end if
      end if

      allocate( cpy, source=mat)

      allocate( tau(k), work(1))
      call dgeqrfp( m, n, cpy, m, tau, work, -1, info)
      lwork = work(1)
      deallocate( work)
      allocate( work( lwork))
      call dgeqrfp( m, n, cpy, m, tau, work, lwork, info)
      if( info .ne. 0) then
        write(*,'("Error (rqr): QR factorization failed. DGEQRFP returned info ",i4)') info
        stop
      end if

      if( dor) then
        r = 0.d0
        do i = 1, k
          do j = i, n
            r(i,j) = cpy(i,j)
          end do
        end do
      end if

      if( doq) then
        call dorgqr( m, k, k, cpy, m, tau, work, -1, info)
        lwork = nint( dble( work(1)))
        deallocate( work)
        allocate( work( lwork))
        call dorgqr( m, k, k, cpy, m, tau, work, lwork, info)
        if( info .ne. 0) then
          write(*,'("Error (rqr): QR factorization failed. DORGQR returned info ",i4)') info
          stop
        end if
        q = cpy
      end if
      deallocate( cpy, tau, work)
      return
    end subroutine rqr

    !***********************************
    !         MATRIX EXPONENTIAL
    !***********************************
    ! exponential of a complex matrix
    subroutine zexpm( mat, expm)
      complex(8), intent( in)  :: mat(:,:)
      complex(8), intent( out) :: expm(:,:)

      integer :: m, s, q, i
      real(8) :: norm, c
      logical :: swap

      integer, allocatable :: ipiv(:)
      complex(8), allocatable :: d(:,:), cpy(:,:), tmp(:,:), tmp2(:,:)

      m = size( mat, 1)

      if( m .le. 0) then
        write(*,'("Error (zexpm): Invalid matrix dimension.")')
        stop
      end if
      if( size( mat, 2) .ne. m) then
        write(*,'("Error (zexpm): The matrix must be squared.")')
        stop
      end if
      if( (size( expm, 1) .ne. m) .or. (size( expm, 2) .ne. m)) then
        write(*,'("Error (zexpm): The exponential must have the same dimension as the matrix.")')
        stop
      end if

      allocate( cpy(m,m), d(m,m), tmp(m,m), tmp2(m,m), ipiv(m))
      norm = sum( abs( mat))
      if( norm .lt. 1.d-16) then
        expm = zzero
        do i = 1, m
          expm(i,i) = zone
        end do
        return
      end if
      s = max( 0, ceiling( log(norm)/log(2.d0)))
      !write(*,'(2i,f13.6)') m, s, norm

      c = 0.5d0; q = 6; swap = .true.
      cpy = mat/2.d0**s
      tmp = cpy
      expm = c*cpy
      d = -c*cpy
      do i = 1, m
        expm(i,i) = zone + expm(i,i)
        d(i,i) = zone + d(i,i)
      end do
      do i = 2, q
        c = c*dble(q-i+1)/dble(i*(2*q-i+1))
        call zgemm( 'n', 'n', m, m, m, zone, cpy, m, tmp, m, zzero, tmp2, m)
        tmp = tmp2
        expm = expm + c*tmp
        if( swap) then
          d = d + c*tmp
        else
          d = d - c*tmp
        end if
        swap = .not. swap
      end do
      call zgesv( m, m, d, m, ipiv, expm, m, i)
      do i = 1, s
        call zgemm( 'n', 'n', m, m, m, zone, expm, m, expm, m, zzero, tmp, m)
        expm = tmp
      end do
      deallocate( d, cpy, tmp, tmp2, ipiv)
      return
    end subroutine zexpm

    ! exponential of a real matrix
    subroutine rexpm( mat, expm)
      real(8), intent( in)  :: mat(:,:)
      real(8), intent( out) :: expm(:,:)

      integer :: m, s, q, i
      real(8) :: norm, c
      logical :: swap

      integer, allocatable :: ipiv(:)
      real(8), allocatable :: d(:,:), cpy(:,:), tmp(:,:), tmp2(:,:)

      m = size( mat, 1)

      if( m .le. 0) then
        write(*,'("Error (rexpm): Invalid matrix dimension.")')
        stop
      end if
      if( size( mat, 2) .ne. m) then
        write(*,'("Error (rexpm): The matrix must be squared.")')
        stop
      end if
      if( (size( expm, 1) .ne. m) .or. (size( expm, 2) .ne. m)) then
        write(*,'("Error (rexpm): The exponential must have the same dimension as the matrix.")')
        stop
      end if

      allocate( cpy(m,m), d(m,m), tmp(m,m), tmp2(m,m), ipiv(m))
      norm = sum( abs( mat))
      if( norm .lt. 1.d-16) then
        expm = 0.d0
        do i = 1, m
          expm(i,i) = 1.d0
        end do
        return
      end if
      s = max( 0, ceiling( log(norm)/log(2.d0)))
      !write(*,'(2i,f13.6)') m, s, norm

      c = 0.5d0; q = 6; swap = .true.
      cpy = mat/2.d0**s
      tmp = cpy
      expm = c*cpy
      d = -c*cpy
      do i = 1, m
        expm(i,i) = 1.d0 + expm(i,i)
        d(i,i) = 1.d0 + d(i,i)
      end do
      do i = 2, q
        c = c*dble(q-i+1)/dble(i*(2*q-i+1))
        call dgemm( 'n', 'n', m, m, m, 1.d0, cpy, m, tmp, m, 0.d0, tmp2, m)
        tmp = tmp2
        expm = expm + c*tmp
        if( swap) then
          d = d + c*tmp
        else
          d = d - c*tmp
        end if
        swap = .not. swap
      end do
      call dgesv( m, m, d, m, ipiv, expm, m, i)
      do i = 1, s
        call dgemm( 'n', 'n', m, m, m, 1.d0, expm, m, expm, 0.d0, tmp, m)
        expm = tmp
      end do
      deallocate( d, cpy, tmp, tmp2, ipiv)
      return
    end subroutine rexpm

    !***********************************
    !   UNITARY/ORTHOGONAL COMPLETION
    !***********************************
    ! complete a set of orthonormal column vectors
    subroutine zucomp( mat, com)
      complex(8), intent( in)  :: mat(:,:)
      complex(8), intent( out) :: com(:,:)

      integer :: m, n, i, j
      real(8), allocatable :: val(:)
      complex(8), allocatable :: tmp(:,:), q(:,:), r(:,:)

      m = size( mat, 1)
      n = size( mat, 2)

      if( m .lt. n) then
        write(*,'("Error (zucomp): The matrix must not have more columns than rows.")')
        stop
      end if
      if( (size( com, 1) .ne. m) .or. (size( com, 2) .ne. m)) then
        write(*,'("Error (zucomp): The completed matrix must be squared with the same number of columns as the input matrix.")')
        stop
      end if

      if( m .eq. n) then
        com = mat
        return
      end if

      allocate( tmp(m,m), val(m), q(m,m), r(m,m))
      call zgemm( 'n', 'c', m, m, n, -zone, mat, m, mat, m, zzero, tmp, m)
      do i = 1, m
        tmp(i,i) = tmp(i,i) + zone
      end do
      com(:,1:n) = mat
      call zqr( tmp, q=q, r=r)
      do i = 1, m
        val(i) = abs( dble( r(i,i)))
      end do
      do i = n+1, m
        j = maxloc( val, 1)
        com(:,i) = q(:,j)
        val(j) = 0.d0
      end do
      deallocate( tmp, val, q, r)

      return
    end subroutine zucomp
end module m_linalg
