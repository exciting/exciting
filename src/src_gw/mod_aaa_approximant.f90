
module mod_aaa_approximant

    implicit none

    type aaa_approximant
        integer(4) :: nj                 ! number of support points
        complex(8), allocatable :: zj(:) ! support point coordinates
        complex(8), allocatable :: fj(:) ! function values
        complex(8), allocatable :: wj(:) ! barycentric weights
        real(8),    allocatable :: ej(:) ! error at support points
    end type

    complex(8), parameter :: zero = cmplx(0.d0, 0.d0, 8)

contains

!--------------------------------------------------------------------------------    
    subroutine test_aaa()
        implicit none
        integer(4) :: n, j
        real(8)    :: a, b, d
        complex(8), allocatable :: z(:), f(:), fr(:)
        type(aaa_approximant) :: aaa

        ! Define support points
        n = 6
        allocate(z(n), f(n))

        a = -1.d0
        b =  1.d0
        d = (b-a) / dble(n-1)
        do j = 1, n
            z(j) = cmplx(a + dble(j-1)*d, 0.d0, 8)
            f(j) = exp(z(j))
            ! print*, j, zj(j), fj(j)
        end do

        ! Compute the approximant
        call set_aaa_approximant(aaa, z, f)
        print*, ''
        print*, 'Number of support points: ', aaa%nj
        print*, 'Error at support points: ', aaa%ej(:)
        print*, ''
        deallocate(z, f)

        ! Compute values of the barycentric rational function
        n = 100
        allocate(z(n), f(n), fr(n))
        a = -3.d0
        b =  3.d0
        d = (b-a) / dble(n-1)
        open(77, file="aaa_approximant.dat")
        do j = 1, n
            z(j) = cmplx(a + dble(j-1)*d, 0.d0, 8)
            f(j) = exp(z(j))
            fr(j) = reval_aaa_approximant(aaa, z(j))
            write(77,'(3f18.6)') dble(z(j)), dble(f(j)), dble(fr(j))
        end do
        close(77)
        deallocate(z, f, fr)
        call delete_aaa_approximant(aaa)

    end subroutine

!--------------------------------------------------------------------------------    
    function reval_aaa_approximant(aaa, z) result(zr)
        implicit none
        type(aaa_approximant), intent(in) :: aaa
        complex(8), intent(in) :: z
        complex(8)             :: zr
        ! local
        real(8), parameter     :: epstol = 1.d-16
        logical    :: ldone
        integer(4) :: j
        complex(8) :: zdiff, N, D
        complex(8), allocatable :: C(:)

        ! check if z coincides with one of the supporting points
        ldone = .false.
        do j = 1, aaa%nj
            zdiff = z-aaa%zj(j)
            if (abs(zdiff) < epstol) then
                zr = aaa%fj(j)
                ldone = .true.
                exit
            end if
        end do

        if (.not.ldone) then
            ! Compute rational function
            ! Build Cauchy matrix
            allocate(C(aaa%nj))
            do j = 1, aaa%nj
                C(j) = 1.d0 / (z - aaa%zj(j))
            end do
            N = dot_product(C(:), aaa%wj(:)*aaa%fj(:))
            D = dot_product(C(:), aaa%wj(:))
            zr = N / D
            deallocate(C)
        end if

        return
    end function


!--------------------------------------------------------------------------------    
    subroutine delete_aaa_approximant(aaa)
        type(aaa_approximant), intent(inout) :: aaa
        if (allocated(aaa%zj)) deallocate(aaa%zj)
        if (allocated(aaa%fj)) deallocate(aaa%fj)
        if (allocated(aaa%wj)) deallocate(aaa%wj)
    end subroutine


!--------------------------------------------------------------------------------    
    subroutine init_aaa_approximant(aaa, n, z, f, w, e)
        type(aaa_approximant), intent(out) :: aaa
        integer(4), intent(in) :: n
        complex(8), intent(in) :: z(:)
        complex(8), intent(in) :: f(:)
        complex(8), intent(in) :: w(:)
        real(8),    intent(in) :: e(:)
        integer(4) :: j
        aaa%nj = n
        if (allocated(aaa%zj)) deallocate(aaa%zj)
        allocate(aaa%zj(aaa%nj))
        if (allocated(aaa%fj)) deallocate(aaa%fj)
        allocate(aaa%fj(aaa%nj))
        if (allocated(aaa%wj)) deallocate(aaa%wj)
        allocate(aaa%wj(aaa%nj))
        if (allocated(aaa%ej)) deallocate(aaa%ej)
        allocate(aaa%ej(aaa%nj))
        do j = 1, aaa%nj
            aaa%zj(j) = z(j)
            aaa%fj(j) = f(j)
            aaa%wj(j) = w(j)
            aaa%ej(j) = e(j)
        end do
    end subroutine


!--------------------------------------------------------------------------------    
    subroutine set_aaa_approximant(aaa, Z, F)
        implicit none
        ! input/output
        type(aaa_approximant), intent(out) :: aaa
        complex(8),            intent(in)  :: Z(:)
        complex(8),            intent(in)  :: F(:)
        ! parameters
        integer(4), parameter :: mmax = 1000
        real(8),    parameter :: tol = 1.d-8
        
        ! local variables
        integer(4) :: npts, j, jj, jmax(1), m, mm
        real(8)    :: reltol
        integer(4), allocatable :: jndx(:)
        real(8),    allocatable :: error(:)
        complex(8), allocatable :: RS(:,:), LS(:,:)
        complex(8), allocatable :: A(:,:), B(:,:), C(:,:)
        complex(8), allocatable :: R(:), N(:), D(:)
        complex(8), allocatable :: zj(:), fj(:), wj(:)
        ! SVD output
        real(8),    allocatable :: sval(:)
        complex(8), allocatable :: rsvec(:,:)

        ! input data size
        npts = size(Z)
        ! print*, 'npts=', npts

        ! relative tolerance <---- check it later
        reltol = tol ! * znorm(F)
        ! print*, 'reltol=', reltol

        ! Left scaling matrix
        allocate(LS(npts,npts))
        LS(:,:) = zero
        do j = 1, npts
            LS(j,j) = F(j)
        end do

        ! Right scaling matrix
        allocate(RS(npts,npts))
        RS(:,:) = zero

        ! Error at sampling points
        allocate(error(mmax))

        ! SVD arrays
        allocate(sval(mmax))
        allocate(rsvec(mmax,mmax))

        ! Initialization for AAA iteration
        allocate(zj(mmax))
        zj(:) = zero
        allocate(fj(mmax))
        fj(:) = zero
        allocate(wj(mmax))
        wj(:) = zero

        allocate(A(npts,mmax))
        allocate(B(npts,mmax))
        allocate(C(npts,mmax))

        ! mean value
        allocate(R(npts), N(npts), D(npts))
        R(:) = sum(F) / dble(npts)

        ! Index vector
        allocate(jndx(npts))
        do j = 1, npts
            jndx(j) = j
        end do

        do m = 1, mmax

            ! Case when all points are used but convergence is not reached
            if ( m == npts ) then
                stop 'ERROR(mod_aaa_approximant) Convergence cannot be reached! Try to increase the input data size.'
            end if

            ! Select next supporting point where error is largest
            jmax = maxloc(abs(F-R))
            ! print*, 'jmax=', jmax

            ! Update index vector
            jndx(jmax(1)) = 0  ! exclude the point

            zj(m) = Z(jmax(1)) ! Update support points
            fj(m) = F(jmax(1)) ! Update data values

            RS(m,m) = fj(m)    ! Next element of right scaling matrix

            ! Next column of Cauchy matrix
            do j = 1, npts
                if ( jmax(1) /= j ) then
                    C(j,m) = 1.d0 / (Z(j) - Z(jmax(1)))
                else
                    C(j,m) = zero
                end if
            end do
            ! print*, 'C=', sum(C(:,1:m))

            !-------------------
            ! Compute weights
            !-------------------

            ! Loewner matrix
            A(1:npts,1:m) = matmul(LS(1:npts,1:npts), C(1:npts,1:m)) - &
                            matmul(C(1:npts,1:m), RS(1:m,1:m))
           
            ! SVD of the reduced matrix
            jj = 0
            do j = 1, npts
                if (jndx(j) > 0) then
                    jj = jj+1
                    B(jj,1:m) = A(j,1:m)
                    ! print*, 'B=', dble(B(jj,1:m))
                end if
            end do
            ! print*, 'jj=', jj

            call mkl_svd(jj, m, B(1:jj,1:m), rsvec(1:m,1:m), sval(1:m))

            ! Weight vector = min singular vector
            wj(1:m) = rsvec(1:m,m)

            ! Rational approximant on Z
            N(1:npts) = matmul(C(1:npts,1:m), wj(1:m)*fj(1:m))
            D(1:npts) = matmul(C(1:npts,1:m), wj(1:m))

            do j = 1, npts
                if (jndx(j) > 0) then
                    R(j) = N(j) / D(j)
                else
                    R(j) = F(j)
                end if
            end do

            ! Error in the sample points
            error(m) = maxval(abs(F-R))

            ! Check if converged
            if (error(m) <= reltol) exit

        end do ! m
        deallocate(sval, rsvec)
        deallocate(A, B, C)
        deallocate(RS, LS)

        ! Save results into an object
        call init_aaa_approximant(aaa, m, zj, fj, wj, error)
        
        deallocate(zj, fj, wj, error)

    end subroutine


!--------------------------------------------------------------------------------
    real(8) function znorm(z)
        implicit none
        complex(8), intent(in) :: z(:)
        integer(4) :: i, n
        n = size(z)
        znorm = 0.d0
        do i = 1, n
            znorm = znorm + abs(z(i))**2
        end do
        znorm = sqrt(znorm)
    end function
  

!--------------------------------------------------------------------------------
    subroutine mkl_svd(m, n, A, evec, eval)
        implicit none
        integer,    intent(in)  :: m
        integer,    intent(in)  :: n
        complex(8), intent(in)  :: A(m,n)
        complex(8), intent(out) :: evec(n,n)
        real(8),    intent(out) :: eval(n)

        ! local
        ! RWORK dimension should be at least MAX( 1, 5*MIN(M,N) )
        integer :: i
        integer :: lmn, lwork, lrwork, info
        real(8),    allocatable :: S(:), rwork(:)
        complex(8), allocatable :: U(:,:), VT(:,:), work(:)
        complex(8), allocatable :: A_(:,:)
        external zgesvd

        allocate(A_(m,n))
        A_(:,:) = A(:,:)

        lmn = min(m,n)
        allocate(S(lmn))
        allocate(U(m,m),VT(n,n))
        lrwork = 5*lmn
        allocate(rwork(lrwork))
        !
        ! Query the optimal workspace
        !
        lwork = -1
        allocate(work(1))
        call zgesvd( 'all', 'all', m, n, A_, m, S, U, m, VT, n, &
        &             work, lwork, rwork, info )
        lwork  = int(work(1))
        deallocate(work)
        !
        ! Compute SVD
        !
        allocate(work(lwork))
        call zgesvd( 'all', 'all', m, n, A_, m, S, U, m, VT, n, &
        &             work, lwork, rwork, info )
        !
        ! check for convergence.
        !
        if( info > 0 ) then
        write(*,*)'The algorithm computing svd failed to converge.'
        stop
        end if
        !
        ! Eigenvalues: \lambda = S^{+}*S
        !
        eval(:) = 0.d0
        do i = 1, m
        eval(i) = S(i)*S(i)
        end do
        !
        ! Eigenvectors
        !
        do i = 1, n
        evec(:,i) = conjg(VT(i,:))
        end do

        deallocate(A_)
        deallocate(S)
        deallocate(U,VT)
        deallocate(rwork,work)

        return
    end subroutine

end module