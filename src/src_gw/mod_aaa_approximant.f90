
module mod_aaa_approximant

    implicit none

    type aaa_approximant
        integer(4) :: nj                  ! number of support points
        complex(8), allocatable :: zj(:)  ! support point coordinates
        complex(8), allocatable :: fj(:)  ! function values
        complex(8), allocatable :: wj(:)  ! barycentric weights
        integer(4) :: npol                ! number of poles
        complex(8), allocatable :: pol(:) ! poles
        complex(8), allocatable :: res(:) ! residues
        complex(8), allocatable :: zer(:) ! zeros of rational function in barycentric form
    end type

    complex(8), parameter :: zero = cmplx(0.d0, 0.d0, 8)
    complex(8), parameter :: zone = cmplx(1.d0, 0.d0, 8)

contains

!--------------------------------------------------------------------------------    
    function linspace(a, b, n)
        implicit none
        complex(8), intent(in) :: a
        complex(8), intent(in) :: b
        integer(4), intent(in) :: n
        complex(8), dimension(n) :: linspace
        integer(4) :: i
        complex(8) :: step
        if (n < 1) then
            stop 'Error(linspace): Wrong number of points!'
        else if (n == 1) then
            linspace(1) = a
        else
            step = (b-a) / dble(n-1)
            do i = 1, n
                linspace(i) = a + dble(i-1)*step
            end do
        end if
    end function

!--------------------------------------------------------------------------------    
    subroutine test_aaa_2()
        use modinput
        implicit none
        integer(4) :: n
        complex(8) :: a, b
        complex(8), allocatable :: z(:), f(:)
        real(8), parameter :: pi = 3.14159265359d0
        integer(4) :: j
        real(8) :: tol
        type(aaa_approximant) :: aaa

        print*, '-----------------------------------------'
        print*, 'Test AAA 2: Poles of function tan(pi*z/2)'
        print*, '-----------------------------------------'

        a = cmplx(-0.5d0, 0.d0, 8)
        b = cmplx(0.5d0, 15d0*pi, 8)
        n = 1000

        allocate(z(n), f(n))
        z = linspace(a, b, n)
        do j = 1, n
            z(j) = exp(z(j))
            f(j) = tan(0.5d0*pi*z(j))
            ! print*, z(i), f(i)
        end do

        ! Compute the approximant
        tol = input%gw%selfenergy%tol
        call set_aaa_approximant(aaa, z, f, tol, .true.)
        print*, ''
        print*, 'Number of support points: ', aaa%nj
        print*, ''
        deallocate(z, f)

        call delete_aaa_approximant(aaa)
        
    end subroutine

!--------------------------------------------------------------------------------    
    subroutine test_aaa_1()
        use modinput
        implicit none
        integer(4) :: n, j
        real(8)    :: tol
        complex(8) :: a, b
        complex(8), allocatable :: z(:), f(:), fr(:)
        type(aaa_approximant) :: aaa

        print*, '------------------------------------------'
        print*, 'Test AAA 1: Rational approximant of exp(z)'
        print*, '------------------------------------------'

        ! Define support points
        a = cmplx(-1.d0, 0.d0, 8)
        b = cmplx( 1.d0, 0.d0, 8)
        n = 100

        allocate(z(n), f(n))
        z = linspace(a, b, n)
        do j = 1, n
            f(j) = exp(z(j))
            ! print*, j, zj(j), fj(j)
        end do

        ! Compute the approximant
        tol = input%gw%selfenergy%tol
        call set_aaa_approximant(aaa, z, f, tol, .true.)
        print*, ''
        print*, 'Number of support points: ', aaa%nj
        print*, ''
        deallocate(z, f)

        ! Compute values of the barycentric rational function
        a = cmplx(-3.d0, 0.d0, 8)
        b = cmplx( 3.d0, 0.d0, 8)
        n = 100
        allocate(z(n), f(n), fr(n))
        open(77, file="aaa_approximant.dat")
        do j = 1, n
            f(j) = exp(z(j))
            fr(j) = get_aaa_approximant(aaa, z(j))
            write(77,'(3f18.6)') dble(z(j)), dble(f(j)), dble(fr(j))
        end do
        close(77)
        deallocate(z, f, fr)
        call delete_aaa_approximant(aaa)

    end subroutine

!--------------------------------------------------------------------------------    
    function get_aaa_approximant(aaa, z) result(fz)
        implicit none
        type(aaa_approximant), intent(in) :: aaa
        complex(8), intent(in) :: z
        complex(8)             :: fz
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
                fz = aaa%fj(j)
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
            fz = N / D
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
        if (allocated(aaa%pol)) deallocate(aaa%pol)
        if (allocated(aaa%res)) deallocate(aaa%res)
        if (allocated(aaa%zer)) deallocate(aaa%zer)
    end subroutine


!--------------------------------------------------------------------------------    
    subroutine init_aaa_approximant(aaa, n, z, f, w)
        type(aaa_approximant), intent(out) :: aaa
        integer(4), intent(in) :: n
        complex(8), intent(in) :: z(:)
        complex(8), intent(in) :: f(:)
        complex(8), intent(in) :: w(:)
        integer(4) :: j
        aaa%nj = n
        if (allocated(aaa%zj)) deallocate(aaa%zj)
        allocate(aaa%zj(aaa%nj))
        if (allocated(aaa%fj)) deallocate(aaa%fj)
        allocate(aaa%fj(aaa%nj))
        if (allocated(aaa%wj)) deallocate(aaa%wj)
        allocate(aaa%wj(aaa%nj))
        do j = 1, aaa%nj
            aaa%zj(j) = z(j)
            aaa%fj(j) = f(j)
            aaa%wj(j) = w(j)
        end do
    end subroutine


!--------------------------------------------------------------------------------    
    subroutine set_aaa_approximant(aaa, Z, F, tol, cleanup_flag, cleanup_tol)
        implicit none
        ! input/output
        type(aaa_approximant), intent(out) :: aaa
        complex(8), intent(in)             :: Z(:)
        complex(8), intent(in)             :: F(:)
        real(8),    intent(in), optional   :: tol
        logical,    intent(in), optional   :: cleanup_flag
        real(8),    intent(in), optional   :: cleanup_tol
        ! parameters
        integer(4),   parameter :: mmax = 1000
        
        ! local variables
        integer(4) :: npts, j, jj, jmax(1), m
        real(8)    :: reltol, error
        integer(4), allocatable :: jndx(:)
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

        ! relative tolerance
        if (present(tol)) then
            reltol = tol
        else
            reltol = 1.d-14
        end if
        reltol = reltol * znorm(F)
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
                stop 'ERROR(mod_aaa_approximant) Convergence cannot be reached! Increase the tolerance or the sampling data size.'
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
            error = maxval(abs(F-R))
            ! write(*,*) 'iter=', m, 'error=', error

            ! Check if converged
            if (error <= reltol) exit

        end do ! m
        deallocate(sval, rsvec)
        deallocate(A, B, C)
        deallocate(RS, LS)
        deallocate(jndx)

        ! Save results into an object
        call init_aaa_approximant(aaa, m, zj, fj, wj)
        deallocate(zj, fj, wj)

        ! Remove spurious pole-zero pairs
        if (present(cleanup_flag)) then
            if (cleanup_flag) then
                if (present(cleanup_tol)) then
                    reltol = cleanup_tol
                else
                    reltol = 1.d-13
                end if
                reltol = reltol * znorm(F)
                print*, 'FD tol=', reltol
                call prz(aaa)
                call cleanup(aaa, Z, F, reltol)
            end if
        end if
        
    end subroutine


!--------------------------------------------------------------------------------
    subroutine prz(aaa)
        implicit none
        type(aaa_approximant), intent(InOut) :: aaa
        ! local
        integer(4) :: m, i, j
        integer(4), allocatable :: idxPol(:)
        
        complex(8), allocatable :: B(:,:), E(:,:)
        complex(8), allocatable :: D(:), N(:)
        complex(8) :: DD, NN

        integer(4) :: ndim
        complex(8), allocatable :: alpha(:), beta(:)

        m = aaa%nj

        ! Compute poles via generalized eigenvalue problem
        allocate(B(m+1,m+1), E(m+1,m+1))
        B(:,:) = zero
        E(:,:) = zero

        E(1,1) = zero
        do i = 2, m+1
            B(i,i) = zone
            E(1,i) = aaa%wj(i-1)
            E(i,1) = zone
            E(i,i) = aaa%zj(i-1)
        end do
        ! do i = 1, m+1
        !     write(*,'(100f12.6)') dble(E(i,:))
        ! end do

        ! Generalized eigenvalue problem
        ndim = m+1
        allocate(alpha(ndim), beta(ndim))
        call mkl_zggev(ndim, E, B, alpha, beta)

        ! Remove zeros of denominator at infinity
        allocate(idxPol(ndim))
        aaa%nPol = 0
        write(*,*)
        do i = 1, ndim
            write(*,*) 'pole=', alpha(i)/beta(i)
            if (abs(beta(i)) > 1.d-16) then
                aaa%nPol = aaa%nPol+1
                idxPol(aaa%nPol) = i
            end if
        end do
        write(*,*)
        write(*,*) "nPol=", aaa%nPol
        write(*,*)
        if (aaa%nPol < 1) stop "Error(prz): No poles found!"

        if (allocated(aaa%pol)) deallocate(aaa%pol)
        allocate(aaa%pol(aaa%nPol))
        do i = 1, aaa%nPol
            j = idxPol(i)
            aaa%pol(i) = alpha(j) / beta(j)
        end do
        deallocate(idxPol)

        ! Compute residues via formula for res of quotient of analytic functions
        if (allocated(aaa%res)) deallocate(aaa%res)
        allocate(aaa%res(aaa%nPol))
        allocate(N(aaa%nj), D(aaa%nj))
        do i = 1, aaa%nPol
            do j = 1, aaa%nj
                N(j) = 1.d0 / (aaa%pol(i) - aaa%zj(j))
                D(j) = -1.d0 / (aaa%pol(i) - aaa%zj(j))**2
            end do
            NN = dot_product(N(:), aaa%wj(:)*aaa%fj(:))
            DD = dot_product(D(:), aaa%wj(:))
            aaa%res(i) = NN / DD
            write(*,*) 'res=', aaa%res(i)
        end do
        write(*,*)
        deallocate(N, D)

        ! Compute zeros via generalized eigenvalue problem
        B(:,:) = zero
        E(:,:) = zero
        do i = 2, m+1
            B(i,i) = zone
            E(1,i) = aaa%wj(i-1)*aaa%fj(i-1)
            E(i,1) = zone
            E(i,i) = aaa%zj(i-1)
        end do
        ! do i = 1, m+1
        !     write(*,'(100f12.6)') dble(E(i,:))
        ! end do
        call mkl_zggev(ndim, E, B, alpha, beta)

        if (allocated(aaa%zer)) deallocate(aaa%zer)
        allocate(aaa%zer(aaa%nPol))
        j = 0
        do i = 1, ndim
            write(*,*) 'zer=', alpha(i) / beta(i)
            if (abs(beta(i)) > 1.d-16) then
                j = j+1
                aaa%zer(j) = alpha(i)/beta(i)
            end if
        end do
        write(*,*)
        deallocate(alpha, beta)

        deallocate(E, B)

    end subroutine


!--------------------------------------------------------------------------------
    subroutine cleanup(aaa, Z, F, cleanup_tol)
        implicit none
        type(aaa_approximant), intent(inout) :: aaa
        complex(8), intent(in) :: Z(:)
        complex(8), intent(in) :: F(:)
        real(8), intent(in) :: cleanup_tol
        ! local
        integer(4) :: i, j, jmin(1), m
        real(8) :: tol
        integer(4) :: nFD
        integer(4), allocatable :: idxFD(:), idxList(:)
        real(8), allocatable :: azp(:)

        integer(4) :: nj_new, m_new
        complex(8), allocatable :: zj_new(:), fj_new(:), wj_new(:)
        complex(8), allocatable :: RS(:,:), LS(:,:), C(:,:), A(:,:)
        real(8),    allocatable :: sval(:)
        complex(8), allocatable :: rsvec(:,:)
        
        ! Find negligible residues (Froissart doublets)
        allocate(idxFD(aaa%npol))
        nFD = 0
        do i = 1, aaa%npol
            if (abs(aaa%res(i)) < tol) then
                nFD = nFD+1
                idxFD(nFD) = i
            end if
        end do
        if (nFD == 0) then
            ! Nothing to do
            return
        else
            print*, 'Found ', nFD, ' Froissart doublets'
            do i = 1, nFD
                j = idxFD(i)
                print*, i, aaa%pol(j), aaa%res(j), aaa%zer(j)
            end do
        end if

        ! For each spurious pole find and remove closest support point
        allocate(azp(aaa%nj))
        allocate(idxList(nFD))
        do i = 1, nFD
            azp(:) = abs(aaa%zj(:) - aaa%pol(idxFD(i)))
            jmin = minloc(azp)
            idxList(i) = jmin(1)
        end do
        deallocate(azp)
        deallocate(idxFD)

        ! Remove corresponding support points
        nj_new = aaa%nj-nFD
        allocate(zj_new(nj_new), fj_new(nj_new))
        j = 0
        do i = 1, aaa%nj
            if (any(idxList /= i)) then
                j = j+1
                zj_new(j) = aaa%zj(i)
                fj_new(j) = aaa%fj(i)
                print*, dble(zj_new(j)), dble(fj_new(j))
            end if
        end do
        print*, 'j=', j
        deallocate(idxList)

        ! Remove support points z from sample set
        m = size(Z)
        m_new = 0
        allocate(idxList(m))
        do i = 1, m
            print*, 'Z=', dble(Z(i))
            if ( all( abs(Z(i)-zj_new) > 1.d-14 ) ) then
                m_new = m_new+1
                idxList(m_new) = i
            end if
        end do
        print*, "m_new=", m_new

        ! Build Loewner matrix
        allocate(LS(m_new,m_new))
        LS(:,:) = zero
        do i = 1, m_new
            LS(i,i) = F(idxList(i))
        end do
        allocate(RS(nj_new,nj_new))
        RS(:,:) = zero
        do i = 1, nj_new
            RS(i,i) = fj_new(i)
        end do
        allocate(C(m_new,nj_new))
        do i = 1, m_new
            C(i,:) = 1.d0 / (Z(idxList(i)) - zj_new(:))
        end do
        allocate(A(m_new,nj_new))
        A = matmul(LS, C) - matmul(C, RS)
        deallocate(LS, RS, C)
        
        ! Weight vector = min singular vector
        allocate(rsvec(nj_new,nj_new), sval(nj_new))
        call mkl_svd(m_new, nj_new, A, rsvec, sval)
        allocate(wj_new(nj_new))
        wj_new(:) = rsvec(:,nj_new)
        deallocate(A, rsvec, sval)

        ! Save results into an object
        call init_aaa_approximant(aaa, nj_new, zj_new, fj_new, wj_new)
        deallocate(zj_new, fj_new, wj_new)
        
    end subroutine

!--------------------------------------------------------------------------------
    subroutine mkl_zggev(ndim, A, B, alpha, beta)
        implicit none
        ! input/output
        integer(4), intent(in) :: ndim
        complex(8), intent(in) :: A(ndim, ndim)
        complex(8), intent(in) :: B(ndim, ndim)
        complex(8), intent(out) :: alpha(ndim)
        complex(8), intent(out) :: beta(ndim)
        ! local
        integer(4) :: ldvl, ldvr, lwork
        real(8), allocatable :: rwork(:)
        complex(8), allocatable :: work(:), vl(:,:), vr(:,:)
        integer(4) :: info 

        ldvl = max(1, ndim)
        ldvr = max(1, ndim)

        lwork = -1
        allocate(rwork(8*ndim))
        allocate(work(1))
        call zggev('N', 'N', ndim, A, ndim, B, ndim, alpha, beta, &
                   vl, ldvl, vr, ldvr, work, lwork, rwork, info)
        if (info /= 0) stop 'Error executing zggev!'
        lwork = work(1)
        deallocate(work)

        allocate(work(lwork))
        call zggev('N', 'N', ndim, A, ndim, B, ndim, alpha, beta, &
                   vl, ldvl, vr, ldvr, work, lwork, rwork, info)
        if (info /= 0) stop 'Error executing zggev!'
        deallocate(rwork, work)

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
        do i = 1, lmn
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
