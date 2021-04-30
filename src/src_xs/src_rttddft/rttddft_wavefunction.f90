! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created May 2019 (Ronaldo)
! Reference: https://arxiv.org/abs/2102.02630

module rttddft_Wavefunction
  use rttddft_GlobalVariables, only: ham_time, ham_past, overlap, &
    & evecfv_time, tstep, method
  use mod_kpoint, only: nkpt
  use mod_eigensystem, only: nmat
  use mod_eigenvalue_occupancy, only: nstfv
  use modinput, only: input
  use constants, only: zone, zzero, zi
  use precision, only: dp
  use modmpi
#ifdef USE_ASSERT
  use asserts, only: assert
#endif

  implicit none

  private

  public :: UpdateWavefunction

contains
  !> subroutine UpdateWavefunction: update wavefunction
  !> Here, we employ a propagator to evolve the Kohn-Sham wavefunctions
  !> Extrapolation scheme for the hamiltonian (predcorr .False.)
  !> \[ \hat{H}(t+f\Delta t) = (1+f)\hat{H}(t) - f\hat{H}(t-\Delta t). \]
  !> \(\hat{H}(t)\) is stored in "ham_time", whereas \(\hat{H}(t -\Delta t)\),
  !> in "ham_past"
  !> Extrapolation scheme for the hamiltonian (predcorr .True.)
  !> \[ \hat{H}(t+f\Delta t) = f\hat{H}(t + \Delta t) + (1-f)\hat{H}(t). \]
  !> \(\hat{H}(t)\) is stored in "ham_past", whereas \(\hat{H}(t +\Delta t)\),
  !> in "ham_time" (which comes from a previous iteration in the predictor corrector
  !> loop
  subroutine UpdateWavefunction( predcorr )

    implicit none
    !> @param[in]   predcorr  tells if we are in the loop of the predictor-Corrector scheme
    logical, intent(in)       :: predcorr
    !> Counter for loops with k-points
    integer                   :: ik
    !> Dimension of the Hamiltonian and Overlap matrices (given a k-point)
    integer                   :: nmatp
    !> indexes of the first and the last k-points
    integer                   :: first_kpt, last_kpt
    !> Factors that multiply the hamiltonian in the following propagator:
    !> Commutator-Free Magnus expansion of 4th order
    real(dp)                  :: f1, f2, a1, a2
    !> auxiliary variable to store the overlap and hamiltonian matrices
    complex(dp), allocatable  :: overl(:,:), ham(:,:), hamold(:,:)

#ifdef MPI
    first_kpt = firstk(rank)
    last_kpt = lastk(rank)
#else
    first_kpt = 1
    last_kpt = nkpt
#endif

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ik,nmatp,overl,ham,hamold), &
!$OMP& SHARED(first_kpt, last_kpt,a1,a2,f1,f2,nstfv,nkpt,method,predcorr), &
!$OMP& SHARED(input,nmat,tstep,ham_time,ham_past,overlap,evecfv_time)
!$OMP DO
#endif
    do ik = first_kpt, last_kpt
      ! Dimension of the Hamiltonian and Overlap matrices for the current k-point
      nmatp = nmat(1,ik)

      allocate(overl(nmatp,nmatp))
      allocate(ham(nmatp,nmatp))
      overl(1:nmatp,1:nmatp) = overlap(1:nmatp,1:nmatp,ik)

      select case(method)
        ! SE (simple exponential)
        ! CN (Crank-Nicolson)
        ! EMR (Exponential at midpoint rule)
        ! AETRS (approximate enforced time-reversal symmetry)
        ! CFM4 (Commutator-Free Magnus expansion of 4th order)
        ! EH (exponential using a basis of the hamiltonian-eigenvectors)
        ! EHM (same as before, but uses the hamiltonian at midpoint)
        ! RK4 (Runge-Kutta of 4th order)
        case ('SE')
          ham(1:nmatp,1:nmatp) = ham_time(1:nmatp,1:nmatp,ik)
          call exponential( ik, nmatp, ham, overl )
        case ('EMR')
          if ( .not. predcorr ) then
            ham(1:nmatp,1:nmatp) = 1.5_dp*ham_time(1:nmatp,1:nmatp,ik) -0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
          else
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_time(1:nmatp,1:nmatp,ik) +0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
          end if
          call exponential(ik,nmatp,ham,overl)
        case ('AETRS')
          if ( .not. predcorr ) then
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_time(1:nmatp,1:nmatp,ik)
            call exponential(ik,nmatp,ham,overl)
            ham(1:nmatp,1:nmatp) = ham_time(1:nmatp,1:nmatp,ik)-0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
            call exponential(ik,nmatp,ham,overl)
          else
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
            call exponential(ik,nmatp,ham,overl)
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_time(1:nmatp,1:nmatp,ik)
            call exponential(ik,nmatp,ham,overl)
          end if
        case ('CFM4')
          f1 =  0.21132486540518713_dp ! 1/2 - sqrt(3)/6
          f2 =  0.78867513459481290_dp ! 1/2 + sqrt(3)/6
          a1 = -0.03867513459481287_dp ! 1/4 - sqrt(3)/6
          a2 =  0.53867513459481290_dp ! 1/4 + sqrt(3)/6
          if ( .not. predcorr ) then
            ham(1:nmatp,1:nmatp) = a1*((1+f2)*ham_time(1:nmatp,1:nmatp,ik)-f2*ham_past(1:nmatp,1:nmatp,ik)) + &
                                   a2*((1+f1)*ham_time(1:nmatp,1:nmatp,ik)-f1*ham_past(1:nmatp,1:nmatp,ik))
            call exponential(ik,nmatp,ham,overl)
            ham(1:nmatp,1:nmatp) = a1*((1+f1)*ham_time(1:nmatp,1:nmatp,ik)-f1*ham_past(1:nmatp,1:nmatp,ik)) + &
                                   a2*((1+f2)*ham_time(1:nmatp,1:nmatp,ik)-f2*ham_past(1:nmatp,1:nmatp,ik))
            call exponential(ik,nmatp,ham,overl)
          else
            ham(1:nmatp,1:nmatp) = a1*((1-f2)*ham_past(1:nmatp,1:nmatp,ik)+f2*ham_time(1:nmatp,1:nmatp,ik)) + &
                                   a2*((1-f1)*ham_past(1:nmatp,1:nmatp,ik)+f1*ham_time(1:nmatp,1:nmatp,ik))
            call exponential(ik,nmatp,ham,overl)
            ham(1:nmatp,1:nmatp) = a1*((1-f1)*ham_past(1:nmatp,1:nmatp,ik)+f1*ham_time(1:nmatp,1:nmatp,ik)) + &
                                   a2*((1-f2)*ham_past(1:nmatp,1:nmatp,ik)+f2*ham_time(1:nmatp,1:nmatp,ik))
            call exponential(ik,nmatp,ham,overl)
          end if
        case ('EH')
          ham(1:nmatp,1:nmatp) = ham_time(1:nmatp,1:nmatp,ik)
          call exphouston(ik,nmatp,ham,overl)
        case ('EHM')
          if ( .not. predcorr ) then
            ham(1:nmatp,1:nmatp) = 1.5_dp*ham_time(1:nmatp,1:nmatp,ik) -0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
          else
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_time(1:nmatp,1:nmatp,ik) +0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
          end if
          call exphouston(ik,nmatp,ham,overl)
        case ('RK4')
          ham(1:nmatp,1:nmatp) = ham_time(1:nmatp,1:nmatp,ik)
          allocate(hamold(1:nmatp,1:nmatp))
          hamold(1:nmatp,1:nmatp) = ham_past(1:nmatp,1:nmatp,ik)
          call rk4(ik,nmatp,ham,hamold,overl,predcorr)
          deallocate(hamold)
      end select

      ! Normalize WFs, if this is the case
      ! The propagator operator should be unitary, so this step would be unnecessary
      ! However, numerically this is almost never possible
      ! This normalization may help to avoid numerical issues
      if ( input%xs%realTimeTDDFT%normalizeWF ) call normalize( ik, nmatp, overl )

      deallocate(overl)
      deallocate(ham)
    end do
#ifdef USEOMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

  end subroutine UpdateWavefunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> subroutine exponential: obtain the exponential of H and applies to the
  !> coefficients that characterize KS wavefunctions in terms of the basis
  !> Because the basis is not orthonormal, we have
  !> \[
  !>    C_{j\mathbf{k}}(t+\Delta t) = \mathrm{exp}[-\mathrm{i}\Delta t
  !>    S_{\mathbf{k}}^{-1}H_{\mathbf{k}}(t)] \; C_{j\mathbf{k}}(t)
	!> \]
  !> \( \C_{j\mathbf{k}}(t) \) is stored in "evecfv_time" (global variable)
  !> \( H_{\mathbf{k}}(t) \) is stored in "ham"
  !> \( S_{\mathbf{k}} \) is stored in "overl"
  !> The exponential here is approximated by a Taylor expansion
  !> up to the order defined by input%xs%rt_tddft%order_taylor
  subroutine exponential( ik, nmatp, ham, overl )
    implicit none

    !> index of the k-point considered
    integer, intent(in)       :: ik
    !> Dimension of the Hamiltonian and Overlap matrices (given a k-point)
    integer, intent(in)       :: nmatp
    !> Hamiltonian matrix
    complex(dp),intent(in)    :: ham(nmatp,nmatp)
    !> Overlap matrix
    complex(dp),intent(inout) :: overl(nmatp,nmatp)

    ! Counter for the main loop
    integer                   :: it
    ! Auxiliary variables for the Lapack routines
    integer                   :: info, lwork
    integer, allocatable      :: ipiv(:)
    ! time step multiplied by the imaginary number
    complex(dp)               :: zit
    ! Auxiliary variables
    complex(dp), allocatable  :: scratch(:, :), store(:, :)
    complex(dp), allocatable  :: work(:)


    allocate( scratch(nmatp, nstfv) )
    allocate( store(nmatp, nstfv) )
    allocate( ipiv(nmatp) )


    zit = (-tstep) * zi
    scratch(:, :) = zzero
    store(1:nmatp, 1:nstfv) = evecfv_time(1:nmatp, 1:nstfv, ik)

    ! Taylor expansion
    do it = 1, input%xs%realTimeTDDFT%TaylorOrder
      ! Matrix multiplication C := alpha*AB+beta*C
      call ZHEMM( 'L', 'U', nmatp, nstfv, zone, ham, nmatp, store, &
        & nmatp, zzero, scratch, nmatp )
      if ( it == 1) then
        ! Gets the best size of lwork
        allocate(work(2))
        lwork = -1
        call ZHESV( 'U', nmatp, nstfv, overl, nmatp, ipiv, scratch, nmatp, &
          & work, lwork, info )
        lwork = int( work(1) )
        deallocate( work )
        allocate( work(lwork) )
      end if
      overl(1:nmatp, 1:nmatp) = overlap(1:nmatp, 1:nmatp,ik)
      ! Solves Ax = B
      call ZHESV( 'U', nmatp, nstfv, overl, nmatp, ipiv, scratch, nmatp, &
        & work, lwork, info )
      store(1:nmatp,1:nstfv) = (zit/it)*scratch(1:nmatp,1:nstfv)
      evecfv_time(1:nmatp,1:nstfv,ik) = evecfv_time(1:nmatp,1:nstfv,ik) + &
        & store(1:nmatp,1:nstfv)
    end do

    deallocate(scratch)
    deallocate(store)
    deallocate(ipiv)
    deallocate(work)

  end subroutine exponential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> subroutine exphouston: exponential using the so-called Houston expansion
  !> Here, we evaluate the exponential operator exactly rather than
  !> Taylor-expanding it. This can be done by taking into account an adiabatic
  !> basis formed by the eigenvectors of \( \hat{H}(t) \),
  !> which means that for each \( t \), we solve
  !>  \[
  !>       \hat{H}(t)| \phi_{j\mathbf{k}}(t)\rangle =
	!>       \varepsilon_{j\mathbf{k}}(t)| \phi_{j\mathbf{k}}(t)\rangle
  !>  \]
  !> and then expand
  !>  \[
  !>      | \psi_{j\mathbf{k}}(t)\rangle = \sum_i \alpha_{ij\mathbf{k}}(t)
  !>      | \phi_{i\mathbf{k}}(t)\rangle,
  !>  \]
  !> where
  !> \( \alpha_{ij\mathbf{k}}(t) = \langle \phi_{i}^{\mathbf{k}}(t)| \psi_{j}^{\mathbf{k}}(t)\rangle \)
  !> Since
  !>  \[
  !>	   \hat{U}(t+\Delta t,t) |\phi_{m\mathbf{k}}(t)\rangle =
	!>     \mathrm{e}^{-\mathrm{i}\varepsilon_{m\mathbf{k}}(t) \Delta t}
	!>     |\phi_{m\mathbf{k}}(t)\rangle,
	!>  \]
  !> when considering \( \hat{U}(t+\Delta t,t) \) in the form of the
  !> simple exponential propagator, the action of the propagator is
  !>  \[
  !>      \hat{U}(t+\Delta t,t)	| \psi_{j\mathbf{k}}(t)\rangle =
	!>       \sum_i \alpha_{ij\mathbf{k}}(t)
  !>      \mathrm{e}^{-\mathrm{i}\varepsilon_{m\mathbf{k}}(t) \Delta t}
	!>       |\phi_{m\mathbf{k}}(t)\rangle.
  !> \]
  subroutine exphouston( ik, nmatp, ham, overl )
    implicit none

    !> index of the k-point considered
    integer, intent(in)       :: ik
    !> Dimension of the Hamiltonian and Overlap matrices (given a k-point)
    integer, intent(in)       :: nmatp
    !> Hamiltonian matrix
    complex(dp),intent(in)    :: ham(nmatp,nmatp)
    !> Overlap matrix
    complex(dp),intent(in)    :: overl(nmatp,nmatp)

    integer                   :: i,m,lwork,info
    integer                   :: ifail(nmatp),iwork(5*nmatp)
    real (dp)                 :: vl,vu
    real (dp)                 :: w(nmatp)
    complex(dp)               :: zit
    complex(dp)               :: rwork(7*nmatp)
    complex(dp), allocatable  :: store(:,:),evecham(:,:)
    complex(dp), allocatable  :: alpha(:,:),scratch(:,:)
    complex(dp), allocatable  :: overlcopy(:,:), hamcopy(:,:)
    complex(dp), allocatable  :: work(:)

    ! Initialization
    allocate(store(nmatp,nstfv))
    allocate(evecham(nmatp,nstfv))
    allocate(scratch(nmatp,nstfv))
    allocate(alpha(nstfv,nstfv))
    allocate(overlcopy(nmatp,nmatp))
    allocate(hamcopy(nmatp,nmatp))
    store(1:nmatp,1:nstfv) = evecfv_time(1:nmatp,1:nstfv,ik)
    overlcopy(1:nmatp,1:nmatp) = overl(1:nmatp,1:nmatp)
    hamcopy(1:nmatp,1:nmatp) = ham(1:nmatp,1:nmatp)

    vl = 0._dp
    vu = 0._dp
    zit = (-tstep) * zi

    ! Obtain the optimum lwork
    lwork = -1
    allocate(work(2))
    call ZHEGVX(1,'V','I','U',nmatp,hamcopy,nmatp,overlcopy,nmatp,vl,vu,1,nstfv, &
        input%groundstate%solver%evaltol,m,w,store,nmatp,work,lwork,rwork,&
        iwork,ifail,info)
    lwork = int( work(1) )
    deallocate(work)
    allocate(work(lwork))

    ! Obtain the eigenvalues/eigenvectors of the hamiltonian at current time
    call ZHEGVX( 1, &    ! Specifies the problem type to be solved: = 1:  A*x = (lambda)*B*x
                'V', &   ! 'V':  Compute eigenvalues and eigenvectors.
                'I', &   ! 'I': the IL-th through IU-th eigenvalues will be found.
                'U', &   ! 'U': Upper triangles of A and B are stored;
                nmatp, & ! N: The order of the matrices A and B.
                hamcopy,&! A(in,out): On entry, the Hermitian matrix A. On exit,
                         ! the upper triangle (if UPLO='U') of A, including the diagonal, is destroyed.
                nmatp, & ! LDA: The leading dimension of the array A.
                overlcopy, & ! B(in,out): On entry, the Hermitian matrix B.
                         ! On exit, if INFO <= N, the part of B containing the matrix is
                         ! overwritten by the triangular factor U or L from the Cholesky factorization B = U**H*U or B = L*L**H.
                nmatp,&  ! LDB: The leading dimension of the array B.
                vl,vu,&  ! Only matters if it had not been chosen 'I' before
                1, &     ! the index of the smallest eigenvalue to be returned.
                nstfv, & ! the index of the largest eigenvalue to be returned
                input%groundstate%solver%evaltol, & ! The absolute error tolerance for the eigenvalues
                m, &     ! M(out): The total number of eigenvalues found.
                w, &     ! W(out): is DOUBLE PRECISION array, dimension (N).
                         ! The first M elements contain the selected eigenvalues in ascending order.
                evecham,&! Z(out): dimension (LDZ, max(1,M)). The orthonormal eigenvectors
                nmatp,&  ! LDZ: The leading dimension of the array Z.
                work,&   ! WORK(out) is COMPLEX*16 array, dimension (MAX(1,LWORK))
                lwork,&  ! LWORK is INTEGER
                rwork, & ! RWORK(out) is DOUBLE PRECISION array, dimension (7*N)
                iwork, & ! IWORK(out) is INTEGER array, dimension (5*N)
                ifail, & ! IFAIL(out) is INTEGER array, dimension (N)
                info )   ! INFO(out) is INTEGER = 0:  successful exit

    ! Check if there were problems with the diagonalization
    if ( info /= 0 ) then
#ifdef USE_ASSERT
      call assert( .false., &
         & 'Problems with ZHEGVX - exphouston' )
#else
      write(*,*) 'Problems with ZHEGVX - exphouston'
      stop
#endif
    end if

    ! Project the current WFs onto these eigenvectors
    ! Matrix multiplication C := alpha*AB+beta*C
    call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
               'N', &           ! TRANSB = 'N'  op( B ) = B.
                nmatp, &          ! M ... rows of op( A ) = rows of C
                nstfv, &           ! N ... cols of op( B ) = cols of C
                nmatp, &          ! K ... cols of op( A ) = rows of op( B )
                zone, &          ! alpha
                overl, &           ! A
                nmatp,&           ! LDA ... leading dimension of A
                store, &           ! B
                nmatp, &          ! LDB ... leading dimension of B
                zzero, &          ! beta
                scratch, &  ! C
                nmatp &      ! LDC ... leading dimension of C
                )
    call zgemm( 'C', &           ! TRANSA = 'C' or 'c',  op( A ) = A**H.
                'N', &           ! TRANSB = 'N'  op( B ) = B.
                nstfv, &          ! M ... rows of op( A ) = rows of C
                nstfv, &           ! N ... cols of op( B ) = cols of C
                nmatp, &          ! K ... cols of op( A ) = rows of op( B )
                zone, &          ! alpha
                evecham, &         ! A
                nmatp,&           ! LDA ... leading dimension of A
                scratch, &           ! B
                nmatp, &          ! LDB ... leading dimension of B
                zzero, &          ! beta
                alpha, &  ! C
                nstfv &      ! LDC ... leading dimension of C
                )
    ! Now, evolve
    do i=1,nstfv
      evecham(1:nmatp,i) = zexp( zit*w(i) )*evecham(1:nmatp,i)
    end do

    call zgemm( 'N', &           ! TRANSA = 'N'  op( A ) = A.
                'N', &           ! TRANSB = 'N'  op( B ) = B.
                nmatp, &          ! M ... rows of op( A ) = rows of C
                nstfv, &           ! N ... cols of op( B ) = cols of C
                nstfv, &          ! K ... cols of op( A ) = rows of op( B )
                zone, &          ! alpha
                evecham, &         ! A
                nmatp,&           ! LDA ... leading dimension of A
                alpha, &           ! B
                nstfv, &          ! LDB ... leading dimension of B
                zzero, &          ! beta
                store, &  ! C
                nmatp &      ! LDC ... leading dimension of C
                )
    evecfv_time(1:nmatp,1:nstfv,ik) = store(1:nmatp,1:nstfv)
    deallocate(store)
    deallocate(evecham)
    deallocate(alpha)
    deallocate(work)
  end subroutine exphouston
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> subroutine rk4: Runge-Kutta of 4th order
  !> The basic formulas employed here are
  !> \[
  !>	   | \psi_{j\mathbf{k}}(t+\Delta t)\rangle = | \psi_{j\mathbf{k}}(t) \rangle
	!>    	- \frac{\mathrm{i}\Delta t}{6}(k_1+2k_2+2k_3+k_4)
  !>  \]
  !>	where
  !> \[
  !>	   k_1 = S_k^{-1}H_k(t)| \psi_{j\mathbf{k}}(t) \rangle,
  !>  \]
  !> \[
	!>     k_2 = S_k^{-1}H_k\left(t + \frac{\Delta t}{2}\right)
  !>    \left[| \psi_{j\mathbf{k}}(t) \rangle + k_1 \frac{\Delta t}{2}\right],
  !>  \]
  !> \[
	!>	   k_3 = S_k^{-1}H_k\left(t + \frac{\Delta t}{2}\right)
  !>    \left[| \psi_{j\mathbf{k}}(t) \rangle + k_2 \frac{\Delta t}{2}\right],
  !>  \]
  !> \[
	!>     k_4 = S_k^{-1}H_k(t+\Delta t)\left[| \psi_{j\mathbf{k}}(t) \rangle
  !>      + k_3\Delta t \right].
  !>  \]
  subroutine rk4( ik, nmatp, ham, hamold, overl, predcorr )
    implicit none
    !> index of the k-point considered
    integer, intent(in)     :: ik
    !> Dimension of the Hamiltonian and Overlap matrices (given a k-point)
    integer, intent(in)     :: nmatp
    !> Hamiltonian matrix at time \( t \)
    complex(dp),intent(in)  :: ham(nmatp,nmatp)
    !> Hamiltonian matrix at time \( t - \Delta t \)
    complex(dp),intent(in)  :: hamold(nmatp,nmatp)
    !> Overlap matrix
    complex(dp),intent(in)  :: overl(nmatp,nmatp)
    !> Predictor Corrector scheme is used (True) or not (False)
    !> If True, then "ham" and "hamold" mean the hamiltonian matrices at times
    !> \( t + \Delta t \) and ( t \), respectively
    logical,intent(in)      :: predcorr

    integer                 :: i, info, lwork, ipiv(nmatp)

    complex(dp)             :: zit
    complex(dp)             :: k(nmatp,nstfv,4)
    complex(dp)             :: psi(nmatp,nstfv)
    complex(dp)             :: hamcopy(nmatp,nmatp)
    complex(dp)             :: overlcopy(nmatp,nmatp)
    complex(dp),allocatable :: work(:)

    zit = -tstep*zi

    psi(1:nmatp,1:nstfv) = evecfv_time(1:nmatp,1:nstfv,ik)
    hamcopy(1:nmatp,1:nmatp) = ham(1:nmatp,1:nmatp)
    do i = 1, 4
      ! ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      ! C := alpha*A*B + beta*C, or alpha*B*A + beta*C,
      call zhemm('L', 'U', nmatp, nstfv, zit,&
        & hamcopy(:,:), nmatp,&
        & psi(:,:), nmatp, zzero, k(:,:,i), nmatp)

      ! ZHESV(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,LWORK,INFO)
      ! Solves A*X = B, A is an N-by-N Hermitian matrix
      overlcopy(1:nmatp,1:nmatp) = overl(1:nmatp,1:nmatp)
      lwork = -1
      allocate(work(2))
      call zhesv('U', nmatp, nstfv, overlcopy, nmatp, ipiv, k(:,:,i), &
        & nmatp, work, lwork, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      call zhesv('U', nmatp, nstfv, overlcopy, nmatp, ipiv, k(:,:,i), &
        & nmatp, work, lwork, info)
      deallocate(work)

      select case(i)
        case(1)
          psi(1:nmatp,1:nstfv) = evecfv_time(1:nmatp,1:nstfv,ik) + &
            & 0.5_dp*k(1:nmatp,1:nstfv,1)
          if (predcorr) then
            hamcopy(1:nmatp,1:nmatp) = 0.5_dp*ham(1:nmatp,1:nmatp) + &
              & 0.5_dp*hamold(1:nmatp,1:nmatp)
          else
            hamcopy(1:nmatp,1:nmatp) = 1.5_dp*ham(1:nmatp,1:nmatp) - &
              & 0.5_dp*hamold(1:nmatp,1:nmatp)
          end if
        case(2)
          psi(1:nmatp,1:nstfv) = evecfv_time(1:nmatp,1:nstfv,ik) + &
            & 0.5_dp*k(1:nmatp,1:nstfv,2)
        case(3)
          psi(1:nmatp,1:nstfv) = evecfv_time(1:nmatp,1:nstfv,ik) + &
            & k(1:nmatp,1:nstfv,3)
          if ( predcorr ) then
            hamcopy(1:nmatp,1:nmatp) = ham(1:nmatp,1:nmatp)
          else
            hamcopy(1:nmatp,1:nmatp) = 2._dp*ham(1:nmatp,1:nmatp) - &
              & hamold(1:nmatp,1:nmatp)
          end if
        case(4)
          exit
      end select
    end do


    evecfv_time(1:nmatp,1:nstfv,ik) = evecfv_time(1:nmatp,1:nstfv,ik) + &
      & (1._dp/6._dp)*( k(1:nmatp,1:nstfv,1) + 2._dp*k(1:nmatp,1:nstfv,2) + &
      & 2._dp*k(1:nmatp,1:nstfv,3) + k(1:nmatp,1:nstfv,4) )

  end subroutine rk4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> subroutine normalize: imposes that the norm of the KS wavefunctions are 1
!> The norm is calculated as
!> \[
!>  C_{j\mathbf{k}}^\dagger S_{\mathbf{k}} C_{j\mathbf{k}}
!>  \]
!> where \( C_{j\mathbf{k}} \) is the array of coefficients of the wavefunction
!> in terms of the basis, and (\ S_{\mathbf{k}} \) is the overlap matrix.
  subroutine normalize( ik, nmatp, overl )
    implicit none

    !> index of the k-point considered
    integer, intent(in)       :: ik
    !> Dimension of the Hamiltonian and Overlap matrices (given a k-point)
    integer, intent(in)       :: nmatp
    !> Overlap matrix
    complex(dp),intent(inout) :: overl(nmatp,nmatp)

    integer                   :: ist
    complex(dp)               :: scratch(nmatp,nstfv)
    complex(dp)               :: store(nmatp,nstfv)
    complex(dp)               :: zdotc
    complex(dp)               :: norm

    store(1:nmatp,1:nstfv) = evecfv_time(1:nmatp,1:nstfv,ik)
    overl(1:nmatp,1:nmatp) = overlap(1:nmatp,1:nmatp,ik)

    ! TODO: replace zgemm with zhemm
    call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
               'N', &           ! TRANSB = 'N'  op( B ) = B.
               nmatp, &          ! M ... rows of op( A ) = rows of C
               nstfv, &           ! N ... cols of op( B ) = cols of C
               nmatp, &          ! K ... cols of op( A ) = rows of op( B )
               zone, &          ! alpha
               overl, &           ! A
               nmatp,&           ! LDA ... leading dimension of A
               store(:,:), &           ! B
               nmatp, &          ! LDB ... leading dimension of B
               zzero, &          ! beta
               scratch(:,:), &  ! C
               nmatp &      ! LDC ... leading dimension of C
               )
    do ist = 1,nstfv
      norm = dsqrt(dble(zdotc(nmatp,store(:,ist),1,scratch(:,ist),1)))
      evecfv_time(:,ist,ik) = evecfv_time(:,ist,ik)/norm
    end do
  end subroutine normalize
end module rttddft_Wavefunction
