!> Unit tests for Hermitian eigensolvers
module herm_eigensolver_test
  use precision, only: dp, i32
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close, identity_complex_dp, fill_random
  use modscl
  use m_dzmatmult, only: dzmatmult
  use herm_eigensolver, only: lapack_eigensolver
  use exciting_mpi, only: xmpi_bcast
#ifdef SCAL
  use herm_eigensolver, only: scalapack_eigensolver_pzheevx, scalapack_eigensolver_pzheevd
  use modmpi, only: MPI_UNDEFINED, MPI_Comm_split
#endif
#ifdef _ELPA_
  use herm_eigensolver, only: elpa_eigensolver
#endif
  
  implicit none
  
  private
  public :: herm_eigensolver_test_driver
  
contains

  !> Run tests for Hermitian eigensolvers
  subroutine herm_eigensolver_test_driver(mpiglobal, kill_on_failure)
    !> mpiglobal:MPI environment
    type(mpiinfo), intent(inout) :: mpiglobal
    !> kill_on_failure:Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure
    
    !> test_report:Test report object
    type(unit_test_type) :: test_report
    
    call test_report%init(mpiglobal)
    
    ! Test LAPACK solver (always available)
    call test_lapack_solver(test_report, mpiglobal)
    
#ifdef SCAL
    ! Test ScaLAPACK pzheevx solver (if compiled with ScaLAPACK)
    call test_scalapack_solver(test_report, mpiglobal, 'pzheevx')
    ! Test ScaLAPACK pzheevd solver (if compiled with ScaLAPACK)
    call test_scalapack_solver(test_report, mpiglobal, 'pzheevd')
#endif
    
#ifdef _ELPA_
    ! Test ELPA 2stage solver (if compiled with ELPA)
    call test_elpa_solver(test_report, mpiglobal, '2')
    ! Test ELPA 1stage solver (if compiled with ELPA)
    call test_elpa_solver(test_report, mpiglobal, '1')
#endif
    
    if (present(kill_on_failure)) then
      call test_report%report('herm_eigensolver', kill_on_failure)
    else
      call test_report%report('herm_eigensolver')
    end if
    
    call test_report%finalise()
  end subroutine herm_eigensolver_test_driver

  !> Test LAPACK solver with a random Hermitian matrix
  subroutine test_lapack_solver(test_report, mpiglobal)
    !> test_report:Test report object
    type(unit_test_type), intent(inout) :: test_report
    !> mpiglobal:MPI environment
    type(mpiinfo), intent(inout) :: mpiglobal
    
    !> N:Matrix dimension
    integer, parameter :: N = 5
    !> ham:Hermitian matrix to diagonalize
    type(dzmat) :: ham
    !> ham_orig:Original Hamiltonian (for verification)
    type(dzmat) :: ham_orig
    !> evec:Eigenvectors
    type(dzmat) :: evec
    !> binfo:BLACS context information
    type(blacsinfo) :: binfo
    !> eval:Eigenvalues
    real(dp) :: eval(N)
    !> H_matrix:Test Hermitian matrix (global array)
    complex(dp) :: H_matrix(N,N)
    !> eigenvalue_eq_ok:Flag for eigenvalue equation check
    logical :: eigenvalue_eq_ok
    !> evecs_orthonormal:Flag for eigenvector orthonormality check
    logical :: evecs_orthonormal
    !> i:Loop index
    integer :: i
    !> tol:Tolerance for numerical comparison
    real(dp), parameter :: tol = 1e-9_dp
    
    ! Create a positive definite Hermitian matrix: H = A^H * A + alpha * I
    ! This guarantees full rank (linear independence of all columns)
    call fill_random(H_matrix)
    H_matrix = matmul(conjg(transpose(H_matrix)), H_matrix)
    do i = 1, N
      H_matrix(i,i) = H_matrix(i,i) + cmplx(1.0_dp, 0.0_dp, dp)
    end do
    
    ! Initialize BLACS context (non-distributed: 1×1 grid for LAPACK)
    call setupblacs(mpiglobal, '0d', binfo, np=1)
    call new_dzmat(ham, N, N, binfo)
    call new_dzmat(ham_orig, N, N, binfo)
    call new_dzmat(evec, N, N, binfo)
    
    ! Copy global matrix to distributed matrices
    call dzmat_copy_global2local(H_matrix, ham, binfo)
    call dzmat_copy_global2local(H_matrix, ham_orig, binfo)
    
    ! Run LAPACK eigensolver directly (destroys ham)
    call lapack_eigensolver(ham, eval, binfo, evec)
    
    ! Test 1: Check eigenvector orthonormality (X^H * X = I)
    call check_orthonormality(evec, binfo, mpiglobal, evecs_orthonormal, tol)
    call test_report%assert(evecs_orthonormal, &
         'Test if eigenvectors are orthonormal for random Hermitian matrix (X^H * X = I)')
    
    call check_eigenvalue_equation(ham_orig, evec, eval, binfo, mpiglobal, eigenvalue_eq_ok, tol)
    call test_report%assert(eigenvalue_eq_ok, &
         'Test eigenvalue equation for random Hermitian matrix (H * X_i = lambda_i * X_i)')
    
    call del_dzmat(ham)
    call del_dzmat(ham_orig)
    call del_dzmat(evec)
    call exitblacs(binfo)
  end subroutine test_lapack_solver

#ifdef SCAL
  !> Test ScaLAPACK solver with a random Hermitian matrix  
  subroutine test_scalapack_solver(test_report, mpiglobal, scala_solver)
    !> test_report:Test report object
    type(unit_test_type), intent(inout) :: test_report
    !> mpiglobal:MPI environment
    type(mpiinfo), intent(inout) :: mpiglobal
    !> which ScaLAPACK solver to use, must be 'pzheevx' or 'pzheevd'
    character(7), intent(in) :: scala_solver

    !> N:Matrix dimension
    integer, parameter :: N = 200
    !> ham:Hermitian matrix to diagonalize
    type(dzmat) :: ham
    !> ham_orig:Original Hamiltonian (for verification)
    type(dzmat) :: ham_orig
    !> evec:Eigenvectors
    type(dzmat) :: evec
    !> binfo:BLACS context information
    type(blacsinfo) :: binfo
    !> eval:Eigenvalues
    real(dp) :: eval(N)
    !> H_matrix:Test Hermitian matrix (global array)
    complex(dp) :: H_matrix(N,N)
    !> eigenvalue_eq_ok:Flag for eigenvalue equation check
    logical :: eigenvalue_eq_ok
    !> evecs_orthonormal:Flag for eigenvector orthonormality check
    logical :: evecs_orthonormal
    !> i:Loop index
    integer :: i
    !> tol:Tolerance for numerical comparison
    real(dp), parameter :: tol = 1e-9_dp
    
    ! Create a positive definite Hermitian matrix: H = A^H * A + alpha * I
    ! This guarantees full rank (linear independence of all columns)
    call fill_random(H_matrix)
    H_matrix = matmul(conjg(transpose(H_matrix)), H_matrix)
    do i = 1, N
      H_matrix(i,i) = H_matrix(i,i) + cmplx(1.0_dp, 0.0_dp, dp)
    end do
    
    ! Initialize BLACS context (parallel: use all available processes for ScaLAPACK)
    call setupblacs(mpiglobal, '2d', binfo, np=mpiglobal%procs)

    ! guard for idle ranks
    if(binfo%isactive) then
      call new_dzmat(ham, N, N, binfo)
      call new_dzmat(ham_orig, N, N, binfo)
      call new_dzmat(evec, N, N, binfo)
      call dzmat_copy_global2local(H_matrix, ham, binfo)
      call dzmat_copy_global2local(H_matrix, ham_orig, binfo)

      ! Run ScaLAPACK eigensolver directly (destroys ham)
      if (scala_solver == 'pzheevx') then
        call scalapack_eigensolver_pzheevx(ham, eval, binfo, evec)
      else if (scala_solver == 'pzheevd') then
        call scalapack_eigensolver_pzheevd(ham, eval, binfo, evec)
      else
        call terminate('Error(test_scalapack_solver): Solver name not recognized, use "pzheevx" or "pzheevd".')
      end if

      ! Test 1: Check eigenvector orthonormality
      call check_orthonormality(evec, binfo, mpiglobal, evecs_orthonormal, tol)
      call test_report%assert(evecs_orthonormal, &
           'Test ScaLAPACK: eigenvectors orthonormal (X^H * X = I)')

      call check_eigenvalue_equation(ham_orig, evec, eval, binfo, mpiglobal, eigenvalue_eq_ok, tol)
      call test_report%assert(eigenvalue_eq_ok, &
           'Test ScaLAPACK: eigenvalue equation (H * X = Lambda * X)')

      call del_dzmat(ham)
      call del_dzmat(ham_orig)
      call del_dzmat(evec)
      call blacsbarrier(binfo)
     end if
    call exitblacs(binfo)
  end subroutine test_scalapack_solver
#endif

#ifdef _ELPA_
  !> Test ELPA solver with a random Hermitian matrix
  subroutine test_elpa_solver(test_report, mpiglobal, elpa_solver)
    !> test_report:Test report object
    type(unit_test_type), intent(inout) :: test_report
    !> mpiglobal:MPI environment
    type(mpiinfo), intent(inout) :: mpiglobal
    !> which ELPA solver to use, must be '1' or '2'
    character(1), intent(in) :: elpa_solver

    !> N:Matrix dimension
    integer, parameter :: N = 200
    !> ham:Hermitian matrix to diagonalize
    type(dzmat) :: ham
    !> ham_orig:Original Hamiltonian (for verification)
    type(dzmat) :: ham_orig
    !> evec:Eigenvectors
    type(dzmat) :: evec
    !> binfo:BLACS context information
    type(blacsinfo) :: binfo
    !> eval:Eigenvalues
    real(dp) :: eval(N)
    !> H_matrix:Test Hermitian matrix (global array)
    complex(dp) :: H_matrix(N,N)
    !> eigenvalue_eq_ok:Flag for eigenvalue equation check
    logical :: eigenvalue_eq_ok
    !> evecs_orthonormal:Flag for eigenvector orthonormality check
    logical :: evecs_orthonormal
    !> i:Loop index
    integer :: i
    !> tol:Tolerance for numerical comparison
    real(dp), parameter :: tol = 1e-9_dp
    !> total number of processes
    integer :: nprocs
    !> number of processes in columns/rows
    integer :: npcols,nprows
    !> number of used processes in the 2d grid
    integer :: nprocs2d
    !> specifier if rank is active or not
    integer :: color
    !> MPI communicator of the active BLACS grid
    integer :: comm_active
    !> error status
    integer :: ierror
    !> mpi instance only for active ranks, used for ELPA
    type(mpiinfo) :: mpicom_active

    ! Create a positive definite Hermitian matrix: H = A^H * A + alpha * I
    ! This guarantees full rank (linear independence of all columns)
    call fill_random(H_matrix)
    H_matrix = matmul(conjg(transpose(H_matrix)), H_matrix)
    do i = 1, N
      H_matrix(i,i) = H_matrix(i,i) + cmplx(1.0_dp, 0.0_dp, dp)
    end do
    
    ! Initialize BLACS context (parallel: use all available processes for ELPA)
    nprocs = mpiglobal%procs
    npcols = int(sqrt(dble(nprocs)))
    nprows = nprocs / npcols
    nprocs2d = npcols*nprows
    if(mpiglobal%rank < nprocs2d) then
      color = 1
    else
      color = MPI_UNDEFINED
    end if
    call MPI_Comm_split(mpiglobal%comm, color, 0, comm_active, ierror)
    if (color /= MPI_UNDEFINED) then
      call mpicom_active%init(comm_active)
      call setupblacs(mpicom_active, 'grid', binfo, np=nprocs2d)
    else
      call setupblacs(mpiglobal, 'xxx', binfo, np=nprocs2d)
    end if
    if(binfo%isactive) then
      call new_dzmat(ham, N, N, binfo)
      call new_dzmat(ham_orig, N, N, binfo)
      call new_dzmat(evec, N, N, binfo)
      call dzmat_copy_global2local(H_matrix, ham, binfo)
      call dzmat_copy_global2local(H_matrix, ham_orig, binfo)

      ! Run ELPA eigensolver directly (destroys ham)
      call elpa_eigensolver(ham, eval, binfo, elpa_solver, evec)

      ! Test 1: Check eigenvector orthonormality
      call check_orthonormality(evec, binfo, mpiglobal, evecs_orthonormal, tol)
      call test_report%assert(evecs_orthonormal, &
           'Test ELPA: eigenvectors orthonormal (X^H * X = I)')

      call check_eigenvalue_equation(ham_orig, evec, eval, binfo, mpiglobal, eigenvalue_eq_ok, tol)
      call test_report%assert(eigenvalue_eq_ok, &
           'Test ELPA: eigenvalue equation (H * X = Lambda * X)')

      call del_dzmat(ham)
      call del_dzmat(ham_orig)
      call del_dzmat(evec)
      call blacsbarrier(binfo)
    end if
    call exitblacs(binfo)
    call barrier(callername='test_elpa_solver')
  end subroutine test_elpa_solver
#endif

  !> Helper: Check if eigenvectors are orthonormal (X^H * X = I)
  subroutine check_orthonormality(evec, binfo, mpiglobal, result, tol)
    !> evec:Eigenvector matrix
    type(dzmat), intent(in) :: evec
    !> binfo:BLACS context
    type(blacsinfo), intent(in) :: binfo
    !> mpiglobal:MPI environment
    type(mpiinfo), intent(inout) :: mpiglobal
    !> result:True if orthonormal
    logical, intent(out) :: result
    !> tol:Tolerance
    real(dp), intent(in) :: tol
    
    !> tmp:Temporary matrix for X^H * X
    type(dzmat) :: tmp
    !> gram:Global Gram matrix (on root)
    complex(dp), allocatable :: gram(:,:)
    !> identity:Identity matrix
    complex(dp), allocatable :: identity(:,:)
    !> result_arr:Array for MPI broadcast
    logical :: result_arr(1)
    
    ! X^H * X has size (ncols x ncols), not (nrows x ncols)
    call new_dzmat(tmp, evec%ncols, evec%ncols, binfo)
    
    ! Compute X^H * X using dzmatmult
    call dzmatmult(evec, evec, tmp, transa='C', transb='N')
    
    ! Gather to root for comparison
    if (binfo%mpi%rank == 0) then
      allocate(gram(evec%ncols, evec%ncols))
      allocate(identity(evec%ncols, evec%ncols))
      identity = identity_complex_dp(evec%ncols)
    else
      allocate(gram(1,1))
      allocate(identity(1,1))
    end if
    
    call dzmat_send2global_root(gram, tmp, binfo)
    
    result = .false.
    if (binfo%mpi%rank == 0) then
      result = all_close(gram, identity, tol)
    end if
    
    result_arr(1) = result
    call xmpi_bcast(mpiglobal, result_arr)
    result = result_arr(1)
    
    call del_dzmat(tmp)
    if(allocated(gram)) deallocate(gram)
    if(allocated(identity)) deallocate(identity)
  end subroutine check_orthonormality

  !> Helper: Check eigenvalue equation (H * X = Lambda * X)
  subroutine check_eigenvalue_equation(ham, evec, eval, binfo, mpiglobal, result, tol)
    !> ham:Hamiltonian matrix
    type(dzmat), intent(in) :: ham
    !> evec:Eigenvector matrix
    type(dzmat), intent(in) :: evec
    !> eval:Eigenvalues
    real(dp), intent(in) :: eval(:)
    !> binfo:BLACS context
    type(blacsinfo), intent(in) :: binfo
    !> mpiglobal:MPI environment
    type(mpiinfo), intent(inout) :: mpiglobal
    !> result:True if equation satisfied
    logical, intent(out) :: result
    !> tol:Tolerance
    real(dp), intent(in) :: tol
    
    !> tmp:H * X
    type(dzmat) :: tmp
    !> lhs,rhs:Left/right hand side (global, on root)
    complex(dp), allocatable :: lhs(:,:), rhs(:,:)
    !> i:Loop index
    integer :: i
    !> residual:Residual norm
    real(dp) :: residual
    !> result_arr:Array for MPI broadcast
    logical :: result_arr(1)
    
    call new_dzmat(tmp, ham%nrows, evec%ncols, binfo)
    
    ! Compute H * X using dzmatmult
    call dzmatmult(ham, evec, tmp)
    
    ! Gather to root for comparison
    if (binfo%mpi%rank == 0) then
      allocate(lhs(ham%nrows, evec%ncols))
      allocate(rhs(ham%nrows, evec%ncols))
    else
      allocate(lhs(1,1))
      allocate(rhs(1,1))
    end if
    
    call dzmat_send2global_root(lhs, tmp, binfo)
    call dzmat_send2global_root(rhs, evec, binfo)
    
    ! Compute residual: ||H*X - X*Lambda||_F / ||H*X||_F
    result = .false.
    if (binfo%mpi%rank == 0) then
      do i = 1, size(eval)
        rhs(:,i) = lhs(:,i) - eval(i) * rhs(:,i)
      end do
      residual = sqrt(sum(abs(rhs)**2)) / max(sqrt(sum(abs(lhs)**2)), 1.0e-10_dp)
      result = (residual < tol)
    end if
    
    result_arr(1) = result
    call xmpi_bcast(mpiglobal, result_arr)
    result = result_arr(1)
    
    call del_dzmat(tmp)
    if(allocated(lhs)) deallocate(lhs)
    if(allocated(rhs)) deallocate(rhs)
  end subroutine check_eigenvalue_equation

end module herm_eigensolver_test

