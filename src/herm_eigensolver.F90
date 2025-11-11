module herm_eigensolver
#ifdef _ELPA_
  use elpa
#endif
  use modmpi, only: terminate, terminate_if_false
  use m_hesolver
  use modscl
  use precision
  use iso_fortran_env, only: error_unit
  
  implicit none

  private
  public :: he_eigensolver_wrapper, lapack_eigensolver
#ifdef SCAL
  public :: scalapack_eigensolver
#endif
#ifdef _ELPA_
  public :: elpa_eigensolver
  integer(i32), parameter :: elpa_api_version = 20250131
#endif

contains

  !> Wrapper for Hermitian eigensolvers (ELPA, ScaLAPACK, or LAPACK)
  !> Automatically selects the appropriate eigensolver based on matrix
  !> distribution and compile-time options. Falls back to LAPACK if
  !> distributed solvers are not available.
  subroutine he_eigensolver_wrapper(ham, eval, binfo, evec, i1, i2, v1, v2, found, eecs)
    !> ham:Hermitian matrix to diagonalize
    type(dzmat), intent(inout) :: ham
    !> eval:Eigenvalues (output)
    real(dp), intent(inout) :: eval(:)
    !> binfo:BLACS context information
    type(blacsinfo), intent(in) :: binfo
    !> evec:Eigenvectors (optional)
    type(dzmat), intent(inout), optional :: evec
    !> i1,i2:Index bounds for eigenvalue subset (optional)
    integer(i32), intent(in), optional :: i1, i2
    !> v1,v2:Value bounds for eigenvalue subset (optional)
    real(dp), intent(in), optional :: v1, v2
    !> found:Number of eigenvalues found (optional)
    integer(i32), intent(out), optional :: found
    !> eecs:Extra parameter for eigensolvers (optional)
    integer(i32), intent(in), optional :: eecs

    logical :: distributed, sane

    if (present(evec)) then
      distributed = ham%isdistributed .and. evec%isdistributed
      sane = (ham%isdistributed .eqv. evec%isdistributed)
    else
      distributed = ham%isdistributed
      sane = .true.
    end if

    call terminate_if_false(sane, &
         'Error(he_eigensolver_wrapper): Inconsistent matrix distribution')

    if (distributed) then
#ifdef _ELPA_
      call elpa_eigensolver(ham, eval, binfo, evec, i1, i2, v1, v2, found, eecs)
#elif defined(SCAL)
      call scalapack_eigensolver(ham, eval, binfo, evec, i1, i2, v1, v2, found, eecs)
#else
      ! Fall back to LAPACK if no distributed solver is available
      call lapack_eigensolver(ham, eval, binfo, evec, i1, i2, v1, v2, found, eecs)
#endif
    else
      call lapack_eigensolver(ham, eval, binfo, evec, i1, i2, v1, v2, found, eecs)
    end if

  end subroutine he_eigensolver_wrapper

  !> Solves a Hermitian eigenvalue problem using ELPA library
  !> Computes eigenvalues and optionally eigenvectors using ELPA
  !> (Eigenvalue SoLvers for Petaflop-Applications).
  subroutine elpa_eigensolver(ham, eval, binfo, evec, i1, i2, v1, v2, found, eecs)
#ifdef _OPENMP
    use omp_lib
#endif
    !> ham:Hermitian matrix to diagonalize
    type(dzmat), intent(inout) :: ham
    !> eval:Eigenvalues (output)
    real(dp), intent(inout) :: eval(:)
    !> binfo:BLACS context information
    type(blacsinfo), intent(in) :: binfo
    !> evec:Eigenvectors (optional)
    type(dzmat), intent(inout), optional :: evec
    !> i1,i2:Index bounds for eigenvalue subset (optional)
    integer(i32), intent(in), optional :: i1, i2
    !> v1,v2:Value bounds for eigenvalue subset (optional)
    real(dp), intent(in), optional :: v1, v2
    !> found:Number of eigenvalues found (optional)
    integer(i32), intent(out), optional :: found
    !> eecs:Extra parameter (currently unused)
    integer(i32), intent(in), optional :: eecs

#ifdef _ELPA_
    class(elpa_t), pointer :: eh
    integer(i32)          :: na, nev
    integer(i32)          :: error, success
    integer(i32)          :: loc_il, loc_iu
    logical               :: subset_i, subset_v
    logical               :: evalsonly
    integer(i32)          :: omp_threads
    ! Variables for subset selection
    integer(i32)          :: k, j
    integer(i32), allocatable :: idx(:)
    logical, allocatable  :: mask(:)
    real(dp), allocatable :: eval_all(:)

    evalsonly = .not. present(evec)

    call terminate_if_false(ham%context == binfo%context, 'Error(elpa_eigensolver): &
         context mismatch: ham%context /= binfo%context.')
    if (present(evec)) then
      call terminate_if_false(evec%context == binfo%context, 'Error(elpa_eigensolver): &
           context mismatch: evec%context /= binfo%context.')
    end if
    na        = ham%nrows
    subset_i  = present(i1) .or. present(i2)
    subset_v  = present(v1) .and. present(v2)
    loc_il = 1 ; loc_iu = na
    if (subset_i) then
      if (present(i1)) loc_il = i1
      if (present(i2)) loc_iu = i2
    end if
    nev = loc_iu - loc_il + 1

    call terminate_if_false(elpa_init(elpa_api_version) == ELPA_OK, &
         'Error(elpa_eigensolver): ELPA API version mismatch')

    eh => elpa_allocate(error)
    call terminate_if_false(error == ELPA_OK, 'Error(elpa_eigensolver): ELPA allocate failed')

    ! Set OpenMP thread count
#ifdef _OPENMP
    omp_threads = omp_get_max_threads()
#else
    omp_threads = 1
#endif

    call eh%set("na",              na,             error)
    call eh%set("nev",             nev,            error)
    call eh%set("local_nrows",     ham%nrows_loc,  error)
    call eh%set("local_ncols",     ham%ncols_loc,  error)
    call eh%set("nblk",            ham%mblck,      error)
    call eh%set("mpi_comm_parent", binfo%mpi%comm, error)
    call eh%set("process_row",     binfo%myprow,   error)
    call eh%set("process_col",     binfo%mypcol,   error)
    call eh%set("omp_threads", omp_threads, error)

    success = eh%setup()
    call terminate_if_false(success == ELPA_OK, 'Error(elpa_eigensolver): ELPA setup failed')

    call eh%set("solver", ELPA_SOLVER_2STAGE, error)
    call terminate_if_false(error == ELPA_OK, 'Error(elpa_eigensolver): setting solver failed')
    
    if (evalsonly) then
      call eh%eigenvalues(ham%za, eval, error)
      call terminate_if_false(error == ELPA_OK, 'Error(elpa_eigensolver): computing eigenvalues failed')
    else
      call eh%eigenvectors(ham%za, eval, evec%za, error)
      call terminate_if_false(error == ELPA_OK, 'Error(elpa_eigensolver): computing eigenvectors failed')
    end if
    
    ! Handle subset selection (safe form to avoid shape mismatch)
    if (subset_i) then
      k = loc_iu - loc_il + 1
      eval(1:k) = eval(loc_il:loc_iu)
      if (present(evec)) then
        do j = 1, k
          evec%za(:, j) = evec%za(:, loc_il + j - 1)
        end do
      end if
      if (present(found)) found = k
    else if (subset_v) then
      allocate(eval_all(size(eval)))
      eval_all = eval
      allocate(mask(size(eval_all)))
      mask = (eval_all >= v1 .and. eval_all <= v2)
      k = count(mask)
      if (k > 0) then
        allocate(idx(k))
        idx = pack([(j, j=1,size(eval_all))], mask)
        eval(1:k) = eval_all(idx)
        if (present(evec)) then
          do j = 1, k
            evec%za(:, j) = evec%za(:, idx(j))
          end do
        end if
        deallocate(idx)
      end if
      if (present(found)) found = k
      deallocate(eval_all, mask)
    else
      if (present(found)) found = size(eval)
    end if
    call elpa_deallocate(eh, error)
    call elpa_uninit()
#else
    call terminate("Error(elpa_eigensolver): elpa_eigensolver called but exciting is &
         not linked to ELPA.")
#endif
  end subroutine elpa_eigensolver

  !> Solves a Hermitian eigenvalue problem using ScaLAPACK library
  !> Computes eigenvalues and optionally eigenvectors using ScaLAPACK's pzheevx routine.
  subroutine scalapack_eigensolver(ham, eval, binfo, evec, i1, i2, v1, v2, found, eecs)
    !> ham:Hermitian matrix to diagonalize
    type(dzmat), intent(inout) :: ham
    !> eval:Eigenvalues (output)
    real(dp), intent(inout) :: eval(:)
    !> binfo:BLACS context information
    type(blacsinfo), intent(in) :: binfo
    !> evec:Eigenvectors (optional)
    type(dzmat), intent(inout), optional :: evec
    !> i1,i2:Index bounds for eigenvalue subset (optional)
    integer(i32), intent(in), optional :: i1, i2
    !> v1,v2:Value bounds for eigenvalue subset (optional)
    real(dp), intent(in), optional :: v1, v2
    !> found:Number of eigenvalues found (optional)
    integer(i32), intent(out), optional :: found
    !> eecs:Expected eigenvalue cluster size (optional)
    integer(i32), intent(in), optional :: eecs

#ifdef SCAL
    logical :: evalsonly
    character(1) :: rangechar, jobzchar
    integer(i32) :: info
    integer(i32) :: lwork, lrwork, liwork
    integer(i32) :: il, iu
    real(dp) :: abstol, vl, vu
    complex(dp), allocatable :: work(:)
    real(dp), allocatable :: rwork(:)
    integer(i32), allocatable :: iwork(:)
    real(dp) :: orfac
    integer(i32) :: ia, ja, iz, jz
    integer(i32) :: nevalfound, nevecfound
    real(dp), allocatable :: gap(:)
    integer(i32), allocatable :: iclustr(:), ifail(:)
    integer(i32) :: clusterguess
    type(dzmat) :: evecdummy
    real(dp), external :: pdlamch

    evalsonly = .not. present(evec)

    if(.not. evalsonly) then
      call terminate_if_false(ham%context == evec%context .and. &
           ham%context == binfo%context .and. evec%context == binfo%context, &
           'Error(scalapack_eigensolver): ham, evec and binfo have differing contexts.')
    end if
    call terminate_if_false(.not. ((present(i1) .or. present(i2)) .and. &
         (present(v1) .or. present(v2))), &
         'Error(scalapack_eigensolver): I and V specified.')
    call terminate_if_false(.not. (present(v1) .and. .not. present(v2)), &
         'Error(scalapack_eigensolver): Specify whole interval.')

    if(present(i1) .or. present(i2)) then
      il = 1
      iu = ham%nrows
      rangechar = 'I'
      if(present(i1)) il = i1
      if(present(i2)) iu = i2
      call terminate_if_false(il >= 1 .and. iu >= 1, &
           'Error(scalapack_eigensolver): iu and il need to be positive.')
      call terminate_if_false(il <= iu, &
           'Error(scalapack_eigensolver): il > iu.')
      call terminate_if_false(iu <= ham%nrows, &
           'Error(scalapack_eigensolver): iu > nrows.')
    end if

    if(present(v1) .and. present(v2)) then
      rangechar = 'V'
      vl = v1
      vu = v2
      call terminate_if_false(vl <= vu, &
           'Error(scalapack_eigensolver): vl > vu')
    end if

    if(.not. (present(i1) .or. present(i2)) .and. .not. present(v1)) then
      rangechar = 'A'
    end if

    if(evalsonly) then
      jobzchar = 'N'
    else
      jobzchar = 'V'
    end if

    orfac = 1.0e-9_dp
    abstol = 2.0_dp * pdlamch(binfo%context, 'S')

    allocate(ifail(ham%nrows))
    allocate(iclustr(2*binfo%nprows*binfo%npcols))
    allocate(gap(binfo%nprows*binfo%npcols))

    ia = 1
    ja = 1
    iz = 1
    jz = 1

    clusterguess = ceiling(real(ham%nrows, dp) * 0.2_dp)
    if(present(eecs)) then
      if(eecs >= 1) then
        clusterguess = eecs
      end if
    end if

    call workspacequery(jobzchar, rangechar)
    allocate(work(lwork), rwork(lrwork), iwork(liwork))

    if(evalsonly) then
      call new_dzmat(evecdummy,1,1,binfo)
      call pzheevx(jobzchar, rangechar, 'U', ham%nrows, ham%za, ia, ja, ham%desc,&
        & vl, vu, il, iu, abstol, nevalfound, nevecfound, eval, orfac,&
        & evecdummy%za, iz, jz, evecdummy%desc, work, lwork, rwork, lrwork, iwork, liwork,&
        & ifail, iclustr, gap, info)
      call del_dzmat(evecdummy)
    else
      call pzheevx(jobzchar, rangechar, 'U', ham%nrows, ham%za, ia, ja, ham%desc,&
        & vl, vu, il, iu, abstol, nevalfound, nevecfound, eval, orfac,&
        & evec%za, iz, jz, evec%desc, work, lwork, rwork, lrwork, iwork, liwork,&
        & ifail, iclustr, gap, info)
    end if

    if(present(found)) then
      found = nevalfound
    end if

    call errorinspect(binfo, info)

    deallocate(ifail)
    deallocate(iclustr)
    deallocate(gap)
    deallocate(work, rwork, iwork)

#else
    call terminate("Error(scalapack_eigensolver): scalapack_eigensolver called but exciting &
         is not linked to ScaLAPACK.")
#endif
  contains

#ifdef SCAL
    subroutine workspacequery(jobtype, rangetype)
      character(1), intent(in) :: jobtype, rangetype

      integer(i32) :: nhetrd_lwork, anb, nps, n
      integer(i32), external :: pjlaenv, iceil, numroc
      integer(i32) :: sqnpc

      n = ham%nrows
      anb = pjlaenv(binfo%context, 3, 'PZHETTRD', 'L', 0, 0, 0, 0)
      sqnpc = int(sqrt(real(binfo%nprows*binfo%npcols, dp)), i32)
      nps = max(numroc(n, 1, 0, 0, sqnpc), 2*anb)
      nhetrd_lwork = n + 2*(anb+1)*(4*nps+2)+(nps+1)*nps

      lwork=-1
      lrwork=-1
      liwork=-1

      allocate(work(3), rwork(3), iwork(3))

      if(jobtype == 'N') then
        call new_dzmat(evecdummy,1,1,binfo)
        call pzheevx(jobtype, rangetype, 'U', ham%nrows, ham%za, ia, ja, ham%desc,&
          & vl, vu, il, iu, abstol, nevalfound, nevecfound, eval, orfac,&
          & evecdummy%za, iz, jz, evecdummy%desc, work, lwork, rwork, lrwork, iwork, liwork,&
          & ifail, iclustr, gap, info)
        call del_dzmat(evecdummy)
      else
        call pzheevx(jobtype, rangetype, 'U', ham%nrows, ham%za, ia, ja, ham%desc,&
          & vl, vu, il, iu, abstol, nevalfound, nevecfound, eval, orfac,&
          & evec%za, iz, jz, evec%desc, work, lwork, rwork, lrwork, iwork, liwork,&
          & ifail, iclustr, gap, info)
      end if

      lwork=max(int(work(1)), nhetrd_lwork)
      if(jobtype == 'N') then
        lrwork=int(rwork(1))
      else
        lrwork=int(rwork(1)) + (clusterguess-1)*n
      end if
      liwork=int(iwork(1))

      deallocate(work, rwork, iwork)
    end subroutine workspacequery

    subroutine errorinspect(binfo, ierror)
      type(blacsinfo), intent(in) :: binfo
      integer(i32), intent(in) :: ierror
      integer(i32) :: i, maxcs, tmp

      if (ierror == 0) return

      if (binfo%mpi%rank == 0) then
        write(error_unit,'("Error(scalapack_eigensolver): pzheevx returned non-zero info:", i6)') ierror
        
        if( ierror < 0) then
          write(error_unit,'("Error(scalapack_eigensolver) cause: Invalid input")')
        else if(mod(ierror,2) /= 0) then
          write(error_unit,'("Error(scalapack_eigensolver) cause: Eigenvectors not converged")')
          write(error_unit,'("scalapack_eigensolver ifail")')
          write(error_unit,'(I8)') ifail
        else if(mod(ierror/2,2) /= 0) then
          maxcs = 0
          do i = 1, size(iclustr)-1
            tmp = iclustr(i+1) - iclustr(i)
            if(tmp > 0) then
              maxcs = max(tmp, maxcs)
            else
              exit
            end if
          end do
          i = i-1
          write(error_unit,'("Warning(scalapack_eigensolver) cause: Reorthogonalization failed,&
            & insufficient workspace. There are", i8," clusters of eigenvalues&
            & and the largest one has size ", i8,". &
            & Increase input%xs%bse%eecs to ", i8," to guarantee orthogonal eigenvectors.&
            & (usually the results are still good)")')&
            & i, maxcs, maxcs
          write(error_unit,'("scalapack_eigensolver iclustr:")')
          write(error_unit,'(I8)') iclustr
        else if(mod(ierror/4,2) /= 0) then
          write(error_unit,'("Error(scalapack_eigensolver) cause:&
            & Not all eigenvectors computed, insufficient workspace.")')
        else if(mod(ierror/8,2) /= 0) then
          write(error_unit,'("Error(scalapack_eigensolver) cause:&
            & Eigenvalue computation failed")')
        end if
      end if
      
      call terminate_if_false(.false., 'Error(scalapack_eigensolver): &
           pzheevx returned non-zero info')
    end subroutine errorinspect
#endif

  end subroutine scalapack_eigensolver

  !> Solves a Hermitian eigenvalue problem using LAPACK library
  !> Computes eigenvalues and optionally eigenvectors using LAPACK routines (non-distributed).
  subroutine lapack_eigensolver(ham, eval, binfo, evec, i1, i2, v1, v2, found, eecs)
    !> ham:Hermitian matrix to diagonalize
    type(dzmat), intent(inout) :: ham
    !> eval:Eigenvalues (output)
    real(dp), intent(inout) :: eval(:)
    !> binfo:BLACS context information
    type(blacsinfo), intent(in) :: binfo
    !> evec:Eigenvectors (optional)
    type(dzmat), intent(inout), optional :: evec
    !> i1,i2:Index bounds for eigenvalue subset (optional)
    integer(i32), intent(in), optional :: i1, i2
    !> v1,v2:Value bounds for eigenvalue subset (optional)
    real(dp), intent(in), optional :: v1, v2
    !> found:Number of eigenvalues found (optional)
    integer(i32), intent(out), optional :: found
    !> eecs:Extra parameter (currently unused)
    integer(i32), intent(in), optional :: eecs

    if (present(evec)) then
      call hesolver(ham%za, eval, evec%za, i1, i2, v1, v2, found)
    else
      call hesolver(ham%za, eval, i1=i1, i2=i2, v1=v1, v2=v2, found=found)
    end if
    
  end subroutine lapack_eigensolver

end module herm_eigensolver

