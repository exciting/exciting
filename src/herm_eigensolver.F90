module herm_eigensolver
#ifdef _ELPA_
  use elpa, only: elpa_t, ELPA_OK, elpa_init, elpa_allocate, ELPA_SOLVER_1STAGE,&
    & ELPA_SOLVER_2STAGE, elpa_deallocate, elpa_uninit
#endif
  use modmpi, only: terminate, terminate_if_false
  use m_hesolver
  use modscl
  use precision
  use to_char_conversion, only: to_char
  use iso_fortran_env, only: error_unit
  
  implicit none

  private
  public :: &
          lapack_eigensolver, &
          scalapack_eigensolver_pzheevx, &
          scalapack_eigensolver_pzheevd, &
          elpa_eigensolver
  ! Used if the ScaLAPACK routine PZHEEVX is requested to solve the BSE Hamiltonian.
  ! Specifies which eigenvectors should be reorthogonalized. See the ScaLAPACK documentation for more details.
  ! A negative value (the default) will use the ScaLAPACK default (which is the legacy behavior).
  real(dp), parameter :: scalapack_pzheevx_orfac = -1.0
  integer(i32), parameter :: elpa_api_version = 20250131

contains
  !> Solves a Hermitian eigenvalue problem using ELPA library
  !> Computes eigenvalues and optionally eigenvectors using ELPA
  !> (Eigenvalue SoLvers for Petaflop-Applications).
  subroutine elpa_eigensolver(ham, eval, binfo, elpa_solver, evec, i2, found)
#ifdef _OPENMP
    use omp_lib
#endif
    !> ham:Hermitian matrix to diagonalize
    type(dzmat), intent(inout) :: ham
    !> eval:Eigenvalues (output)
    real(dp), intent(inout) :: eval(:)
    !> binfo:BLACS context information
    type(blacsinfo), intent(in) :: binfo
    !> which ELPA solver to use, must be '1' or '2'
    character(1), intent(in) :: elpa_solver
    !> evec:Eigenvectors (optional)
    type(dzmat), intent(inout), optional :: evec
    !> i2:Upper index bound for eigenvalue subset (optional)
    integer(i32), intent(in), optional :: i2
    !> found:Number of eigenvalues found (optional)
    integer(i32), intent(out), optional :: found

#ifdef _ELPA_
    class(elpa_t), pointer :: eh
    integer(i32)          :: na, nev
    integer(i32)          :: error, success
    logical               :: evalsonly
    integer(i32)          :: omp_threads

    evalsonly = .not. present(evec)

    call terminate_if_false(ham%context == binfo%context, 'Error(elpa_eigensolver): &
         context mismatch: ham%context /= binfo%context.')
    if (present(evec)) then
      call terminate_if_false(evec%context == binfo%context, 'Error(elpa_eigensolver): &
           context mismatch: evec%context /= binfo%context.')
    end if

    na = ham%nrows
    nev = na
    if (present(i2)) nev = i2
    call terminate_if_false(1 <= nev .and. nev <= na, 'Error(elpa_eigensolver): Upper index selection out of range.')

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
    call eh%set("omp_threads",     omp_threads,    error)

    success = eh%setup()
    call terminate_if_false(success == ELPA_OK, 'Error(elpa_eigensolver): ELPA setup failed')
    if (elpa_solver == '1') then
      call eh%set("solver", ELPA_SOLVER_1STAGE, error)
    else if (elpa_solver == '2') then
      call eh%set("solver", ELPA_SOLVER_2STAGE, error)
    else
      call terminate('Error(elpa_eigensolver): Solver name not recognized, use "1" or "2".')
    end if
    call terminate_if_false(error == ELPA_OK, 'Error(elpa_eigensolver): setting solver failed')

    if (evalsonly) then
      call eh%eigenvalues(ham%za, eval, error)
      call terminate_if_false(error == ELPA_OK, 'Error(elpa_eigensolver): computing eigenvalues failed')
    else
      call eh%eigenvectors(ham%za, eval, evec%za, error)
      call terminate_if_false(error == ELPA_OK, 'Error(elpa_eigensolver): computing eigenvectors failed')
    end if

    if (present(found)) found = nev

    call elpa_deallocate(eh, error)
    call elpa_uninit()
#else
    call terminate("Error(elpa_eigensolver): elpa_eigensolver called but exciting is &
         not linked to ELPA.")
#endif
  end subroutine elpa_eigensolver

  !> Solves a Hermitian eigenvalue problem using ScaLAPACK library
  !> Computes eigenvalues and optionally eigenvectors using ScaLAPACK's pzheevx routine.
  subroutine scalapack_eigensolver_pzheevx(ham, eval, binfo, evec, i1, i2, v1, v2, found, eecs)
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
           'Error(scalapack_eigensolver_pzheevx): ham, evec and binfo have different contexts.')
    end if
    call terminate_if_false(.not. ((present(i1) .or. present(i2)) .and. &
         (present(v1) .or. present(v2))), &
         'Error(scalapack_eigensolver_pzheevx): I and V specified.')
    call terminate_if_false(.not. (present(v1) .and. .not. present(v2)), &
         'Error(scalapack_eigensolver_pzheevx): Specify whole interval.')

    if(present(i1) .or. present(i2)) then
      il = 1
      iu = ham%nrows
      rangechar = 'I'
      if(present(i1)) il = i1
      if(present(i2)) iu = i2
      call terminate_if_false(il >= 1 .and. iu >= 1, &
           'Error(scalapack_eigensolver_pzheevx): iu and il need to be positive.')
      call terminate_if_false(il <= iu, &
           'Error(scalapack_eigensolver_pzheevx): il > iu.')
      call terminate_if_false(iu <= ham%nrows, &
           'Error(scalapack_eigensolver_pzheevx): iu > nrows.')
    end if

    if(present(v1) .and. present(v2)) then
      rangechar = 'V'
      vl = v1
      vu = v2
      call terminate_if_false(vl <= vu, &
           'Error(scalapack_eigensolver_pzheevx): vl > vu')
    end if

    if(.not. (present(i1) .or. present(i2)) .and. .not. present(v1)) then
      rangechar = 'A'
    end if

    if(evalsonly) then
      jobzchar = 'N'
    else
      jobzchar = 'V'
    end if

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

    call workspacequery_pzheevx(jobzchar, rangechar)
    allocate(work(lwork), rwork(lrwork), iwork(liwork))

    if(evalsonly) then
      call new_dzmat(evecdummy,1,1,binfo)
      call pzheevx(jobzchar, rangechar, 'U', ham%nrows, ham%za, ia, ja, ham%desc,&
        & vl, vu, il, iu, abstol, nevalfound, nevecfound, eval, scalapack_pzheevx_orfac,&
        & evecdummy%za, iz, jz, evecdummy%desc, work, lwork, rwork, lrwork, iwork, liwork,&
        & ifail, iclustr, gap, info)
      call del_dzmat(evecdummy)
    else
      call pzheevx(jobzchar, rangechar, 'U', ham%nrows, ham%za, ia, ja, ham%desc,&
        & vl, vu, il, iu, abstol, nevalfound, nevecfound, eval, scalapack_pzheevx_orfac,&
        & evec%za, iz, jz, evec%desc, work, lwork, rwork, lrwork, iwork, liwork,&
        & ifail, iclustr, gap, info)
    end if

    if(present(found)) then
      found = nevalfound
    end if

    call errorinspect_pzheevx(binfo, info)

    deallocate(ifail)
    deallocate(iclustr)
    deallocate(gap)
    deallocate(work, rwork, iwork)

#else
    call terminate("Error(scalapack_eigensolver_pzheevx): scalapack_eigensolver_pzheevx called but exciting &
         is not linked to ScaLAPACK.")
#endif
  contains

#ifdef SCAL
    !> Call ScaLAPACK to obtain the optimal workspace sizes
    subroutine workspacequery_pzheevx(jobtype, rangetype)
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
          & vl, vu, il, iu, abstol, nevalfound, nevecfound, eval, scalapack_pzheevx_orfac,&
          & evecdummy%za, iz, jz, evecdummy%desc, work, lwork, rwork, lrwork, iwork, liwork,&
          & ifail, iclustr, gap, info)
        call del_dzmat(evecdummy)
      else
        call pzheevx(jobtype, rangetype, 'U', ham%nrows, ham%za, ia, ja, ham%desc,&
          & vl, vu, il, iu, abstol, nevalfound, nevecfound, eval, scalapack_pzheevx_orfac,&
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
    end subroutine workspacequery_pzheevx

    !> Analyse the return code from ScaLAPACK, kill the code if a fatal code is detected, print error messages
    subroutine errorinspect_pzheevx(binfo, ierror)
      type(blacsinfo), intent(in) :: binfo
      integer(i32), intent(in) :: ierror
      integer(i32) :: i, maxcs, tmp

      if (ierror == 0) return

      write(error_unit,'("Error(scalapack_eigensolver_pzheevx): pzheevx returned non-zero info:", i6, " on rank:", i6)')&
        & ierror, binfo%mpi%rank

      if( ierror < 0) then
        write(error_unit,'("Error(scalapack_eigensolver_pzheevx): ", a)') &
          "pzheevx argument no. " // to_char(ierror) // " is invalid."
      else if(mod(ierror,2) /= 0) then
        write(error_unit,'("Error(scalapack_eigensolver_pzheevx): ", a)') &
          "Eigenvectors did not converge. pzheevx returned ifail = " // to_char(ifail)
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

        write(error_unit,'("Warning(scalapack_eigensolver_pzheevx): ", a)') &
          "Reorthogonalization failed because of insufficient workspace. There are " &
       // to_char(i) // " clusters of eigenvalues. The largest has size " // to_char(maxcs) // "." &
       // "To guarantee orthogonal eigenvectors set <input><xs><bse><eecs> at least to " // to_char(maxcs) // "." &
       // "Usually the results are still good."
        write(error_unit,'("scalapack_eigensolver_pzheevx iclustr:")')
        write(error_unit,'(I8)') iclustr
        return
      else if(mod(ierror/4,2) /= 0) then
        write(error_unit,'("Error(scalapack_eigensolver_pzheevx): ", a)') &
          "Could not compute all eigenvectors because of insufficient workspace."
      else if(mod(ierror/8,2) /= 0) then
        write(error_unit,'("Error(scalapack_eigensolver_pzheevx): ", a)') &
          "Eigenvalue computation failed"
      end if

      call terminate_if_false(.false., 'Error(scalapack_eigensolver_pzheevx): &
           pzheevx returned non-zero info')
    end subroutine errorinspect_pzheevx
#endif

  end subroutine scalapack_eigensolver_pzheevx

  !> Solves a Hermitian eigenvalue problem using ScaLAPACK library
  !> Computes eigenvalues and eigenvectors using ScaLAPACK's pzheevd routine.
  subroutine scalapack_eigensolver_pzheevd(ham, eval, binfo, evec, found)
    !> ham:Hermitian matrix to diagonalize
    type(dzmat), intent(inout) :: ham
    !> eval:Eigenvalues (output)
    real(dp), intent(inout) :: eval(:)
    !> binfo:BLACS context information
    type(blacsinfo), intent(in) :: binfo
    !> evec:Eigenvectors
    type(dzmat), intent(inout) :: evec
    !> found:Number of eigenvalues found (optional)
    integer(i32), intent(out), optional :: found

#ifdef SCAL
    integer(i32) :: info
    integer(i32) :: lwork, lrwork, liwork
    complex(dp), allocatable :: work(:)
    real(dp), allocatable :: rwork(:)
    integer(i32), allocatable :: iwork(:)
    integer(i32) :: ia, ja, iz, jz
    character(16) :: info_str

    call terminate_if_false(ham%context == evec%context .and. &
         ham%context == binfo%context .and. evec%context == binfo%context, &
         'Error(scalapack_eigensolver_pzheevd): ham, evec and binfo have different contexts.')

    ia = 1
    ja = 1
    iz = 1
    jz = 1

    lwork=-1
    lrwork=-1
    liwork=-1

    allocate(work(3), rwork(3), iwork(3))

    ! workspace query
    call pzheevd('V', 'U', ham%nrows, ham%za, ia, ja, ham%desc,&
      & eval,&
      & evec%za, iz, jz, evec%desc, work, lwork, rwork, lrwork, iwork, liwork,&
      & info)

    write(info_str,'(I5)') info
    call terminate_if_false(info == 0, 'Error(scalapack_eigensolver_pzheevd): &
         pzheevd workspacequery returned non-zero info: '//info_str)

    lwork  = int(real(work(1)))
    lrwork = int(rwork(1))
    liwork = iwork(1)

    deallocate(work, rwork, iwork)
    allocate(work(lwork), rwork(lrwork), iwork(liwork))

    call pzheevd('V', 'U', ham%nrows, ham%za, ia, ja, ham%desc, eval,&
      & evec%za, iz, jz, evec%desc, work, lwork, rwork, lrwork, iwork, liwork, info)

    call errorinspect_pzheevd(binfo, info)

    if (present(found)) found = size(eval)

    deallocate(work, rwork, iwork)

#else
    call terminate("Error(scalapack_eigensolver_pzheevd): scalapack_eigensolver_pzheevd called but exciting &
         is not linked to ScaLAPACK.")
#endif
    end subroutine scalapack_eigensolver_pzheevd

#ifdef SCAL
    !> Analyse the return code from ScaLAPACK, kill the code if a fatal code is detected, print error messages
    subroutine errorinspect_pzheevd(binfo, ierror)
      type(blacsinfo), intent(in) :: binfo
      integer(i32), intent(in) :: ierror

      if (ierror == 0) return

    write(error_unit,'("Error(scalapack_eigensolver_pzheevd): pzheevd returned non-zero info:", i6, " on rank:", i6)')&
        & ierror, binfo%mpi%rank

    if( ierror < 0) then
      write(error_unit,'("Error(scalapack_eigensolver_pzheevd): ", a)') &
          "pzheevd argument no. " // to_char(ierror) // " is invalid."
    else
      write(error_unit,'("Error(scalapack_eigensolver_pzheevd): ", a)') &
          "Eigenvalue no. " // to_char(ierror) // " did not converged."
    end if

      call terminate_if_false(.false., 'Error(scalapack_eigensolver_pzheevd): &
           pzheevd returned non-zero info')
    end subroutine errorinspect_pzheevd
#endif

  !> Solves a Hermitian eigenvalue problem using LAPACK library
  !> Computes eigenvalues and optionally eigenvectors using LAPACK routines (non-distributed).
  subroutine lapack_eigensolver(ham, eval, binfo, evec, i1, i2, v1, v2, found)
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

    if (present(evec)) then
      call hesolver(ham%za, eval, evec%za, i1, i2, v1, v2, found)
    else
      call hesolver(ham%za, eval, i1=i1, i2=i2, v1=v1, v2=v2, found=found)
    end if
    
  end subroutine lapack_eigensolver

end module herm_eigensolver

