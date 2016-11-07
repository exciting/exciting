module m_dhesolver
  use modmpi 
  use modscl
  use m_hesolver

  implicit none

  contains

    !BOP
    ! !ROUTINE: dhesolver
    ! !INTERFACE:
    subroutine dhesolver(ham, evec, eval, binfo, i1, i2, v1, v2, found, eecs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   type(blacsinfo) :: binfo       ! Info type describing the BLACS grid
    !   integer(4), optional :: i1, i2 ! Index range of eigen-solutions 
    !   real(8), optional :: v1, v2    ! Range of eigenvalues to search for
    !   integer(4), optional :: eecs   ! Estimate for eigenvalue clustering.
    !                                  ! Apriori not known but needed for propper
    !                                  ! orthogonalization of eigenvectors (default = 3)
    ! IN/OUT:
    !   type(dzmat) :: ham         ! 2D block cyclic distributed hermitian matrix 
    !   type(dzmat) :: evec        ! 2D block cyclic distributed eigenvector matrix
    !   real(8) :: eval(ham%nrows) ! Real valued eigenvalues in ascending order 
    !                              ! (the first solsize elements are set)
    ! OUT:
    !   integer(4), optional :: found ! How many solutions were found
    !
    ! !DESCRIPTION:
    !   Takes the upper triangular part of an distributed complex matrix matrix,
    !   assumed to be hermitian and finds eigenvalues and eigenvectors using 
    !   the scalapack routine {\tt pzheevx}.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      type(dzmat), intent(inout) :: ham
      type(dzmat), intent(inout) :: evec
      real(8), intent(inout) :: eval(:)
      type(blacsinfo), intent(in) :: binfo
      integer(4), intent(in), optional :: eecs
      integer(4), intent(in), optional :: i1, i2
      real(8), intent(in), optional :: v1, v2
      integer(4), intent(out), optional :: found

      ! Local variables
#ifdef SCAL
      character(1) :: rangechar
      integer(4) :: info
      integer(4) :: lwork, lrwork, liwork
      integer(4) :: il, iu, solsize
      real(8) :: abstol, vl, vu
      complex(8), allocatable :: work(:)
      real(8), allocatable :: rwork(:)
      integer(4), allocatable :: iwork(:)
      real(8) :: orfac
      integer(4) :: ia, ja, iz, jz
      integer(4) :: nevalfound, nevecfound
      real(8), allocatable :: gap(:)
      integer(4), allocatable :: iclustr(:), ifail(:)
      integer(4) :: clusterguess

      ! External functions
      real(8), external :: pdlamch


      ! Distributed EVP
      if(ham%isdistributed .and. evec%isdistributed) then

        ! Sanity checking input
        if(ham%context /= evec%context&
          &.or. ham%context /= binfo%context&
          &.or. evec%context /= binfo%context) then
          if(binfo%mpi%rank == 0) then
            write(*,*) "dhesolver (ERROR):&
              & ham,evec and binfo have differing contexts.",&
              & ham%context, evec%context, binfo%context
            call terminate
          end if
        end if
        if((present(i1) .or. present(i2)) .and. (present(v1) .or. present(v2))) then
          if(binfo%mpi%rank == 0) then
            write(*,*) "dhesolver (ERROR): I and V specivied."
            call terminate
          end if
        end if
        if(present(v1) .and. .not. present(v2)) then 
          if(binfo%mpi%rank == 0) then
            write(*,*) "dhesolver (ERROR): Specify whole interval."
            call terminate
          end if
        end if
        if(present(i1) .or. present(i2)) then
          il = 1
          iu = ham%nrows
          rangechar = 'I'
          if(present(i1)) il = i1
          if(present(i2)) iu = i2
          ! Sanity check
          if(il < 1 .or. iu < 1) then 
            if(binfo%mpi%rank == 0) then
              write(*,*) "dhesolver (ERROR): iu and il need to be positive."
              call terminate
            end if
          end if
          if(il > iu) then 
            if(binfo%mpi%rank == 0) then
              write(*,*) "dhesolver (ERROR): il > iu."
              call terminate
            end if
          end if
          if(iu > ham%nrows) then 
            if(binfo%mpi%rank == 0) then
              write(*,*) "dhesolver (ERROR): iu > nrows."
              write(*,*) "range:", il, iu
              write(*,*) "hamsize:", ham%nrows
              call terminate
            end if
          end if
        end if
        if(present(v1) .and. present(v2)) then
          rangechar = 'V'
          vl = v1
          vu = v2
          if(vl > vu) then
            if(binfo%mpi%rank == 0) then
              write(*,*) "dhesolver (ERROR): vl > vu"
              write(*,*) "range:", vl, vu
              call terminate
            end if
          end if
        end if
        if(.not. (present(i1) .or. present(i2)) .or. .not. present(v1)) then
          rangechar = 'A'
        end if

        ! Tolerance parameter
        ! Ortho-normality criterion (1d-3 is default)
        orfac = 1d0-9
        ! Numeric tolerance optimized for eigenvalues
        abstol = 2.0d0 * pdlamch(binfo%context, 'S') 
        ! Numeric tolerance optimized for orthogonality
        !abstol = pdlamch(binfo%context, 'U') 

        ! Allocate supports
        allocate(ifail(ham%nrows))
        allocate(iclustr(2*binfo%nprows*binfo%npcols))
        allocate(gap(binfo%nprows*binfo%npcols))

        ! Global coordinates of first element of sub-matrix of ham
        ! (We have no sub-selection here we operate on all of ham)
        ia = 1
        ja = 1
        ! Global coordinates of first element of sub-matrix of evec
        ! (We have no sub-selection here we operate on all of evec)
        iz = 1
        jz = 1

        ! Set guess for eigenvalue clustering
        if(present(eecs)) then
          if(eecs >= 1) then
            clusterguess = eecs
          else
            clusterguess = 3
          end if
        else
          clusterguess = 3
        end if
      

        ! Get optimal work array lengths
        call workspacequery(rangechar)
        allocate(work(lwork), rwork(lrwork), iwork(liwork))

        ! Compute
        call pzheevx('V', rangechar, 'U', ham%nrows, ham%za, ia, ja, ham%desc,&
          & vl, vu, il, iu, abstol, nevalfound, nevecfound, eval, orfac,&
          & evec%za, iz, jz, evec%desc, work, lwork, rwork, lrwork, iwork, liwork,&
          & ifail, iclustr, gap, info)

        if(present(found)) then
          found = nevalfound
        end if

        ! Error inspection
        if(info .ne. 0) then
          if(binfo%mpi%rank == 0) then
            write(*, '("ERROR (dhesolver):&
              & pzheevx returned non-zero info:", i6)') info
            call errorinspect(info)
          end if
        end if

        ! Clean up
        deallocate(ifail)
        deallocate(iclustr)
        deallocate(gap)
        deallocate(work, rwork, iwork)

      ! Non distributed EVP
      else if(.not. ham%isdistributed .and. .not. evec%isdistributed) then

        call hesolver(ham%za, evec%za, eval, i1, i2, v1, v2, found)

      ! Mix -> Error
      else 

        write(*, '("dhesolver (ERROR):&
          & Distributed matrix mixed with non distributed:", l, l)')&
          & ham%isdistributed, evec%isdistributed
        call terminate

      end if
#else
      call hesolver(ham%za, evec%za, eval, i1, i2, v1, v2, found)
#endif

      contains

#ifdef SCAL
        subroutine workspacequery(rangetype)

          character(1), intent(in) :: rangetype

          integer(4) :: nhetrd_lwork, anb, nps, n
          integer(4), external :: pjlaenv, iceil
          integer(4) :: sqnpc

          ! If automatic lwork is smaller than nhetrd_lwork, set it to this value
          n = ham%nrows
          anb = pjlaenv(binfo%context, 3, 'PZHETTRD', 'L', 0, 0, 0, 0)
          sqnpc = int(sqrt(dble(binfo%nprows*binfo%npcols)))
          nps = max(numroc(n, 1, 0, 0, sqnpc), 2*anb)
          nhetrd_lwork = n + 2*(anb+1)*(4*nps+2)+(nps+1)*nps

          ! lrwork
          ! to guarantee orthonormal vectors add (clustersize-1)*n to lrwork

          ! Automatic calculation of work array lengths
          lwork=-1
          lrwork=-1
          liwork=-1

          allocate(work(3), rwork(3), iwork(3))

          call pzheevx('V', rangetype, 'U', ham%nrows, ham%za, ia, ja, ham%desc,&
            & vl, vu, il, iu, abstol, nevalfound, nevecfound, eval, orfac,&
            & evec%za, iz, jz, evec%desc, work, lwork, rwork, lrwork, iwork, liwork,&
            & ifail, iclustr, gap, info)

          ! Adjust workspace
          lwork=max(int(work(1)), nhetrd_lwork)
          lrwork=int(rwork(1)) + (clusterguess-1)*n
          liwork=int(iwork(1))

          deallocate(work, rwork, iwork)
        end subroutine workspacequery

        subroutine errorinspect(ierror)
          integer(4), intent(in) :: ierror

          integer(4) :: i, maxcs, tmp

          if( ierror < 0) then
            write(*,'("Error (dhesolver) cause: Invalid input")')

          else if(mod(ierror,2) /= 0) then
            write(*,'("Error (dhesolver) cause: Eigenvectors not converged")')
            write(*,'("dhesolver ifail")')
            write(*,'(I8)') ifail

          else if(mod(ierror/2,2) /= 0) then
            ! Inspect iclstr
            do i = 1, size(iclustr)-1
              tmp = iclustr(i+1) - iclustr(i)
              if(tmp > 0) then 
                maxcs = max(tmp, maxcs)
              else 
                exit
              end if
            end do
            i = i-1
            write(*,'("Warning (dhesolver) cause: Reorthogonalization failed,&
              & insufficent workspace. There are", i4," clusters of eignevalues&
              & and the lagest one has size ", i4,". &
              & Increase eecs to ", i4," to garantee orthogonal eigenvectors")')&
              & i, maxcs, maxcs
            write(*,'("dhesolver iclustr:")')
            write(*,'(I8)') iclustr

          else if(mod(ierror/4,2) /= 0) then
            write(*,'("Error (dhesolver) cause:&
              & Not all eigenvectors computed, insufficent workspace.")')

          else if(mod(ierror/8,2) /= 0) then
            write(*,'("Error (dhsolver) cause:&
              & Eigenvalue comptaion failed")')
          end if
        end subroutine errorinspect
#endif

    end subroutine dhesolver
    !EOC
    
end module m_dhesolver
