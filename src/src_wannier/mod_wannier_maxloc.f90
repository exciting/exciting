module mod_wannier_maxloc
  use mod_wannier_variables
  use m_linalg

  use mod_constants,            only: zone, zzero, zi, pi, twopi
  use modinput

  implicit none

! module variables
  ! parameters
  integer :: minit, maxit
  integer :: maxcg
  integer :: maxls
  integer :: newls
  integer :: nwrite
  real(8) :: sl0
  real(8) :: noise_level
  real(8) :: etahz = 0.01d0

  ! gradients
  complex(8), allocatable :: mlwf_grad(:,:,:)
  complex(8), allocatable :: mlwf_grad_last(:,:,:)
  real(8) :: mlwf_lingrad
  real(8) :: mlwf_gradnorm

  ! update directions
  complex(8), allocatable :: mlwf_dir(:,:,:)
  complex(8), allocatable :: mlwf_dir_last(:,:,:)
  complex(8), allocatable :: mlwf_dir_evec(:,:,:)
  real(8), allocatable :: mlwf_dir_eval(:,:)

  ! counters
  integer :: cntit, cntls, cntcg, cntomega, cntgrad
  
  ! others
  real(8) :: cgparam
  real(8) :: step, last_step
  real(8) :: mlwf_omega

! methods
  contains

    !BOP
    ! !ROUTINE: wannier_gen_max
    ! !INTERFACE:
    !
    subroutine wannier_maxloc
      ! !USES:
      use m_plotmat
      use mod_wannier_filehandling
      ! !INPUT/OUTPUT PARAMETERS:
      !   bandstart : n1, lowest band index of the band range used for generation of
      !               Wannier functions (in,integer)
      !   nband     : N, number of bands used for generation (in,integer)
      !   loproj    : indices of local-orbitals used as projection functions
      !               (in, integer(nband))
      ! !DESCRIPTION:
      !   Does the same thing as {\tt genwf} does but the transformation
      !   matrices are used for the generation of maximally localized Wannier
      !   functions. The matrices are computed in a self consistent loop. 
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (SeTi)
      !EOP
      !BOC

      ! local variables
      integer :: ik, convun
      real(8) :: grad, nwgt
      real(8) :: omegastart, omegamean, uncertainty, omega_old
      real(8) :: t0, t1, t2
      complex(8) :: auxmat( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf)
      logical :: success
      character(64) :: convfname

      ! allocatable arrays
      integer, allocatable :: integerlist(:)
      real(8), allocatable :: omegahist(:), gradhist(:)

      ! external
      complex(8) :: zdotc

      minit = input%properties%wannier%grouparray( wf_group)%group%minit
      if( minit .le. 0) minit = min( 100, wf_groups( wf_group)%nwf*(wf_groups( wf_group)%nwf+1)/2)
      maxit = input%properties%wannier%grouparray( wf_group)%group%maxit
      maxcg = min( 100, wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf)
      maxcg = min( 20, wf_groups( wf_group)%nwf*(wf_groups( wf_group)%nwf+1)/2)
      maxls = 20
      newls = 1
      sl0 = input%properties%wannier%grouparray( wf_group)%group%step
      noise_level = input%properties%wannier%grouparray( wf_group)%group%noise
      nwrite = input%properties%wannier%grouparray( wf_group)%group%nwrite
      if( nwrite .le. 0) nwrite = maxit

      allocate( mlwf_grad(      wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_kset%nkpt))
      allocate( mlwf_grad_last( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_kset%nkpt))
      allocate( mlwf_dir(       wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_kset%nkpt))
      allocate( mlwf_dir_last(  wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_kset%nkpt))
      allocate( mlwf_dir_evec(  wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_kset%nkpt))
      allocate( mlwf_dir_eval(  wf_groups( wf_group)%nwf, wf_kset%nkpt))

!#ifdef MPI
!      if( rank .eq. 0) then
!#endif
      if( input%properties%wannier%grouparray( wf_group)%group%writeconv) then
        call getunit( convun)
        write( convfname, '("maxloc_conv_",i3.3,".dat")'), wf_group
        open( convun, file=trim( convfname), action='write', form='formatted')
        write( convun, '("# Iteration   Time   Omega   DeltaOmega   dOmega   Step   Gradient   Uncertainty   CGstep   OmegaEvaluations")')
      end if

      write( wf_info, '(" minimize localization functional Omega...")')
      call timesec( t0)

      !********************************************************************
      ! initialize M matrices and Wannier centers
      !********************************************************************
      call wannier_loc
        
      !********************************************************************
      ! minimize localization functional
      !********************************************************************
      step = sl0
      last_step = step

      ! start minimization loop
      allocate( omegahist( minit), gradhist( minit), integerlist( minit))
      do cntit = 1, minit
        integerlist( cntit) = cntit
      end do
      cntit = 0
      cntcg = 0
      omegastart = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      mlwf_omega = omegastart
      omega_old = mlwf_omega
      omegahist = 0.d0
      gradhist = 1.d0
      grad = 1.d0
      success = .false.
      uncertainty = 1.d0

      ! combined neighbor weights
      nwgt = 4.d0*sum( wf_n_wgt)

      !************************************************************
      !* start of minimization
      !************************************************************
      call timesec( t2)
      do while( .not. success)
        cntit = cntit + 1
        cntomega = 0
        cntgrad = 0
    
        call wannier_phases
        !write(*,'(2(3f13.6,6x),10f9.5)') wf_rguide( :, 1), wf_centers( :, 1), wf_sheet( 1, :)
        !write(*,'(12(2f8.3,5x))') wf_m( 1, 1, 1, :)
        !write(*,'(12(2f8.3,5x))') wf_transform( 1, 1:wf_n_ntot, 1)
        !write(*,*)

        call wannier_gradient( mlwf_grad)
        mlwf_gradnorm = sqrt( dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_grad, 1, mlwf_grad, 1)))

        !-------------------------------------------------------
        !- calculate CG parameter
        !-------------------------------------------------------
        if( (cntit .gt. 1)) then              
          call wannier_getcgparam( cgparam)
          cntcg = cntcg + 1
        else
          cgparam = 0.d0
          cntcg = 0
        end if

        !-------------------------------------------------------
        !- add direction
        !-------------------------------------------------------
        mlwf_dir = -mlwf_grad + cgparam*mlwf_dir_last

        mlwf_lingrad = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_grad, 1, mlwf_dir, 1))/nwgt

        ! reset if necessary
        if( (mod( cntcg, maxcg) .eq. 0) .or. (mlwf_lingrad .ge. 0.d0)) then
          mlwf_dir = -mlwf_grad
          mlwf_lingrad = -mlwf_gradnorm*mlwf_gradnorm/nwgt
          cntcg = 0
        end if
        
        mlwf_grad_last = mlwf_grad
        mlwf_dir_last = mlwf_dir

#ifdef USEOMP
!!$omp parallel default( shared) private( ik, auxmat)
!!$omp do
#endif
        ! diagonalize update directions
        do ik = 1, wf_kset%nkpt
          auxmat = zi*mlwf_dir( :, :, ik)/nwgt
          call zhediag( auxmat, mlwf_dir_eval( :, ik), mlwf_dir_evec( :, :, ik))
        end do
#ifdef USEOMP
!!$omp end do
!!$omp end parallel
#endif
        omega_old = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
        step = max( sl0, last_step)
        !step = sl0
        if( input%properties%wannier%grouparray( wf_group)%group%ls .and. (mod( cntit-1, newls) .eq. 0)) then
          ! perform line-search
          call wannier_linesearch( omega_old, mlwf_lingrad)
          !call wannier_linesearch_naiv( omega_old)
        else
          call wannier_update( step, mlwf_dir_eval, mlwf_dir_evec, new=.true.)
          cntomega = cntomega + 1
          call wannier_loc( totonly=.true.)
          mlwf_omega = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
        end if

        if( (noise_level .gt. input%properties%wannier%grouparray( wf_group)%group%epsmaxloc) .and. (abs( omega_old - mlwf_omega) .lt. 1.d-3)) then
          cntcg = -1
          call wannier_addnoise( noise_level)
        end if
        
        !write(*,'(3F13.6)') wf_centers( :, 12)
        
        !-------------------------------------------------------
        !- convergence analysis
        !-------------------------------------------------------
        if( cntit .eq. 1) omegahist(:) = 0.d0
        omegahist = cshift( omegahist, -1)
        gradhist = cshift( gradhist, -1)
        omegahist(1) = mlwf_omega
        omegamean = sum( omegahist(:))/min( minit, cntit)
        if( cntit .eq. 1) then
          uncertainty = 1.d0
          grad = 1.d0
        else
          uncertainty = sqrt( sum( (omegahist(:)-omegamean)**2)/(min( minit, cntit)-1))/abs( mlwf_omega)
          grad = dot_product( dble( integerlist( 1:min( minit, cntit))), omegahist( 1:min( minit, cntit))-omegamean) - (min( minit, cntit)+1)*0.5d0*sum( omegahist( 1:min( minit, cntit))-omegamean)
          grad = grad/sum( (dble( integerlist( 1:min( minit, cntit)))-(min( minit, cntit)+1)*0.5d0)**2)/abs( mlwf_omega)
          gradhist(1) = grad
        end if

        if( (mod( cntit, nwrite) .eq. 0) .and. (cntit .ne. maxit)) call wannier_writetransform

        if( mlwf_omega .gt. omegastart+abs( omegastart)) success = .true.
        if( uncertainty .le. input%properties%wannier%grouparray( wf_group)%group%uncertainty) success = .true.
        if( maxval( gradhist) .le. input%properties%wannier%grouparray( wf_group)%group%uncertainty) success = .true.
        if( cntit .lt. minit) success = .false.
        if( mlwf_gradnorm .lt. input%properties%wannier%grouparray( wf_group)%group%epsmaxloc) success = .true.
        if( cntit .ge. maxit) success = .true.
        call timesec( t1)
        if( input%properties%wannier%grouparray( wf_group)%group%writeconv) then
          write( convun, '(i4,7f23.16,2i4)') cntit, t1-t2, mlwf_omega, omega_old-mlwf_omega, mlwf_gradnorm, step, mlwf_lingrad, uncertainty, cntcg, cntomega!, maxval( gradhist)!, omega2 !omegai+omegad+omegaod
        end if
      end do
      !************************************************************
      !* end of minimization
      !************************************************************

      if( mlwf_omega .gt. omegastart+abs( omegastart)) then
        write(*,*)
        write(*, '("Error (wannier_maxloc): Localization functional diverged. Procedure aborted after ",I4," loops.")') cntit
      else if( cntit .ge. maxit) then
        write(*,*)
        write(*, '("Error (wannier_maxloc): Not converged after ",I6," cycles.")') maxit
      end if
        
      call wannier_loc
  
      deallocate( omegahist, gradhist, integerlist)
      deallocate( mlwf_grad, mlwf_grad_last, mlwf_dir, mlwf_dir_last, mlwf_dir_eval, mlwf_dir_evec)

      call timesec( t1)
      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"minimum/maximum iterations: ",T40,I6,"/",I6)') minit, maxit
      write( wf_info, '(5x,"iterations: ",T40,7x,I6)') cntit
      write( wf_info, '(5x,"gradient cutoff: ",T40,E13.6)') input%properties%wannier%grouparray( wf_group)%group%epsmaxloc
      write( wf_info, '(5x,"norm of gradient: ",T40,E13.6)') mlwf_gradnorm
      if( input%properties%wannier%grouparray( wf_group)%group%uncertainty .gt. 0.d0) then
        write( wf_info, '(5x,"aimed uncertainty: ",T40,E13.6)') input%properties%wannier%grouparray( wf_group)%group%uncertainty
        write( wf_info, '(5x,"achieved uncertainty: ",T40,E13.6)') uncertainty
      end if
      write( wf_info, '(5x,"Omega: ",T40,F13.6)') mlwf_omega
      write( wf_info, '(5x,"reduction: ",T40,7x,I5,"%")') nint( 100d0*(omegastart-mlwf_omega)/omegastart)
      write( wf_info, *)
      call flushifc( wf_info)
!#ifdef MPI
!        call barrier
!      else
!        call barrier
!      end if
!#endif
      return
      !EOC
    end subroutine wannier_maxloc
    !EOP

    subroutine wannier_update( step, direval, direvec, new)
      real(8), intent( in) :: step, direval( wf_groups( wf_group)%nwf, wf_kset%nkpt)
      complex(8), intent( in) :: direvec( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_kset%nkpt)
      logical, optional, intent( in) :: new

      integer :: ik, ist
      real(8) :: used_step, min_step, max_step
      real(8), save :: last_step = 0.d0
      complex(8) :: auxmat1( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf), auxmat2( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf)
      logical :: update_step

      min_step = 1.d-3
      max_step = 1.d3

      update_step = .true.
      if( present( new)) update_step = .not. new
      if( update_step) then
        used_step = max( min( step, max_step), min_step) - last_step
      else
        used_step = max( min( step, max_step), min_step)
      end if
      last_step = max( min( step, max_step), min_step)

#ifdef USEOMP
!$omp parallel default( shared) private( ik, ist, auxmat1, auxmat2)
!$omp do
#endif
      do ik = 1, wf_kset%nkpt
        do ist = 1, wf_groups( wf_group)%nwf
          auxmat1( :, ist) = exp( -zi*used_step*direval( ist, ik))*direvec( :, ist, ik) 
        end do
        call zgemm( 'n', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
               wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik), wf_groups( wf_group)%nst, &
               auxmat1, wf_groups( wf_group)%nwf, zzero, &
               auxmat2, wf_groups( wf_group)%nst)
        call zgemm( 'n', 'c', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
               auxmat2, wf_groups( wf_group)%nst, &
               direvec( :, :, ik), wf_groups( wf_group)%nwf, zzero, &
               wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik), wf_groups( wf_group)%nst)
      end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
      return
    end subroutine wannier_update

    subroutine wannier_gradient( grad)
      complex(8), intent( out) :: grad( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_kset%nkpt)

      integer :: ik, idxn, ikb, ist
      real(8) :: a
      complex(8) :: r( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf), &
                    t( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf)

      grad = zzero
#ifdef USEOMP
!$omp parallel default( shared) private( ik, idxn, ist, ikb, a, r, t)
!$omp do
#endif
      do ik = 1, wf_kset%nkpt
        !-------------------------------------------------------
        !- calculate gradient
        !-------------------------------------------------------
        do idxn = 1, wf_n_ntot
          ikb = wf_n_ik2( idxn, ik)
          ! calculating R and T
          r = wf_m( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik, idxn)
          t = wf_m( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik, idxn)
          do ist = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
            a = 0.d0
            if( abs( wf_m( ist, ist, ik, idxn)) .gt. 1.d-10) then
              a = aimag( log( wf_m( ist, ist, ik, idxn))) - twopi*wf_sheet( ist, idxn)
              r( :, ist) = r( :, ist)*conjg( wf_m( ist, ist, ik, idxn))
              t( :, ist) = t( :, ist)/wf_m( ist, ist, ik, idxn)*(a + dot_product( wf_n_vc( :, idxn), wf_centers( :, ist)))
            else
              r( :, ist) = zzero
              t( :, ist) = zzero
            end if
          end do
          r = r - conjg( transpose( r))
          t = t + conjg( transpose( t))
          ! calculating gradient
          grad( :, :, ik) = grad( :, :, ik) - wf_n_wgt( idxn)*( 0.5*r(:,:) + 0.5*zi*t(:,:))

          ! calculating R and T
          r = conjg( transpose( wf_m( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ikb, idxn)))
          t = conjg( transpose( wf_m( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ikb, idxn)))
          do ist = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
            a = 0.d0
            if( abs( wf_m( ist, ist, ikb, idxn)) .gt. 1.d-10) then
              a = aimag( log( wf_m( ist, ist, ikb, idxn))) - twopi*wf_sheet( ist, idxn)
              r( :, ist) = r( :, ist)*wf_m( ist, ist, ikb, idxn)
              t( :, ist) = -t( :, ist)/conjg( wf_m( ist, ist, ikb, idxn))*(a + dot_product( wf_n_vc( :, idxn), wf_centers( :, ist)))
            else
              r( :, ist) = zzero
              t( :, ist) = zzero
            end if
          end do
          r = r - conjg( transpose( r))
          t = t + conjg( transpose( t))
          ! calculating gradient
          grad( :, :, ik) = grad( :, :, ik) - wf_n_wgt( idxn)*( 0.5*r(:,:) + 0.5*zi*t(:,:))
        end do
        cntgrad = cntgrad + 1
      end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
      grad = 4.0d0*grad/wf_kset%nkpt

      return
    end subroutine wannier_gradient

    subroutine wannier_getcgparam( cgp)
      real(8), intent( out) :: cgp

      ! external
      complex(8) :: zdotc

      real(8) :: r1, r2, r3, cgpmin
      cgpmin = 0.d0
      ! Hestenes-Stiefel
      if( input%properties%wannier%grouparray( wf_group)%group%cg .eq. "hs") then
        r1 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_grad, 1, mlwf_grad - mlwf_grad_last, 1))
        r2 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_dir_last, 1, mlwf_grad - mlwf_grad_last, 1))
        cgpmin = 0.d0

      ! Fletcher-Reeves
      else if( input%properties%wannier%grouparray( wf_group)%group%cg .eq. "fr") then
        r1 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_grad, 1, mlwf_grad, 1))
        r2 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_grad_last, 1, mlwf_grad_last, 1))
        cgpmin = 0.d0

      ! Polak-Ribier
      else if( input%properties%wannier%grouparray( wf_group)%group%cg .eq. "pr") then
        r1 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_grad, 1, mlwf_grad - mlwf_grad_last, 1))
        r2 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_grad_last, 1, mlwf_grad_last, 1))
        cgpmin = 0.d0
      
      ! Hager-Zhang
      else if( input%properties%wannier%grouparray( wf_group)%group%cg .eq. "hz") then
        r1 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_grad - mlwf_grad_last, 1, mlwf_grad - mlwf_grad_last, 1))
        r2 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_dir_last, 1, mlwf_grad - mlwf_grad_last, 1))
        r3 = 2.d0*r1/r2
        r1 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_dir_last, 1, mlwf_dir_last, 1))
        r2 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_grad_last, 1, mlwf_grad_last, 1))
        cgpmin = -1.d0/(r1*min( etahz, r2))
        r1 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_grad - mlwf_grad_last - r3*mlwf_dir_last, 1, mlwf_grad, 1))
        r2 = dble( zdotc( wf_kset%nkpt*wf_groups( wf_group)%nwf*wf_groups( wf_group)%nwf, mlwf_dir_last, 1, mlwf_grad - mlwf_grad_last, 1))
      
      ! steepest descent
      else
        r1 = 0.d0
        r2 = 1.d0
      end if

      cgp = max( cgpmin, r1/r2)
      return
    end subroutine wannier_getcgparam

    subroutine wannier_linesearch_naiv( omega0)
      real(8), intent( in) :: omega0

      integer :: it, m
      real(8) :: l, o(3), x(3)

      it = 0
      l = step
      o(1) = omega0
      x(1) = 0.d0
      call wannier_update( step, mlwf_dir_eval, mlwf_dir_evec, new=.true.)
      call wannier_loc( totonly=.true.)
      o(2) = sum( wf_omega)
      x(2) = step
      cntomega = cntomega + 1
      call wannier_update( step+l, mlwf_dir_eval, mlwf_dir_evec)
      call wannier_loc( totonly=.true.)
      o(3) = sum( wf_omega)
      x(3) = step + l
      cntomega = cntomega + 1

      !write(*,*)
      !write(*,'(i4,2f13.6,6x,3f13.6,6x,3f13.6)') it, step, l, o, x
      do while( (l .gt. input%properties%wannier%grouparray( wf_group)%group%epsmaxloc) .and. (it .lt. 20))
        it = it + 1
        m = minloc( o, 1)
        if( (m .eq. 2) .or. ( (m .eq. 1) .and. (abs( step-l) .gt. 1.d-10))) then
          o(2) = o(m)
          x(2) = x(m)
          step = x(2)
          l = l/2.d0
          call wannier_update( step-l, mlwf_dir_eval, mlwf_dir_evec)
          call wannier_loc( totonly=.true.)
          o(1) = sum( wf_omega)
          x(1) = step - l
          cntomega = cntomega + 1
          call wannier_update( step+l, mlwf_dir_eval, mlwf_dir_evec)
          call wannier_loc( totonly=.true.)
          o(3) = sum( wf_omega)
          x(3) = step + l
          cntomega = cntomega + 1
        else
          if( m .eq. 3) then
            o(1) = o(3)
            x(1) = x(3)
            step = step + 2.d0*l
          else
            l = l/3.d0
            step = l
          end if
          call wannier_update( step, mlwf_dir_eval, mlwf_dir_evec)
          call wannier_loc( totonly=.true.)
          o(2) = sum( wf_omega)
          x(2) = step
          cntomega = cntomega + 1
          call wannier_update( step+l, mlwf_dir_eval, mlwf_dir_evec)
          call wannier_loc( totonly=.true.)
          o(3) = sum( wf_omega)
          x(3) = step + l
          cntomega = cntomega + 1
        end if
        !write(*,'(i4,2f13.6,6x,3f13.6,6x,3f13.6)') it, step, l, o, x
      end do
      !write(*,*)
      m = minloc( o, 1)
      step = x(m)
      call wannier_update( step, mlwf_dir_eval, mlwf_dir_evec)
      call wannier_loc( totonly=.true.)
      mlwf_omega = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      cntomega = cntomega + 1
    end subroutine wannier_linesearch_naiv

    subroutine wannier_linesearch( omega0, lingrad0)
      real(8), intent( in) :: omega0, lingrad0

      integer :: n
      real(8) :: lspar(2,3), m3(3,3), m3i(3,3), v3(3), eps1, eps2, a, b

      real(8) :: r3mdet
      
      eps1 = 0.2
      eps2 = 2.0

      lspar( 1, 1) = 0.d0
      lspar( 2, 1) = omega0

      call wannier_update( step*eps2, mlwf_dir_eval, mlwf_dir_evec, new=.true.)
      call wannier_loc( totonly=.true.)
      lspar( 1, 3) = step*eps2
      lspar( 2, 3) = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      cntomega = cntomega + 1

      call wannier_update( step, mlwf_dir_eval, mlwf_dir_evec)
      call wannier_loc( totonly=.true.)
      lspar( 1, 2) = step
      lspar( 2, 2) = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      mlwf_omega = lspar( 2, 2)
      cntomega = cntomega + 1

      ! minimum of parabolic fit
      m3( 1, :) = (/lspar( 1, 1)**2, lspar( 1, 1), 1.d0/)
      m3( 2, :) = (/lspar( 1, 2)**2, lspar( 1, 2), 1.d0/)
      m3( 3, :) = (/lspar( 1, 3)**2, lspar( 1, 3), 1.d0/)
      if( abs( r3mdet( m3)) .gt. 1.d-16) then 
        call r3minv( m3, m3i)
        call r3mv( m3i, lspar( 2, :)-minval( lspar( 2, :)), v3)
        a = -0.5d0*v3(2)/v3(1)
      else
        a = 0.d0
      end if
      !write(*,'(I,2F13.6,4(2F13.6,"  --  "))') cntomega, step, lspar( 2, 2), lspar( :, 1), lspar( :, 2), lspar( :, 3), a, -0.25d0*v3(2)**2/v3(1)+v3(3)

      ! condition 1 and 2 fulfilled
      if( (lspar( 2, 2) .le. omega0 + lspar( 1, 2)*eps1*lingrad0) .and. (lspar( 2, 3) .gt. omega0 + lspar( 1, 3)*eps1*lingrad0)) then
        ! do nothing

      ! condition 1 violated
      ! step too large
      else if( lspar( 2, 2) .gt. omega0 + lspar( 1, 2)*eps1*lingrad0) then
        b = maxval( lspar( 1, :))     !maximum step length
        mlwf_omega = lspar( 2, 2)
        step = lspar( 1, 2)

        do while( (cntomega .lt. 20) .and. (mlwf_omega .gt. omega0 + step*eps1*lingrad0))
          if( (a .le. 0.d0) .or. (a .gt. b) .or. (v3(1) .lt. 0.d0)) then
            step = step/eps2
          else
            step = a
          end if
          call wannier_update( step, mlwf_dir_eval, mlwf_dir_evec)
          call wannier_loc( totonly=.true.)
          mlwf_omega = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
          cntomega = cntomega + 1
          n = maxloc( lspar( 2, :), 1)  !position of Omega_max
          if( mlwf_omega .le. lspar( 2, n)) then
            lspar( 1, n) = step
            lspar( 2, n) = mlwf_omega
          end if
          m3( 1, :) = (/lspar( 1, 1)**2, lspar( 1, 1), 1.d0/)
          m3( 2, :) = (/lspar( 1, 2)**2, lspar( 1, 2), 1.d0/)
          m3( 3, :) = (/lspar( 1, 3)**2, lspar( 1, 3), 1.d0/)
          a = 0.d0
          if( abs( r3mdet( m3)) .gt. 1.d-16) then 
            call r3minv( m3, m3i)
            call r3mv( m3i, lspar( 2, :)-minval( lspar( 2, :)), v3)
            if( abs( v3(1)) .gt. 1.d-16) a = -0.5d0*v3(2)/v3(1)
          end if
          b = maxval( lspar( 1, :))     !maximum step length
          !write(*,'(I,2F13.6,4(2F13.6,"  --  "))') cntomega, step, mlwf_omega, lspar( :, 1), lspar( :, 2), lspar( :, 3), a, -0.25d0*v3(2)**2/v3(1)+v3(3)
          if( abs( a/step-1.d0) .lt. 1.d-3) exit
        end do
        last_step = step

      ! condition 2 violated
      ! step too short
      else if( lspar( 2, 3) .le. omega0 + lspar( 1, 3)*eps1*lingrad0) then
        b = minval( lspar( 1, :))     !minimum step length
        mlwf_omega = lspar( 2, 3)
        step = lspar( 1, 3)

        do while( (cntomega .lt. 20) .and. (mlwf_omega .le. omega0 + step*eps1*lingrad0))
          if( (a .le. b) .or. (v3(1) .lt. 0.d0)) then
            step = step*eps2
          else
            step = a
          end if
          call wannier_update( step, mlwf_dir_eval, mlwf_dir_evec)
          call wannier_loc( totonly=.true.)
          mlwf_omega = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
          cntomega = cntomega + 1
          n = maxloc( lspar( 2, :), 1)  !position of Omega_max
          if( mlwf_omega .le. lspar( 2, n)) then
            lspar( 1, n) = step
            lspar( 2, n) = mlwf_omega
          end if
          m3( 1, :) = (/lspar( 1, 1)**2, lspar( 1, 1), 1.d0/)
          m3( 2, :) = (/lspar( 1, 2)**2, lspar( 1, 2), 1.d0/)
          m3( 3, :) = (/lspar( 1, 3)**2, lspar( 1, 3), 1.d0/)
          a = 0.d0
          if( abs( r3mdet( m3)) .gt. 1.d-16) then 
            call r3minv( m3, m3i)
            call r3mv( m3i, lspar( 2, :)-minval( lspar( 2, :)), v3)
            if( abs( v3(1)) .gt. 1.d-16) a = -0.5d0*v3(2)/v3(1)
          end if
          b = minval( lspar( 1, :))     !minimum step length
          !write(*,'(I,2F13.6,4(2F13.6,"  --  "))') cntomega, step, mlwf_omega, lspar( :, 1), lspar( :, 2), lspar( :, 3), a, -0.25d0*v3(2)**2/v3(1)+v3(3)
          if( abs( a/step-1.d0) .lt. 1.d-3) exit
        end do
        last_step = step
      else
        write(*,*)
        write(*, '("Error (wannier_linesearch): Oops! This was not supposed to happen.")')
        write(*,*)
        !stop
      end if
      ! eventually check whether we have realy decreased omega. 
      ! Otherwhise go back to the smalest value tested. But change something.
      if( mlwf_omega .ge. omega0) then
        !write(*,'(3f13.6)') lspar( 1, :)
        !write(*,'(3f13.6)') lspar( 2, :)
        n = minloc( lspar( 2, :), 1)
        if( lspar( 1, n) .eq. 0.d0) then
          lspar( 2, n) = 1.d100
          n = minloc( lspar( 2, :), 1)
        end if
        step = lspar( 1, n)
        last_step = step
        call wannier_update( step, mlwf_dir_eval, mlwf_dir_evec)
        call wannier_loc( totonly=.true.)
        mlwf_omega = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
        cntomega = cntomega + 1
        if( mlwf_omega .ge. omega0) then
          step = 0d0
          last_step = step
          call wannier_update( step, mlwf_dir_eval, mlwf_dir_evec)
          call wannier_loc( totonly=.true.)
          mlwf_omega = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
          cntomega = cntomega + 1
          cntcg = 0
          !write(*,*) mlwf_omega
        end if
        !write(*,'("CAUTION: ",3F13.6)') lspar( 1, n), lspar( 2, n), mlwf_omega
      end if

    end subroutine wannier_linesearch

    subroutine wannier_linesearch_hz( val0, grad0)
      real(8), intent( in) :: val0, grad0
      
      ! parameters
      real(8) :: delta = 0.1d0
      real(8) :: sigma = 0.9d0
      real(8) :: eps0 = 1.d-6
      real(8) :: theta = 0.5d0
      real(8) :: gamma = 0.66d0

      integer :: k
      real(8) :: bracket(3), val(3), grad(3), it(3), vt(3), gt(3), eps, nwgt
      logical :: goon = .true.

      ! external
      complex(8) :: zdotc
      
      nwgt = 4.d0*sum( wf_n_wgt)
      k = 0
      eps = eps0*abs( val0)

      bracket(1) = 0.d0
      val(1) = val0
      grad(1) = grad0
      bracket(2) = 2.d0*last_step
      call wannier_update( bracket(2), mlwf_dir_eval, mlwf_dir_evec, new=.true.)
      call wannier_loc( totonly=.true.)
      val(2) = sum( wf_omega)
      cntomega = cntomega + 1
      call wannier_gradient( mlwf_grad)
      grad(2) = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad, 1, mlwf_dir, 1))/nwgt
      do while( grad(2) .lt. 0.d0)
        bracket(2) = 2.d0*bracket(2)
        call wannier_update( bracket(2), mlwf_dir_eval, mlwf_dir_evec, new=.true.)
        call wannier_loc( totonly=.true.)
        val(2) = sum( wf_omega)
        cntomega = cntomega + 1
        call wannier_gradient( mlwf_grad)
        grad(2) = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad, 1, mlwf_dir, 1))/nwgt
      end do  

      !write(*,'(i,3(2f13.6,5x))') k, bracket(1:2), val(1:2), grad(1:2)
      do while( goon)
        ! strong Wolfe condition fulfilled
        if( (val(1) - val0 .le. delta*bracket(1)*grad0) .and. (grad(1) .gt. sigma*grad0)) exit
        ! approximate Wolfe condition fulfilled
        if( ((2.d0*delta - 1.d0)*grad0 .gt. grad(1)) .and. (grad(1) .gt. sigma*grad0)) exit

        call linesearch_hz_secant2( bracket, val, grad, it(1:2), vt(1:2), gt(1:2))
        if( it(2) - it(1) .gt. gamma*( bracket(2) - bracket(1))) then
          it(3) = 0.5d0*(it(1) + it(2))
          call wannier_update( it(3), mlwf_dir_eval, mlwf_dir_evec)
          call wannier_loc( totonly=.true.)
          vt(3) = sum( wf_omega)
          cntomega = cntomega + 1
          call wannier_gradient( mlwf_grad)
          gt(3) = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad, 1, mlwf_dir, 1))/nwgt
          call linesearch_hz_update_interval( it, vt, gt, bracket(1:2), val(1:2), grad(1:2))
        else
          bracket(1:2) = it(1:2)
          val(1:2) = vt(1:2)
          grad(1:2) = gt(1:2)
        end if
        k = k + 1
        eps = eps0*abs( val(1))
        !write(*,'(i,3(2f13.6,5x))') k, bracket(1:2), val(1:2), grad(1:2)
      end do
      
      step = bracket(1)
      call wannier_update( step, mlwf_dir_eval, mlwf_dir_evec)
      call wannier_loc( totonly=.true.)
      mlwf_omega = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      cntomega = cntomega + 1
      last_step = step
      return
      contains

        subroutine linesearch_hz_update_interval( iin, vin, gin, iout, vout, gout)
          real(8), intent( in) :: iin(3), vin(3), gin(3)
          real(8), intent( out) :: iout(2), vout(2), gout(2)

          real(8) :: a, b, d, df, f

          ! external
          complex(8) :: zdotc

          if( (iin(3) .lt. iin(1)) .or. (iin(3) .gt. iout(2))) then
            iout = iin( 1:2)
            vout = vin( 1:2)
            gout = gin( 1:2)
            return
          end if
          if( gin(3) .ge. 0.d0) then
            iout(1) = iin(1)
            iout(2) = iin(3)
            vout(1) = vin(1)
            vout(2) = vin(3)
            gout(1) = gin(1)
            gout(2) = gin(3)
            return
          end if
          if( (gin(3) .lt. 0.d0) .and. (vin(3) .le. val0 + eps)) then
            iout(1) = iin(3)
            iout(2) = iin(2)
            vout(1) = vin(3)
            vout(2) = vin(2)
            gout(1) = gin(3)
            gout(2) = gin(2)
            return
          else
            a = iin(1)
            b = iin(3)
            do while( .true.)
              d = (1.d0 - theta)*a + theta*b
              call wannier_update( d, mlwf_dir_eval, mlwf_dir_evec)
              call wannier_loc( totonly=.true.)
              f = sum( wf_omega)
              cntomega = cntomega + 1
              call wannier_gradient( mlwf_grad)
              df = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad, 1, mlwf_dir, 1))/nwgt
              write(*,'(3f13.6)') d, f, df
              if( df .ge. 0.d0) then
                iout(2) = d
                vout(2) = f
                gout(2) = df
                iout(1) = a
                exit
              else
                if( f .le. val0 + eps) then
                  a = d
                  vout(1) = f
                  gout(1) = df
                else
                  b = d
                end if
              end if
            end do
          end if
          return
        end subroutine linesearch_hz_update_interval

        subroutine linesearch_hz_secant( pos, val, c, vc, gc)
          real(8), intent( in) :: pos(2), val(2)
          real(8), intent( out) :: c, vc, gc

          ! external
          complex(8) :: zdotc

          c = (pos(1)*val(2) - pos(2)*val(1))/(val(2) - val(1))
          call wannier_update( c, mlwf_dir_eval, mlwf_dir_evec)
          call wannier_loc( totonly=.true.)
          vc = sum( wf_omega)
          cntomega = cntomega + 1
          call wannier_gradient( mlwf_grad)
          gc = dble( zdotc( wf_kset%nkpt*wf_nwf*wf_nwf, mlwf_grad, 1, mlwf_dir, 1))/nwgt
          return
        end subroutine linesearch_hz_secant

        subroutine linesearch_hz_secant2( iin, vin, gin, iout, vout, gout)
          real(8), intent( in) :: iin(3), vin(3), gin(3)
          real(8), intent( out) :: iout(2), vout(2), gout(2)

          real(8) :: iAB(2), vAB(2), gAB(2), it(3), vt(3), gt(3), c, vc, gc, cp
          
          call linesearch_hz_secant( iin(1:2), gin(1:2), c, vc, gc)
          it = iin
          vt = vin
          gt = gin
          it(3) = c
          vt(3) = vc
          gt(3) = gc
          call linesearch_hz_update_interval( it, vt, gt, iAB, vAB, gAB)
          if( abs( c - iAB(2)) .lt. 1.d-10) then
            it(1) = iin(2)
            it(2) = iAB(2)
            vt(1) = vin(2)
            vt(2) = vAB(2)
            gt(1) = gin(2)
            gt(2) = gAB(2)
            call linesearch_hz_secant( it(1:2), gt(1:2), cp, vc, gc)
          end if
          if( abs( c - iAB(1)) .lt. 1.d-10) then
            it(1) = iin(1)
            it(2) = iAB(1)
            vt(1) = vin(1)
            vt(2) = vAB(1)
            gt(1) = gin(1)
            gt(2) = gAB(1)
            call linesearch_hz_secant( it(1:2), gt(1:2), cp, vc, gc)
          end if
          if( (abs( c - iAB(1)) .lt. 1.d-10) .or. (abs( c - iAB(2)) .lt. 1.d-10)) then
            it = (/iAB(1), iAB(2), cp/)
            vt = (/vAB(1), vAB(2), vc/)
            gt = (/gAB(1), gAB(2), gc/)
            call linesearch_hz_update_interval( it, vt, gt, iout, vout, gout)
          else
            iout = iAB
            vout = vAB
            gout = gAB
          end if

          return
        end subroutine linesearch_hz_secant2
    end subroutine wannier_linesearch_hz

    subroutine wannier_addnoise( nl)
      real(8), intent( inout) :: nl

      integer :: ik, ist, jst
      real(8) :: rand, eval( wf_groups( wf_group)%nwf)
      complex(8) :: z, noise( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf)
      complex(8) :: evec( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf)
      complex(8) :: auxmat1( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf), auxmat2( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf)

      nl = abs( nl)
      if( nl .gt. 1.d-10) then
        noise = zone
        call random_seed()
        do ist = 2, wf_groups( wf_group)%nwf
          do jst = 1, ist-1
            call random_number( rand)
            z = cmplx( cos( twopi*rand), sin( twopi*rand), 8)
            noise( ist, jst) = z
            noise( jst, ist) = -conjg( z)
          end do
        end do
        call zhediag( zi*noise, eval, evec)
        do ist = 1, wf_nwf
          auxmat1( :, ist) = exp( -zi*nl*eval( ist))*evec( :, ist)
        end do
        call zgemm( 'n', 'c', wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
               auxmat1, wf_groups( wf_group)%nwf, &
               evec, wf_groups( wf_group)%nwf, zzero, &
               noise, wf_groups( wf_group)%nwf)
#ifdef USEOMP
!$omp parallel default( shared) private( ik, auxmat2)
!$omp do
#endif
        do ik = 1, wf_kset%nkpt
          call zgemm( 'n', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
                 wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik), wf_groups( wf_group)%nst, &
                 noise, wf_groups( wf_group)%nwf, zzero, &
                 auxmat2, wf_groups( wf_group)%nst)
          wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik) = auxmat2
        end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
      end if
      return
    end subroutine wannier_addnoise

    subroutine wannier_loc( totonly)
      logical, optional, intent( in) :: totonly

      integer :: iknr, j, k, idxn
      complex(8), allocatable :: auxmat(:,:)
      real(8), allocatable :: logsum(:,:), log2sum(:,:), abssum(:,:), abs2sum(:,:)
      real(8) :: tmp
      logical :: tot

      tot = .false.
      if( present( totonly)) tot = totonly
      
      allocate( auxmat(  wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf))
      allocate( logsum(  wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_n_ntot))
      allocate( log2sum( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_n_ntot))
      allocate( abssum(  wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_n_ntot))
      allocate( abs2sum( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_n_ntot))

      !if( .not. wf_initialized) call wannier_init
      if( .not. allocated( wf_m0)) then
        write(*,*)
        write(*, '("Error (wannier_loc): Matrix elements not available.")')
        stop
      end if
      if( .not. allocated( wf_m)) allocate( wf_m( wf_nwf, wf_nwf, wf_kset%nkpt, wf_n_ntot))
      if( .not. allocated( wf_sheet)) then
        allocate( wf_sheet( wf_nwf, wf_n_ntot))
        wf_sheet = 0
      end if
      if( .not. allocated( wf_centers)) then
        allocate( wf_centers( 3, wf_nwf))
        wf_centers = 0.d0
      end if
      if( .not. allocated( wf_omega)) allocate( wf_omega( wf_nwf))
      if( .not. tot) then
        if( .not. allocated( wf_omega_i)) allocate( wf_omega_i( wf_nwf))
        if( .not. allocated( wf_omega_d)) allocate( wf_omega_d( wf_nwf))
        if( .not. allocated( wf_omega_od)) allocate( wf_omega_od( wf_nwf))
      end if

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, idxn, auxmat)
!$OMP DO COLLAPSE( 2)
#endif
      do iknr = 1, wf_kset%nkpt
        do idxn = 1, wf_n_ntot 
          call zgemm( 'n', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nst, zone, &
                 wf_m0( wf_groups( wf_group)%fst, wf_groups( wf_group)%fst, iknr, idxn), wf_nst, &
                 wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, wf_n_ik( idxn, iknr)), wf_nst, zzero, &
                 auxmat, wf_groups( wf_group)%nst)
          call zgemm( 'c', 'n', wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nst, zone, &
                 wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, iknr), wf_nst, &
                 auxmat, wf_groups( wf_group)%nst, zzero, &
                 wf_m( wf_groups( wf_group)%fwf, wf_groups( wf_group)%fwf, iknr, idxn), wf_nwf)
        end do
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

      logsum = 0.d0
      log2sum = 0.d0
      abssum = 0.d0
      abs2sum = 0.d0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( j, idxn, iknr, tmp, k)
!$OMP DO COLLAPSE(2)
#endif
      do idxn = 1, wf_n_ntot
        do j = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
          do iknr = 1, wf_kset%nkpt
            tmp = 0.d0
            if( abs( wf_m( j, j, iknr, idxn)) .gt. 1.d-10) then
              tmp = aimag( log( wf_m( j, j, iknr, idxn))) - twopi*wf_sheet( j, idxn)
              logsum( j, idxn) = logsum( j, idxn) + tmp
              log2sum( j, idxn) = log2sum( j, idxn) + tmp*tmp
              abssum( j, idxn) = abssum( j, idxn) + dble( wf_m( j, j, iknr, idxn)*conjg( wf_m( j, j, iknr, idxn)))
            end if
            if( .not. tot) then
              do k = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
                abs2sum( j, idxn) = abs2sum( j, idxn) + 0.5d0*( dble( wf_m( j, k, iknr, idxn)*conjg( wf_m( j, k, iknr, idxn))) + &
                                                                dble( wf_m( k, j, iknr, idxn)*conjg( wf_m( k, j, iknr, idxn))))
              end do
            end if
          end do
        end do
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      logsum = logsum/wf_kset%nkpt
      log2sum = log2sum/wf_kset%nkpt
      abssum = abssum/wf_kset%nkpt
      abs2sum = abs2sum/wf_kset%nkpt

      wf_centers( :, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf) = 0.d0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( j, idxn)
!$OMP DO 
#endif
      do j = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
        do idxn = 1, wf_n_ntot 
          wf_centers( :, j) = wf_centers( :, j) - 2.d0*wf_n_wgt( idxn)*logsum( j, idxn)*wf_n_vc( :, idxn)
        end do
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      !if( .not. tot) then
      !  do j = 1, wf_nwf
      !    !call r3mv( ainv, wf_centers( :, j), v3)
      !    write(*,'(i,3f13.5)') j, wf_centers( :, j)
      !  end do
      !end if

      wf_omega(    wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf) = 0.d0
      wf_omega_i(  wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf) = 0.d0
      wf_omega_d(  wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf) = 0.d0
      wf_omega_od( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf) = 0.d0

      if( tot) then
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( idxn, j)
!$OMP DO
#endif
        do j = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
          do idxn = 1, wf_n_ntot 
            ! total spread
            wf_omega( j) = wf_omega( j) + 2.d0*wf_n_wgt( idxn)*( 1.d0 - abssum( j, idxn) + log2sum( j, idxn))
          end do
          wf_omega( j) = wf_omega( j) - dot_product( wf_centers( :, j), wf_centers( :, j))
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      else
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( idxn, j)
!$OMP DO
#endif
        do j = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
          do idxn = 1, wf_n_ntot 
            ! total spread
            wf_omega( j) = wf_omega( j) + 2.d0*wf_n_wgt( idxn)*( 1.d0 - abssum( j, idxn) + log2sum( j, idxn))
            ! gauge independent spread
            wf_omega_i( j) = wf_omega_i( j) + 2.d0*wf_n_wgt( idxn)*( 1.d0 - abs2sum( j, idxn))
            ! diagonal spread
            wf_omega_d( j) = wf_omega_d( j) + 2.d0*wf_n_wgt( idxn)*log2sum( j, idxn)
            ! off-diagonal spread
            wf_omega_od( j) = wf_omega_od( j) + 2.d0*wf_n_wgt( idxn)*( abs2sum( j, idxn) - abssum( j, idxn))
          end do
          wf_omega( j) = wf_omega( j) - dot_product( wf_centers( :, j), wf_centers( :, j))
          wf_omega_d( j) = wf_omega_d( j) - dot_product( wf_centers( :, j), wf_centers( :, j))
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      end if

      deallocate( auxmat, logsum, log2sum, abssum, abs2sum)
      return
    end subroutine wannier_loc

    subroutine wannier_phases
      integer :: ik, ist, idxn, i

      integer :: ipiv1( wf_n_ntot)
      real(8) :: m1( wf_n_ntot, wf_n_ntot)
      real(8) :: m2( wf_n_ntot, wf_n_ntot)
      real(8) :: rhs( wf_n_ntot), l(wf_n_ntot), c(3)

      if( .not. allocated( wf_sheet)) allocate( wf_sheet( wf_nwf, wf_n_ntot))

      ! set up coefficient matrix
      do idxn = 1, wf_n_ntot
        do i = 1, wf_n_ntot
          m1( idxn, i) = wf_n_wgt( idxn)*dot_product( wf_n_vc( :, idxn), wf_n_vc( :, i))
          if( idxn .eq. i) then
            m1( idxn, i) = m1( idxn, i) - 1.d0
          end if
          m1( idxn, i) = m1( idxn, i)*wf_n_wgt( i)
        end do
      end do
      m2 = matmul( m1, transpose( m1))
      ! make a copy
      !mt1 = m1
      !! check for singularity (LU decomposition with pivoting)
      !call dgetrf2( wf_n_ntot, wf_n_ntot, m1, wf_n_ntot, ipiv1, i)
      !! set up permutation matrix
      !p = 0.d0
      !do i = 1, wf_n_ntot
      !  p(i,i) = 1.d0
      !end do
      !do i = 1, wf_n_ntot
      !  mt2 = 0.d0
      !  do n = 1, wf_n_ntot
      !    mt2(n,n) = 1.d0
      !  end do
      !  mt2( i, i) = 0.d0
      !  mt2( ipiv1(i), ipiv1(i)) = 0.d0
      !  mt2( ipiv1( i), i) = 1.d0
      !  mt2( i, ipiv1( i)) = 1.d0
      !  p = matmul( p, mt2)
      !end do
      !! find rank of coefficient matrix and linear independent columns and rows
      !n = 0
      !ipiv1 = 1
      !do i = 1, wf_n_ntot
      !  if( abs( m1(i,i)) .gt. 1.d-12) then
      !    n = n + 1
      !    ipiv1(n) = i
      !  end if
      !end do
      !! rotate coefficient matrix
      !m1 = matmul( transpose( p), matmul( mt1, p))
      !! and select linear independent columns / rows
      !m1 = m1( ipiv1, ipiv1)
      !! LU decomposition
      call dgetrf( wf_n_ntot, wf_n_ntot, m2, wf_n_ntot, ipiv1, i)

      do ist = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
        c = 0.d0
        do idxn = 1, wf_n_ntot
          l( idxn) = 0.d0
          do ik = 1, wf_kset%nkpt
            l( idxn) = l( idxn) + aimag( log( wf_m( ist, ist, ik, idxn)))
          end do
          l( idxn) = l( idxn)/wf_kset%nkpt
          c = c - 2*wf_n_wgt( idxn)*l( idxn)*wf_n_vc( :, idxn)
        end do
        do idxn = 1, wf_n_ntot
          rhs( idxn) = l( idxn) + dot_product( c, wf_n_vc( :, idxn))
        end do

        rhs = -matmul( m1, rhs)
        !rhs = rhs( ipiv1)
        call dgetrs( 'n', wf_n_ntot, 1, m2, wf_n_ntot, ipiv1, rhs, wf_n_ntot, i)
        wf_sheet( ist, :) = 0
        do i = 1, wf_n_ntot
          wf_sheet( ist, i) = nint( rhs( i)/twopi)
        end do
      end do
      
      !wf_sheet = 0
      return
    end subroutine wannier_phases

end module mod_wannier_maxloc
