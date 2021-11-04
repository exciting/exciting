module mod_wannier_opf
  use mod_wannier_variables
  use mod_wannier_helper
  use mod_wannier_omega
  !use m_linalg
  use xlapack, only: svd_divide_conquer

  implicit none

  private

! module variables
  real(8), parameter :: epssvd = 1.d-6

  integer :: N0, N  ! external and internal number of states
  integer :: J      ! number of WFs
  integer :: P      ! number of projectors
  integer :: dxo(2) ! dimension of OPF matrix
  logical :: sub    ! start from subspace
  integer, allocatable    :: dx(:,:)
  complex(8), allocatable :: OPF(:,:,:) ! OPF matrix
  complex(8), allocatable :: U0(:,:,:)  ! original transformation matrices 
  complex(8), allocatable :: A(:,:,:)   ! overlap matrices

! methods
  public :: wfopf_gen
    
  contains
    subroutine wfopf_gen( subspace)
      use m_getunit
      use mod_manopt, only: manopt_stiefel_cg, manopt_stiefel_lbfgs

      logical, optional, intent( in) :: subspace
      
      ! local variables
      integer :: convun, minit, maxit, memlen
      real(8) :: gradnorm, minstep

      integer :: ik, i
      real(8) :: t0, t1
      character( 64) :: convfname

      ! allocatable arrays
      real(8), allocatable :: sval(:)
      complex(8), allocatable :: auxmat(:,:), lvec(:,:), rvec(:,:)

      if( mpiglobal%rank .eq. 0) write( wf_info, '(" calculate improved optimized projection functions (iOPF)...")')
      call timesec( t0)
      minit    = input%properties%wannier%grouparray( wf_group)%group%minitopf
      maxit    = input%properties%wannier%grouparray( wf_group)%group%maxitopf
      gradnorm = input%properties%wannier%grouparray( wf_group)%group%epsopf
      minstep  = input%properties%wannier%grouparray( wf_group)%group%minstepopf
      memlen   = input%properties%wannier%grouparray( wf_group)%group%memlenopf

      !****************************
      !* PREPARATION
      !****************************
      sub = .false.
      if( present( subspace)) sub = subspace

      J = wf_groups( wf_group)%nwf
      P = wf_groups( wf_group)%nproj
      N0 = wf_groups( wf_group)%nst
      N = N0
      if( sub) N = J

      dxo = (/P,J/)
      allocate( dx(2,1))
      dx(:,1) = dxo

      allocate( U0(N0,J,wf_kset%nkpt))
      allocate( A(N,P,wf_kset%nkpt))
      allocate( OPF(P,J,1))

      call wfopf_init

      !****************************
      !* MINIMIZATION
      !****************************
      convun = 0
      if( input%properties%wannier%grouparray( wf_group)%group%writeconv) then
        call getunit( convun)
        if( sub) then
          write( convfname, '("opf_sub_conv_",i3.3,".dat")') wf_group
        else
          write( convfname, '("opf_conv_",i3.3,".dat")') wf_group
        end if
        open( convun, file=trim( convfname), action='write', form='formatted')
      end if

      if( input%properties%wannier%grouparray( wf_group)%group%optim .eq. 'cg') then
        call manopt_stiefel_cg( OPF, dxo, 1, dx, &
               cost=wfopf_omega, &
               grad=wfopf_gradient, &
               epsgrad=gradnorm, minit=minit, maxit=maxit, stdout=convun, minstep=minstep)
      else
        call manopt_stiefel_lbfgs( OPF, dxo, 1, dx, &
               cost=wfopf_omega, &
               !grad=wfopf_gradient, &
               costgrad=wfopf_omegagradient, &
               epsgrad=gradnorm, minit=minit, maxit=maxit, stdout=convun, minstep=minstep, memlen=memlen)
      end if

      if( input%properties%wannier%grouparray( wf_group)%group%writeconv) close( convun)

      if( allocated( wf_opf)) deallocate( wf_opf)
      allocate( wf_opf, source=OPF(:,:,1))
      call wfopf_write_opf

      !****************************
      ! find transformation matrices from OPFs
      !****************************
#ifdef USEOMP                
!$omp parallel default( shared) private( ik, auxmat, sval, lvec, rvec)
#endif
      allocate( auxmat(N,J), sval(J), lvec(N,N), rvec(J,J))
#ifdef USEOMP
!$omp do
#endif
      do ik = 1, wf_kset%nkpt
        call zgemm( 'n', 'n', N, J, P, zone, A(1,1,ik), N, wf_opf, P, zzero, auxmat, N)
        call svd_divide_conquer( auxmat, sval, lvec, rvec)
        call zgemm( 'n', 'n', N, J, J, zone, lvec, N, rvec, J, zzero, auxmat, N)
        if( sub) then
          call zgemm( 'n', 'n', N0, J, N, zone, U0(1,1,ik), N0, auxmat, N, zzero, &
               wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, ik), wf_nst)
        else
          wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, &
                        wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik) = auxmat
        end if
      end do
#ifdef USEOMP
!$omp end do
#endif
      deallocate( auxmat, sval, lvec, rvec)
#ifdef USEOMP
!$omp end parallel
#endif

      !****************************
      !* FINALIZATION
      !****************************
      ! determine phases for logarithms
      do i = 1, 5
        call wfomega_m
        call wfomega_diagphases( wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, 1), wf_nst, wf_nwf, wf_groups( wf_group)%nst)
      end do

      ! calculate spread
      call wfomega_gen

      call timesec( t1)
      if( mpiglobal%rank .eq. 0) then
        write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
        write( wf_info, '(5x,"iterations: ",T40,7x,I6)') maxit
        write( wf_info, '(5x,"gradient cutoff: ",T40,E13.6)') input%properties%wannier%grouparray( wf_group)%group%epsopf
        write( wf_info, '(5x,"norm of gradient: ",T40,E13.6)') gradnorm
        write( wf_info, '(5x,"Omega: ",T40,F13.6)') sum( wf_omega ( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
        write( wf_info, *)
        call flushifc( wf_info)
      end if

      call wfopf_destroy
      return
      !EOC
    end subroutine wfopf_gen
    !EOP

    subroutine wfopf_init
      integer :: ik

      real(8), allocatable :: eval(:)
      complex(8), allocatable :: AA(:,:), evec(:,:)

      ! copy original transformation matrices
      U0 = wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, &
                         wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, :)

      ! build projection matrices
      do ik = 1, wf_kset%nkpt
        if( sub) then
          call zgemm( 'c', 'n', J, P, N0, zone, &
                 U0(1,1,ik), N0, &
                 wf_groups( wf_group)%projection(:,:,ik), N0, zzero, &
                 A(1,1,ik), N)
        else
          A(:,:,ik) = wf_groups( wf_group)%projection(:,:,ik)
        end if
      end do

      ! initialize OPF matrix
      allocate( AA(P,P), eval(P), evec(P,P))
      AA = zzero
      do ik = 1, wf_kset%nkpt
        call zgemm( 'c', 'n', P, P, N, zone, A(1,1,ik), N, A(1,1,ik), N, zone, AA, P)
      end do
      AA = AA/dble( wf_kset%nkpt)
      call zhediag( AA, eval, evec=evec)
      do ik = 1, J
        OPF(:,ik,1) = evec(:,P-ik+1)
      end do
      deallocate( AA, eval, evec)

    end subroutine wfopf_init

    subroutine wfopf_omega( X, dxo, kx, dx, omega)
      integer, intent( in)    :: dxo(2), kx, dx(2,kx)
      complex(8), intent( in) :: X(dxo(1),dxo(2),*)
      real(8), intent( out)   :: omega

      integer :: ik, K

      real(8), allocatable :: s(:)
      complex(8), allocatable :: U(:,:), V(:,:), W(:,:)

#ifdef USEOMP
!$omp parallel default( shared) private( ik, s, V, W, U, K)
#endif
      allocate( s(J), V(N,J), W(J,J), U(N,J))
#ifdef USEOMP
!$omp do
#endif

      do ik = 1, wf_kset%nkpt
        call wfopf_X2U( X, dxo, dx(:,1), &
               A(:,:,ik), [ N, P], &
               U, [ N, J], &
               s, V, W, K)
        if( sub) then
          call zgemm( 'n', 'n', N0, J, J, zone, &
                 U0(1,1,ik), N0, &
                 U, N, zzero, &
                 wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, ik), wf_nst)
        else
          wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, &
                        wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik) = U
        end if
      end do
#ifdef USEOMP
!$omp end do
#endif
      deallocate( s, U, V, W)
#ifdef USEOMP
!$omp end parallel
#endif
      call wfomega_gen( totonly=.true.)
      omega = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      return
    end subroutine wfopf_omega

    subroutine wfopf_gradient( X, dxo, kx, dx, GX, dgo)
      integer, intent( in)     :: dxo(2), kx, dx(2,kx), dgo(2)
      complex(8), intent( in)  :: X(dxo(1),dxo(2),*)
      complex(8), intent( out) :: GX(dgo(1),dgo(2),*)

      integer :: ik, i, K
      real(8), allocatable :: s(:), F(:,:)
      complex(8), allocatable :: G(:,:), U(:,:), V(:,:), W(:,:), GU0(:,:), GU(:,:), AV(:,:), VGUW(:,:)
      complex(8), allocatable :: CNN(:,:), CJJ1(:,:), CJJ2(:,:), CNJ1(:,:), CNJ2(:,:)

      allocate( G(P,J))
      G = zzero
#ifdef USEOMP
!$omp parallel default( shared) private( ik, i, K, s, U, F, V, W, GU0, GU, AV, VGUW, CNN, CJJ1, CJJ2, CNJ1, CNJ2) reduction(+:G)
#endif
      allocate( s(J), U(N,J), F(J,J), V(N,J), W(J,J), GU0(N0,J), GU(N,J), AV(P,J), VGUW(J,J))
      allocate( CNN(N,N), CJJ1(J,J), CJJ2(J,J), CNJ1(N,J), CNJ2(N,J))
#ifdef USEOMP
!$omp do
#endif
      do ik = 1, wf_kset%nkpt
        call wfopf_X2U( X(:,:,1), dxo, dx(:,1), &
               A(:,:,ik), [ N, P ], &
               U, [ N, J ], &
               s, V, W, K)
        if( sub) then
          call zgemm( 'n', 'n', N0, J, J, zone, &
                 U0(1,1,ik), N0, &
                 U, N, zzero, &
                 wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, ik), wf_nst)
        else
          wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, &
                        wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik) = U
        end if

        call wfomega_gradu( ik, GU0, N0)
        if( sub) then
          call zgemm( 'c', 'n', N, J, N0, zone, U0(1,1,ik), N0, GU0, N0, zzero, GU, N)
        else
          GU = GU0
        end if
        call wfopf_getF( s, J, K, F)

        call zgemm( 'c', 'n', J, J, N, zone, V, N, GU, N, zzero, CJJ1, J)
        call zgemm( 'n', 'n', J, J, J, zone, CJJ1, J, W, J, zzero, VGUW, J)
        call zgemm( 'c', 'n', P, J, N, zone, A(1,1,ik), N, V, N, zzero, AV, P)

        CJJ1 = cmplx( F, 0.d0, 8)*VGUW
        CJJ1 = CJJ1 - conjg( transpose( CJJ1))
        call zgemm( 'n', 'c', J, J, J, zone, CJJ1, J, W, J, zzero, CJJ2, J)
        call zgemm( 'n', 'n', P, J, J, zone, AV, P, CJJ2, J, zone, G, P)

        if( J .lt. N) then
          CJJ1 = zzero
          do i = 1, K
            CJJ1(:,i) = W(:,i)/cmplx( s(i), 0.d0, 8)
          end do
          call zgemm( 'n', 'c', J, J, J, zone, CJJ1, J, W, J, zzero, CJJ2, J)
          call zgemm( 'n', 'n', N, J, J, zone, GU, N, CJJ2, J, zzero, CNJ1, N)
          call zgemm( 'n', 'c', N, N, J, -zone, V, N, V, N, zzero, CNN, N)
          do i = 1, N
            CNN(i,i) = CNN(i,i) + zone
          end do
          call zgemm( 'n', 'n', N, J, N, zone, CNN, N, CNJ1, N, zzero, CNJ2, N)
          call zgemm( 'c', 'n', P, J, N, zone, A(1,1,ik), N, CNJ2, N, zone, G, P)
        end if

      end do
#ifdef USEOMP
!$omp end do
#endif
      deallocate( s, U, F, V, W, GU0, GU, AV, VGUW)
      deallocate( CNN, CJJ1, CJJ2, CNJ1, CNJ2)
#ifdef USEOMP
!$omp end parallel
#endif
      GX(:,:,1) = zzero
      GX(1:P,1:J,1) = G
      deallocate( G)

      return
    end subroutine wfopf_gradient

    subroutine wfopf_omegagradient( X, dxo, kx, dx, omega, GX, dgo, funonly)
      integer, intent( in)       :: dxo(2), kx, dx(2,kx), dgo(2)
      complex(8), intent( in)    :: X(dxo(1),dxo(2),*)
      real(8), intent( out)      :: omega
      complex(8), intent( inout) :: GX(dgo(1),dgo(2),*)
      logical, optional, intent( in) :: funonly

      integer :: ik, i, K
      logical :: dograd
      real(8), allocatable :: s(:), F(:,:)
      complex(8), allocatable :: G(:,:), U(:,:), V(:,:), W(:,:), GU0(:,:), GU(:,:), AV(:,:), VGUW(:,:)
      complex(8), allocatable :: CNN(:,:), CJJ1(:,:), CJJ2(:,:), CNJ1(:,:), CNJ2(:,:)

      complex(8), external :: zdotc

      call wfomega_m

      dograd = .true.
      if( present( funonly)) dograd = .not. funonly

#ifdef USEOMP
!$omp parallel default( shared) private( ik, K, s, U, V, W)
#endif
      allocate( s(J), U(N,J), V(N,J), W(J,J))
#ifdef USEOMP
!$omp do
#endif
      do ik = 1, wf_kset%nkpt
        call wfopf_X2U( X, dxo, dx(:,1), &
               A(:,:,ik), [ N, P], &
               U, [ N, J], &
               s, V, W, K)
        if( sub) then
          call zgemm( 'n', 'n', N0, J, J, zone, &
                 U0(1,1,ik), N0, &
                 U, N, zzero, &
                 wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, ik), wf_nst)
        else
          wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, &
                        wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik) = U
        end if
      end do
#ifdef USEOMP
!$omp end do
#endif
      deallocate( s, V, W, U)
#ifdef USEOMP
!$omp end parallel
#endif

      call wfomega_gen( totonly=.true.)
      omega = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))

      if( dograd) then
        allocate( G(P,J))
        G = zzero
#ifdef USEOMP
!$omp parallel default( shared) private( ik, i, K, s, U, F, V, W, GU0, GU, AV, VGUW, CNN, CJJ1, CJJ2, CNJ1, CNJ2) reduction(+:G)
#endif
        allocate( s(J), U(N,J), V(N,J), W(J,J))
        allocate( F(J,J), GU0(N0,J), GU(N,J), AV(P,J), VGUW(J,J))
        allocate( CNN(N,N), CJJ1(J,J), CJJ2(J,J), CNJ1(N,J), CNJ2(N,J))
#ifdef USEOMP
!$omp do
#endif
        do ik = 1, wf_kset%nkpt
          call wfopf_X2U( X, dxo, dx(:,1), &
                 A(:,:,ik), [ N, P], &
                 U, [ N, J], &
                 s, V, W, K)
          call wfomega_gradu( ik, GU0, N0)
          if( sub) then
            call zgemm( 'c', 'n', N, J, N0, zone, U0(1,1,ik), N0, GU0, N0, zzero, GU, N)
          else
            GU = GU0
          end if
          call wfopf_getF( s, J, K, F)

          call zgemm( 'c', 'n', J, J, N, zone, V, N, GU, N, zzero, CJJ1, J)
          call zgemm( 'n', 'n', J, J, J, zone, CJJ1, J, W, J, zzero, VGUW, J)
          call zgemm( 'c', 'n', P, J, N, zone, A(1,1,ik), N, V, N, zzero, AV, P)

          CJJ1 = cmplx( F, 0.d0, 8)*VGUW
          CJJ1 = CJJ1 - conjg( transpose( CJJ1))
          call zgemm( 'n', 'c', J, J, J, zone, CJJ1, J, W, J, zzero, CJJ2, J)
          call zgemm( 'n', 'n', P, J, J, zone, AV, P, CJJ2, J, zone, G, P)

          if( J .lt. N) then
            CJJ1 = zzero
            do i = 1, K
              CJJ1(:,i) = W(:,i)/cmplx( s(i), 0.d0, 8)
            end do
            call zgemm( 'n', 'c', J, J, J, zone, CJJ1, J, W, J, zzero, CJJ2, J)
            call zgemm( 'n', 'n', N, J, J, zone, GU, N, CJJ2, J, zzero, CNJ1, N)
            call zgemm( 'n', 'c', N, N, J, -zone, V, N, V, N, zzero, CNN, N)
            do i = 1, N
              CNN(i,i) = CNN(i,i) + zone
            end do
            call zgemm( 'n', 'n', N, J, N, zone, CNN, N, CNJ1, N, zzero, CNJ2, N)
            call zgemm( 'c', 'n', P, J, N, zone, A(1,1,ik), N, CNJ2, N, zone, G, P)
          end if

        end do
#ifdef USEOMP
!$omp end do
#endif
        deallocate( s, V, W, U)
        deallocate( F, GU, AV, VGUW)
        deallocate( CNN, CJJ1, CJJ2, CNJ1, CNJ2)
#ifdef USEOMP
!$omp end parallel
#endif
        GX(1:P,1:J,1) = G
        deallocate( G)
      end if

      return
    end subroutine wfopf_omegagradient

    subroutine wfopf_X2U( X, dxo, dx, A, da, U, du, s, V, W, K)
      integer, intent( in)     :: dxo(2), dx(2), da(2), du(2)
      complex(8), intent( in)  :: X(dxo(1),dxo(2))
      complex(8), intent( in)  :: A(da(1),da(2))
      complex(8), intent( out) :: U(du(1),du(2))
      real(8), intent( out)    :: s(du(2))
      complex(8), intent( out) :: V(du(1),du(2)), W(du(2),du(2))
      integer, intent( out)    :: K

      real(8) :: maxs
      complex(8), allocatable :: AX(:,:), tmp(:,:)

      allocate( AX( da(1), du(2)))
      allocate( tmp( da(1), da(1)))
      call zgemm( 'n', 'n', da(1), du(2), dx(1), zone, A, da(1), X, dxo(1), zzero, AX, da(1))
      call svd_divide_conquer( AX, s, tmp, W)
      maxs = maxval( s)
      V = tmp( 1:N, 1:J)
      W = conjg( transpose( W))
      do K = J, 1, -1
        if( s(K) .gt. epssvd*maxs) exit
        s(K) = 0.d0
        !V(:,K) = zzero
        !W(:,K) = zzero
      end do
      call zgemm( 'n', 'c', du(1), du(2), J, zone, V, du(1), W, du(2), zzero, U, du(1))
      deallocate( AX, tmp)
      return
    end subroutine wfopf_X2U

    subroutine wfopf_getF( s, J, K, F)
      integer, intent( in)  :: J, K
      real(8), intent( in)  :: s(J)
      real(8), intent( out) :: F(J,J)

      integer :: m, n
      real(8) :: maxs
       
      F = 0.d0
      maxs = maxval( s)
      do m = 1, K
        do n = m, K
          if( (abs( s(m) - s(n)) .le. epssvd*maxs) .and. &
              (s(m) + s(n) .gt. epssvd*maxs)) then
            F(m,n) = 1.d0/(s(m) + s(n))
          else
            F(m,n) = (s(n) - s(m))/(s(n)*s(n) - s(m)*s(m))
          end if
          F(n,m) = F(m,n)
        end do
      end do
      return
    end subroutine wfopf_getF

    subroutine wfopf_write_opf
      use m_getunit

      integer :: i, j, n, un
      logical :: exist

      if( .not. input%properties%wannier%printproj .or. (mpiglobal%rank /= 0)) return
      inquire( file=trim( wf_filename)//"_PROJECTION"//trim(filext), exist=exist)
      if( .not. exist) return
      call getunit( un)
      open( un, file=trim( wf_filename)//"_PROJECTION"//trim(filext), status='old', position='append', action='write')

      if( (wf_group .eq. 1) .and. .not. sub) then
        n = 0
        do i = 1, wf_ngroups
          n = max( n, wf_groups( i)%nproj)
        end do
        call printbox( un, '*', "Optimized projection functions")
        write( un, *)
        write( un, '(80("-"))')
        write( un, '(7x,"#",4x)', advance='no')
        do i = 1, n
          write( un, '(i3,2x)', advance='no') i
        end do
        write( un, *)
      end if
      if( sub) then
        write( un, '(" from subspace")')
      else
        write( un, '(80("-"))')
      end if
      do i = 1, wf_groups( wf_group)%nwf
        write( un, '(4x,i4,4x)', advance='no') i
        do j = 1, wf_groups( wf_group)%nproj
          write( un, '(i3,2x)', advance='no') nint( 1.d2*abs( wf_opf(j,i))**2)
        end do
        write( un, *)
      end do
      if( (wf_group .eq. wf_ngroups) .and. (sub .or. (wf_groups( wf_group)%method == 'opf' .or. wf_groups( wf_group)%method == 'opfmax'))) then
        write( un, '(80("-"))')
        write( un, *)
      end if

      close( un)
      return
    end subroutine wfopf_write_opf

    subroutine wfopf_destroy
      if( allocated( OPF)) deallocate( OPF)
      if( allocated( U0))  deallocate( U0)
      if( allocated( A))   deallocate( A)
      if( allocated( dx))  deallocate( dx)
    end subroutine wfopf_destroy

end module mod_wannier_opf
