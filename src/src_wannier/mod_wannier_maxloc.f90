module mod_wannier_maxloc
  use mod_wannier_variables
  use mod_wannier_omega
  use m_linalg

  implicit none

  private

! module variables
  integer                 :: DXYO(2)
  integer, allocatable    :: DXY(:,:)
  complex(8), allocatable :: XY(:,:,:)

! methods
  public :: wfmax_gen, wfmax_fromfile
  contains

    subroutine wfmax_gen
      use mod_manopt, only: manopt_stiefel_cg, manopt_stiefel_lbfgs
      use m_getunit
      use m_plotmat
      !use mod_wannier_filehandling
      
      integer :: convun, minit, maxit, memlen
      real(8) :: gradnorm, minstep

      integer :: ik, kxy
      real(8) :: t0, t1, omega, omegastart
      character(256) :: convfname

      write( wf_info, '(" minimize localization functional Omega...")')
      call timesec( t0)
      minit    = input%properties%wannier%grouparray( wf_group)%group%minitmax
      maxit    = input%properties%wannier%grouparray( wf_group)%group%maxitmax
      gradnorm = input%properties%wannier%grouparray( wf_group)%group%epsmax
      minstep  = input%properties%wannier%grouparray( wf_group)%group%minstepmax
      memlen   = input%properties%wannier%grouparray( wf_group)%group%memlenmax

      !****************************
      !* PREPARATION
      !****************************
      ! initialize M matrices
      call wfomega_m
      ! try to reduce off-diagonal part of spread by phase correction
      call wfomega_diagphases( wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, 1), wf_nst, wf_nwf, wf_groups( wf_group)%nst)
      call wfomega_m
      call wfomega_diagphases( wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, 1), wf_nst, wf_nwf, wf_groups( wf_group)%nst)
      ! reorder matrices such that frozen bands are at the end
      call wfomega_shuffle(1)
      ! get sizes of X matrices
      kxy = wf_kset%nkpt
      if( wf_groups( wf_group)%method .eq. 'disFull') then
        if( (sum( wf_groups( wf_group)%win_ni) .gt. 0) .and. &
            (sum( wf_groups( wf_group)%win_ni) .ne. wf_kset%nkpt*wf_groups( wf_group)%nwf)) kxy = 2*wf_kset%nkpt
      end if
      allocate( DXY(2,kxy))
      DXY(:,1:wf_kset%nkpt) = wf_groups( wf_group)%nwf
      ! get sizes of Y matrices
      if( wf_groups( wf_group)%method .eq. 'disFull') then
        if( kxy .ne. wf_kset%nkpt) then
          do ik = 1, wf_kset%nkpt
            DXY(1,wf_kset%nkpt+ik) = wf_groups( wf_group)%win_no(ik)
            DXY(2,wf_kset%nkpt+ik) = wf_groups( wf_group)%nwf - wf_groups( wf_group)%win_ni(ik)
          end do
        else
          if( sum( wf_groups( wf_group)%win_ni) .eq. 0 .and. .not. allocated( wf_subspace)) &
            DXY(1,1:wf_kset%nkpt) = wf_groups( wf_group)%nst
        end if
      end if
      ! get initial X and Y matrices
      DXYO(1) = wf_groups( wf_group)%nst
      DXYO(2) = wf_groups( wf_group)%nwf
      allocate( XY( DXYO(1), DXYO(2), kxy))
      call wfmax_U2XY( wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, 1), (/wf_nst,wf_nwf/), &
             XY, DXYO, kxy, DXY)
      ! get initial spread
      call wfmax_omega( XY, DXYO, kxy, DXY, omegastart)

      !****************************
      !* MINIMIZATION
      !****************************
      convun = 0
      if( input%properties%wannier%grouparray( wf_group)%group%writeconv) then
        call getunit( convun)
        write( convfname, '("maxloc_conv_",i3.3,".dat")'), wf_group
        open( convun, file=trim( convfname), action='write', form='formatted')
      end if
      if( input%properties%wannier%grouparray( wf_group)%group%optim .eq. 'cg') then
        call manopt_stiefel_cg( XY, DXYO, kxy, DXY, &
               cost=wfmax_omega, &
               grad=wfmax_gradient, &
               update=wfmax_update, &
               epsgrad=gradnorm, minit=minit, maxit=maxit, stdout=convun, minstep=minstep)
      else
        call manopt_stiefel_lbfgs( XY, DXYO, kxy, DXY, &
               cost=wfmax_omega, &
               grad=wfmax_gradient, &
               update=wfmax_update, &
               epsgrad=gradnorm, minit=minit, maxit=maxit, stdout=convun, minstep=minstep, memlen=memlen)
      end if
      if( input%properties%wannier%grouparray( wf_group)%group%writeconv) close( convun)
      call wfomega_diagphases( wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, 1), wf_nst, wf_nwf, wf_groups( wf_group)%nst)

      !****************************
      !* FINALIZATION
      !****************************
      ! get U from (X,Y)
      call wfmax_XY2U( XY, DXYO, kxy, DXY, &
             wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, 1), (/wf_nst,wf_nwf/))
      ! reorder matrices back to original order
      call wfomega_shuffle(-1)
      ! get final spread
      call wfomega_gen
      omega = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
  
      call timesec( t1)
      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      !write( wf_info, '(5x,"minimum/maximum iterations: ",T40,I6,"/",I6)') minit, maxit
      write( wf_info, '(5x,"iterations: ",T40,7x,I6)') maxit
      write( wf_info, '(5x,"gradient cutoff: ",T40,E13.6)') input%properties%wannier%grouparray( wf_group)%group%epsmax
      write( wf_info, '(5x,"norm of gradient: ",T40,E13.6)') gradnorm
      !if( input%properties%wannier%grouparray( wf_group)%group%uncertainty .gt. 0.d0) then
      !  write( wf_info, '(5x,"aimed uncertainty: ",T40,E13.6)') input%properties%wannier%grouparray( wf_group)%group%uncertainty
      !  write( wf_info, '(5x,"achieved uncertainty: ",T40,E13.6)') uncertainty
      !end if
      write( wf_info, '(5x,"Omega: ",T40,F13.6)') omega
      write( wf_info, '(5x,"reduction: ",T40,7x,I5,"%")') nint( 100d0*(omegastart-omega)/omegastart)
      write( wf_info, *)
      call flushifc( wf_info)

      deallocate( DXY, XY)
      ! delete subspace if it was used. we don't need it anymore
      if( allocated( wf_subspace)) deallocate( wf_subspace)
      return
      !EOC
    end subroutine wfmax_gen
    !EOP

    subroutine wfmax_fromfile
      if( wf_groups( wf_group)%method .eq. 'disSMV') then
        if( allocated( wf_subspace)) deallocate( wf_subspace)
        allocate( wf_subspace( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_kset%nkpt))
        wf_subspace = wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, :)
      end if
      call wfmax_gen

      return
    end subroutine wfmax_fromfile

    subroutine wfmax_omega( xy, dxyo, kxy, dxy, omega)
      integer, intent( in)    :: dxyo(2), kxy, dxy(2,kxy)
      complex(8), intent( in) :: xy(dxyo(1),dxyo(2),*)
      real(8), intent( out)   :: omega

      call wfmax_XY2U( xy, dxyo, kxy, dxy, &
             wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, 1), (/wf_nst,wf_nwf/))
      call wfomega_gen( totonly=.true.)
      omega = sum( wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      return
    end subroutine wfmax_omega

    subroutine wfmax_update( xy, dxyo, kxy, dxy, it, change)
      use mod_wannier_filehandling, only: wffile_writetransform
      integer, intent( in)       :: dxyo(2), kxy, dxy(2,kxy), it
      complex(8), intent( inout) :: xy(dxyo(1),dxyo(2),*)
      logical, intent( out)      :: change
      
      integer :: nw

      change = .false.

      nw = input%properties%wannier%grouparray( wf_group)%group%nwrite
      if( (nw > 0) .and. (mod( it, nw) == 0)) call wffile_writetransform
      if( mod( it, 1) /= -1) return
      change = .true.
      call wfmax_XY2U( xy, dxyo, kxy, dxy, &
             wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, 1), (/wf_nst,wf_nwf/))
      call wfomega_m
      call wfomega_diagphases( xy, dxyo(1), dxyo(2), dxy(1,1))
      return
    end subroutine wfmax_update

    subroutine wfmax_gradient( xy, dxyo, kxy, dxy, gxy, dgo)
      integer, intent( in)     :: dxyo(2), kxy, dxy(2,kxy), dgo(2)
      complex(8), intent( in)  :: xy(dxyo(1),dxyo(2),*)
      complex(8), intent( out) :: gxy(dgo(1),dgo(2),*)

      integer :: ikx, iky, nik, nok
      complex(8), allocatable :: gu(:,:)

      call wfmax_XY2U( xy, dxyo, kxy, dxy, &
             wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, 1), (/wf_nst,wf_nwf/))

      ! update M matrices
      call wfomega_m

      ! compute gradient
#ifdef USEOMP
!$omp parallel default( shared) private( ikx, iky, gu, nik, nok)
#endif
      allocate( gu( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf))
#ifdef USEOMP
!$omp do
#endif
      do ikx = 1, wf_kset%nkpt
        ! get the gradient w.r.t. U
        call wfomega_gradu( ikx, gu, wf_groups( wf_group)%nst)
        ! apply chain rule to get gradient w.r.t. X and Y
        if( wf_groups( wf_group)%method .ne. 'disSMV' .and. wf_groups( wf_group)%method .ne. 'disFull') then
          gxy(:,:,ikx) = zzero
          gxy( 1:dxy(1,ikx), 1:dxy(2,ikx), ikx) = gu( 1:dxy(1,ikx), 1:dxy(2,ikx))
        else
          ! used fixed subspace if available
          if( wf_groups( wf_group)%method .eq. 'disSMV') then
            call zgemm( 'c', 'n', wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nst, zone, &
                   wf_subspace(1,1,ikx), wf_groups( wf_group)%nst, &
                   gu(1,1), wf_groups( wf_group)%nst, zzero, &
                   gxy(1,1,ikx), dgo(1))
          else if( sum( wf_groups( wf_group)%win_ni) .eq. 0) then
            gxy(:,:,ikx) = zzero
            gxy( 1:dxy(1,ikx), 1:dxy(2,ikx), ikx) = gu( 1:dxy(1,ikx), 1:dxy(2,ikx))
          else
            iky = wf_kset%nkpt + ikx
            gxy(:,:,ikx) = zzero
            gxy(:,:,iky) = zzero
            nik = wf_groups( wf_group)%win_ni( ikx)
            nok = wf_groups( wf_group)%win_no( ikx)
            if( nik .le. 0) then
              gxy( 1:dxy(1,iky), 1:dxy(2,iky), iky) = gu( 1:dxy(1,iky), 1:dxy(2,iky))
            else if( nik .eq. wf_groups( wf_group)%nwf) then
              gxy( 1:dxy(1,ikx), 1:dxy(2,ikx), ikx) = gu( (nok+1):(nok+dxy(1,ikx)), 1:dxy(2,ikx))
            else
              call zgemm( 'c', 'n', dxy(2,iky), dxy(2,ikx), dxy(1,iky), zone, &
                     xy(1,1,iky), dxyo(1), &
                     gu(1,1), wf_groups( wf_group)%nst, zzero, &
                     gxy(1,1,ikx), dgo(1))
              gxy( (dxy(1,ikx)-nik+1):dxy(1,ikx), 1:dxy(2,ikx), ikx) = gu( (nok+1):(nok+nik), 1:dxy(2,ikx))
              call zgemm( 'n', 'c', dxy(1,iky), dxy(2,iky), dxy(2,ikx), zone, &
                     gu(1,1), wf_groups( wf_group)%nst, &
                     xy(1,1,ikx), dxyo(1), zzero, &
                     gxy(1,1,iky), dgo(1))
            end if
          end if
        end if
      end do
#ifdef USEOMP
!$omp end do
#endif
      deallocate( gu)
#ifdef USEOMP
!$omp end parallel
#endif

      return
    end subroutine wfmax_gradient

    subroutine wfmax_U2XY( u, du, xy, dxyo, kxy, dxy)
      integer, intent( in)     :: du(2), kxy, dxy(2,kxy), dxyo(2)
      complex(8), intent( in)  :: u(du(1),du(2),*)
      complex(8), intent( out) :: xy(dxyo(1),dxyo(2),*)

      integer :: ikx, iky, ist, nik, nok
      real(8), allocatable :: val(:)
      complex(8), allocatable :: auxmat(:,:), vec(:,:)

      if( wf_groups( wf_group)%method .ne. 'disSMV' .and. wf_groups( wf_group)%method .ne. 'disFull') then
        do ikx = 1, wf_kset%nkpt
          xy(:,:,ikx) = zzero
          xy(1:dxy(1,ikx),1:dxy(2,ikx),ikx) = u(1:dxy(1,ikx),1:dxy(2,ikx),ikx)
        end do
      else
        allocate( auxmat( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst))
        allocate( vec( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst))
        allocate( val( wf_groups( wf_group)%nst))
        do ikx = 1, wf_kset%nkpt
          ! use fixed subspace if available
          if( wf_groups( wf_group)%method .eq. 'disSMV') then
            xy(:,:,ikx) = zzero
            call zgemm( 'c', 'n', dxy(1,ikx), dxy(2,ikx), wf_groups( wf_group)%nst, zone, &
                   wf_subspace(1,1,ikx), wf_groups( wf_group)%nst, &
                   u(1,1,ikx), du(1), zzero, &
                   xy(1,1,ikx), dxyo(1))
          else if( sum( wf_groups( wf_group)%win_ni) .eq. 0) then
            xy(:,:,ikx) = zzero
            xy( 1:dxy(1,ikx), 1:dxy(2,ikx), ikx) = u( 1:dxy(1,ikx), 1:dxy(2,ikx), ikx)
          else
            iky = wf_kset%nkpt + ikx
            xy(:,:,ikx) = zzero
            xy(:,:,iky) = zzero
            nik = wf_groups( wf_group)%win_ni( ikx)
            nok = wf_groups( wf_group)%win_no( ikx)
            if( nik .le. 0) then
              xy( 1:dxy(1,iky), 1:dxy(2,iky), iky) = u( 1:dxy(1,iky), 1:dxy(2,iky), ikx)
              do ist = 1, dxy(1,ikx)
                xy( ist, ist, ikx) = zone
              end do
            else if( nik .eq. wf_groups( wf_group)%nwf) then
              xy( 1:dxy(1,ikx), 1:dxy(2,ikx), ikx) = u( (nok+1):(nok+dxy(1,ikx)), 1:dxy(2,ikx), ikx)
            else
              call zgemm( 'n', 'c', nok, nok, wf_groups( wf_group)%nwf, zone, &
                     u(1,1,ikx), du(1), &
                     u(1,1,ikx), du(1), zzero, &
                     auxmat(1,1), wf_groups( wf_group)%nst)
              call zhediag( auxmat( 1:nok, 1:nok), val(1:nok), evec=vec( 1:nok, 1:nok))
              do ist = 1, wf_groups( wf_group)%nwf - nik
                xy( 1:nok, ist, iky) = vec( 1:nok, nok-ist+1)
              end do
              xy( (dxy(1,ikx)-nik+1):dxy(1,ikx), 1:dxy(2,ikx), ikx) = u( (nok+1):(nok+nik), 1:dxy(2,ikx), ikx)
              call zgemm( 'c', 'n', wf_groups( wf_group)%nwf-nik, wf_groups( wf_group)%nwf, nok, zone, &
                     xy(1,1,iky), dxyo(1), &
                     u(1,1,ikx), du(1), zzero, &
                     xy(1,1,ikx), dxyo(1))
            end if
          end if
        end do
        deallocate( auxmat, vec, val)
      end if        
      return 
    end subroutine wfmax_U2XY

    subroutine wfmax_XY2U( xy, dxyo, kxy, dxy, u, du)
      integer, intent( in)     :: dxyo(2), kxy, dxy(2,kxy), du(2)
      complex(8), intent( in)  :: xy(dxyo(1),dxyo(2),*)
      complex(8), intent( out) :: u(du(1),du(2),*)

      integer :: ikx, iky, nik, nok

      if( wf_groups( wf_group)%method .ne. 'disSMV' .and. wf_groups( wf_group)%method .ne. 'disFull') then
        do ikx = 1, wf_kset%nkpt
          u( 1:dxy(1,ikx), 1:dxy(2,ikx), ikx) = xy( 1:dxy(1,ikx), 1:dxy(2,ikx), ikx)
        end do
      else
        do ikx = 1, wf_kset%nkpt
          if( wf_groups( wf_group)%method .eq. 'disSMV') then
            call zgemm( 'n', 'n', wf_groups( wf_group)%nst, dxy(2,ikx), dxy(1,ikx), zone, &
                   wf_subspace(1,1,ikx), wf_groups( wf_group)%nst, &
                   xy(1,1,ikx), dxyo(1), zzero, &
                   u(1,1,ikx), du(1))
          else if( sum( wf_groups( wf_group)%win_ni) .eq. 0) then
            u( 1:dxyo(1), 1:dxyo(2), ikx) = zzero
              u( 1:dxy(1,ikx), 1:dxy(2,ikx), ikx) = xy( 1:dxy(1,ikx), 1:dxy(2,ikx), ikx)
          else
            iky = wf_kset%nkpt + ikx
            u( 1:dxyo(1), 1:dxyo(2), ikx) = zzero
            nik = wf_groups( wf_group)%win_ni( ikx)
            nok = wf_groups( wf_group)%win_no( ikx)
            if( nik .le. 0) then
              u( 1:dxy(1,iky), 1:dxy(2,iky), ikx) = xy( 1:dxy(1,iky), 1:dxy(2,iky), iky)
            else if( nik .eq. wf_groups( wf_group)%nwf) then
              u( (nok+1):(nok+nik), 1:dxy(2,ikx), ikx) = xy( 1:nik, 1:dxy(2,ikx), ikx)
            else
              u( (nok+1):(nok+nik), 1:dxy(2,ikx), ikx) = xy( (dxy(1,ikx)-nik+1):dxy(1,ikx), 1:dxy(2,ikx), ikx)
              call zgemm( 'n', 'n', nok, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf-nik, zone, &
                     xy(1,1,iky), dxyo(1), &
                     xy(1,1,ikx), dxyo(1), zzero, &
                     u(1,1,ikx), du(1))
            end if
          end if
        end do
      end if
      return
    end subroutine wfmax_XY2U

end module mod_wannier_maxloc
