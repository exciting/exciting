module mod_wannier_disentangle
  use mod_wannier_variables
  use mod_wannier_opf
  use mod_wannier_helper
  use mod_wannier_omega
  !use m_linalg
  use xlapack, only: svd_divide_conquer

  implicit none

! module variables
    
contains

    subroutine wfdis_gen
      use m_getunit
      use mod_manopt, only: manopt_stiefel_lbfgs, manopt_stiefel_cg

      integer :: convun, minit, maxit, memlen
      real(8) :: gradnorm, minstep

      integer :: ik, n, nik, nok, dyo(2)
      real(8) :: omegai0, omegai, t0, t1
      character(256) :: convfname

      integer, allocatable :: dy(:,:), map(:,:)
      real(8), allocatable :: sval(:)
      complex(8), allocatable :: projm(:,:), lsvec(:,:), rsvec(:,:), auxmat(:,:), auxmat2(:,:), Y(:,:,:)

      minit    = input%properties%wannier%grouparray( wf_group)%group%minitdis
      maxit    = input%properties%wannier%grouparray( wf_group)%group%maxitdis
      gradnorm = input%properties%wannier%grouparray( wf_group)%group%epsdis
      minstep  = input%properties%wannier%grouparray( wf_group)%group%minstepdis
      memlen   = input%properties%wannier%grouparray( wf_group)%group%memlendis

      write( wf_info, '(" disentangle optimal subspace...")')
      call timesec( t0)

      !****************************
      !* PREPARATION
      !****************************
      ! reorder matrices such that frozen bands are at the end
      call wfomega_shuffle(1)
      ! get initial subspace
      allocate( Y( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_kset%nkpt))
      allocate( projm( wf_groups( wf_group)%nst, wf_groups( wf_group)%nproj))
      allocate( map( wf_groups( wf_group)%nwf, wf_kset%nkpt))
      allocate( lsvec( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst), &
                rsvec( wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf), &
                sval( wf_groups( wf_group)%nst))
      allocate( auxmat( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf))
      allocate( auxmat2( wf_groups( wf_group)%nst, wf_groups( wf_group)%nst))
      allocate( dy(2,wf_kset%nkpt))
      wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, &
                    wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, :) = zzero
      Y = zzero; map = 0
      do ik = 1, wf_kset%nkpt
        nok = wf_groups( wf_group)%win_no( ik)
        nik = wf_groups( wf_group)%win_ni( ik)
        dy(1,ik) = nok
        dy(2,ik) = wf_groups( wf_group)%nwf - nik
        do n = 1, nik
          map(n,ik) = dy(2,ik) + n + wf_groups( wf_group)%fwf - 1
        end do
        do n = 1, dy(2,ik)
          map(nik+n,ik) = n + wf_groups( wf_group)%fwf - 1
        end do

        projm = zzero
        do n = 0, nik-1
          wf_transform( wf_groups( wf_group)%fst+dy(1,ik)+n, wf_groups( wf_group)%fwf+dy(2,ik)+n, ik) = zone
          Y( dy(1,ik)+n+1, dy(2,ik)+n+1, ik) = zone
        end do
        
        if( nok .eq. 0) cycle

        do n = 1, wf_groups( wf_group)%win_no( ik)
          projm( n, :) = wf_groups( wf_group)%projection( wf_groups( wf_group)%win_io( n, ik), :, ik)
        end do
        do n = 1, wf_groups( wf_group)%win_ni( ik)
          projm( nok+n, :) = wf_groups( wf_group)%projection( wf_groups( wf_group)%win_ii( n, ik), :, ik)
        end do
        
        auxmat = zzero
        call zgemm( 'n', 'n', nok+nik, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nproj, zone, &
               projm, wf_groups( wf_group)%nst, &
               wf_opf, wf_groups( wf_group)%nproj, zzero, & 
               auxmat, wf_groups( wf_group)%nst)
        if( nik .eq. 0) then
          call svd_divide_conquer( auxmat( 1:nok, :), &
                 sval( 1:wf_groups( wf_group)%nwf), &
                 lsvec( 1:nok, 1:nok), &
                 rsvec)
          call zgemm( 'n', 'n', nok, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, zone, &
                 lsvec, wf_groups( wf_group)%nst, &
                 rsvec, wf_groups( wf_group)%nwf, zzero, &
                 Y(1,1,ik), wf_groups( wf_group)%nst)
        else
          call zgemm( 'n', 'c', nok+nik, nok+nik, wf_groups( wf_group)%nwf, zone, &
                 auxmat, wf_groups( wf_group)%nst, &
                 auxmat, wf_groups( wf_group)%nst, zzero, &
                 auxmat2, wf_groups( wf_group)%nst)

          lsvec = zzero
          call zhediag( auxmat2( 1:nok, 1:nok), sval( 1:nok), evec=lsvec( 1:nok, 1:nok))
          do n = 0, wf_groups( wf_group)%nwf - nik - 1
            Y(1:nok,n+1,ik) = lsvec( 1:nok, nok-n)
          end do
        end if
      end do
      deallocate( projm, lsvec, rsvec, sval, auxmat, auxmat2)
      dyo(1) = wf_groups( wf_group)%nst
      dyo(2) = wf_groups( wf_group)%nwf
      ! get initial spread
      call wfdis_omegai( Y, dyo, wf_kset%nkpt, dy, omegai0)

      !****************************
      !* MINIMIZATION
      !****************************
      convun = 0
      if( input%properties%wannier%grouparray( wf_group)%group%writeconv) then
        call getunit( convun)
        write( convfname, '("dis_conv_",i3.3,".dat")'), wf_group
        open( convun, file=trim( convfname), action='write', form='formatted')
      end if

      if( input%properties%wannier%grouparray( wf_group)%group%optim .eq. 'cg') then
        call manopt_stiefel_cg( Y, dyo, wf_kset%nkpt, dy, &
               cost=wfdis_omegai, &
               grad=wfdis_gradient, &
               epsgrad=gradnorm, minit=minit, maxit=maxit, stdout=convun, minstep=minstep)
      else
        call manopt_stiefel_lbfgs( Y, dyo, wf_kset%nkpt, dy, &
               cost=wfdis_omegai, &
               grad=wfdis_gradient, &
               epsgrad=gradnorm, minit=minit, maxit=maxit, stdout=convun, minstep=minstep, memlen=memlen)
      end if

      if( input%properties%wannier%grouparray( wf_group)%group%writeconv) close( convun)

      !****************************
      !* FINALIZATION
      !****************************
      ! get final spread
      call wfdis_omegai( Y, dyo, wf_kset%nkpt, dy, omegai)
      ! reorder matrices back to original order
      call wfomega_shuffle(-1)
      do ik = 1, wf_kset%nkpt
        wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, ik) = &
          wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, map(:,ik), ik)
      end do
      ! store subspace if needed
      if( trim( wf_groups( wf_group)%method) .eq. 'disSMV') then
        if( allocated( wf_subspace)) deallocate( wf_subspace)
        allocate( wf_subspace( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_kset%nkpt))
        wf_subspace = wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, :)
      end if

      call timesec( t1)

      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"iterations: ",T40,7x,I6)') maxit
      write( wf_info, '(5x,"convergence cutoff: ",T40,E13.6)') input%properties%wannier%grouparray( wf_group)%group%epsdis
      write( wf_info, '(5x,"norm of gradient: ",T40,E13.6)') gradnorm
      write( wf_info, '(5x,"Omega_I before: ",T40,F13.6)') omegai0
      write( wf_info, '(5x,"Omega_I after: ",T40,F13.6)') omegai
      write( wf_info, '(5x,"reduction: ",T40,7x,I5,"%")') nint( 100.d0*(omegai0-omegai)/omegai0)
      write( wf_info, *)
      call flushifc( wf_info)
      deallocate( Y, dy, map)

    end subroutine wfdis_gen

    subroutine wfdis_omegai( y, dyo, ky, dy, omega)
      integer, intent( in)    :: dyo(2), ky, dy(2,ky)
      complex(8), intent( in) :: y(dyo(1),dyo(2),*)
      real(8), intent( out)   :: omega

      integer :: ik

      ! Y to U
      do ik = 1, wf_kset%nkpt
        wf_transform( wf_groups( wf_group)%fst:(wf_groups( wf_group)%fst+dy(1,ik)-1), &
                      wf_groups( wf_group)%fwf:(wf_groups( wf_group)%fwf+dy(2,ik)-1), ik) = y(1:dy(1,ik),1:dy(2,ik),ik)
      end do
      ! compute spread
      call wfomega_gen
      omega = sum( wf_omega_i( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      return
    end subroutine wfdis_omegai

    subroutine wfdis_gradient( y, dyo, ky, dxy, gy, dgo)
      integer, intent( in)     :: dyo(2), ky, dxy(2,ky), dgo(2)
      complex(8), intent( in)  :: y(dyo(1),dyo(2),*)
      complex(8), intent( out) :: gy(dgo(1),dgo(2),*)

      integer :: ik
      complex(8), allocatable :: gu(:,:)

      ! Y to U
      do ik = 1, wf_kset%nkpt
        wf_transform( wf_groups( wf_group)%fst:(wf_groups( wf_group)%fst+dxy(1,ik)-1), &
                      wf_groups( wf_group)%fwf:(wf_groups( wf_group)%fwf+dxy(2,ik)-1), ik) = y(1:dxy(1,ik),1:dxy(2,ik),ik)
      end do
      ! update M matrices
      call wfomega_m
      ! compute gradient
#ifdef USEOMP
!$omp parallel default( shared) private( ik, gu)
#endif
      allocate( gu( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf))
#ifdef USEOMP
!$omp do
#endif
      do ik = 1, wf_kset%nkpt
        ! get the gradient w.r.t. U
        call wfomega_gradiu( ik, gu, wf_groups( wf_group)%nst)
        ! copy relevant part
        gy( 1:dxy(1,ik), 1:dxy(2,ik), ik) = gu( 1:dxy(1,ik), 1:dxy(2,ik))
      end do
#ifdef USEOMP
!$omp end do
#endif
      deallocate( gu)
#ifdef USEOMP
!$omp end parallel
#endif

      return
    end subroutine wfdis_gradient

end module mod_wannier_disentangle
