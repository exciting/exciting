module mod_wannier_omega
  use mod_wannier_variables
  use m_linalg

  implicit none

! methods
  contains

    !=====================================================================================
    ! initializes/updates the M matrices
    subroutine wfomega_m
      integer :: ik, idxn
      complex(8), allocatable :: auxmat(:,:)

      if( .not. allocated( wf_m)) allocate( wf_m( wf_nwf, wf_nwf, wf_kset%nkpt, wf_n_ntot))

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, idxn, auxmat)
#endif
      allocate( auxmat(  wf_nst, wf_nwf))
#ifdef USEOMP
!$OMP DO COLLAPSE( 2)
#endif
      do ik = 1, wf_kset%nkpt
        do idxn = 1, wf_n_ntot 
          call zgemm( 'n', 'n', wf_nst, wf_nwf, wf_nst, zone, &
                 wf_m0( :, :, ik, idxn), wf_nst, &
                 wf_transform( :, :, wf_n_ik( idxn, ik)), wf_nst, zzero, &
                 auxmat, wf_nst)
          call zgemm( 'c', 'n', wf_nwf, wf_nwf, wf_nst, zone, &
                 wf_transform( :, :, ik), wf_nst, &
                 auxmat, wf_nst, zzero, &
                 wf_m( 1, 1, ik, idxn), wf_nwf)
        end do
      end do
#ifdef USEOMP
!$OMP END DO
#endif
      deallocate( auxmat)
#ifdef USEOMP
!$OMP END PARALLEL
#endif
      return
    end subroutine wfomega_m

    !=====================================================================================
    ! initializes/updates the group diagonal M matrices
    subroutine wfomega_m_diag
      integer :: ik, idxn
      complex(8), allocatable :: auxmat(:,:)

      if( .not. allocated( wf_m)) allocate( wf_m( wf_nwf, wf_nwf, wf_kset%nkpt, wf_n_ntot))

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, idxn, auxmat)
#endif
      allocate( auxmat(  wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf))
#ifdef USEOMP
!$OMP DO COLLAPSE( 2)
#endif
      do ik = 1, wf_kset%nkpt
        do idxn = 1, wf_n_ntot 
          call zgemm( 'n', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nst, zone, &
                 wf_m0( wf_groups( wf_group)%fst, wf_groups( wf_group)%fst, ik, idxn), wf_nst, &
                 wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, wf_n_ik( idxn, ik)), wf_nst, zzero, &
                 auxmat, wf_groups( wf_group)%nst)
          call zgemm( 'c', 'n', wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nst, zone, &
                 wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, ik), wf_nst, &
                 auxmat, wf_groups( wf_group)%nst, zzero, &
                 wf_m( wf_groups( wf_group)%fwf, wf_groups( wf_group)%fwf, ik, idxn), wf_nwf)
        end do
      end do
#ifdef USEOMP
!$OMP END DO
#endif
      deallocate( auxmat)
#ifdef USEOMP
!$OMP END PARALLEL
#endif
      return
    end subroutine wfomega_m_diag

    !=====================================================================================
    ! calculates the spread and WF centers
    subroutine wfomega_gen( totonly)
      logical, optional, intent( in) :: totonly

      integer :: ik, j, k, idxn
      real(8), allocatable :: logsum(:,:), log2sum(:,:), abssum(:,:), abs2sum(:,:)
      real(8) :: tmp, p1, p2
      complex(8) :: m
      logical :: tot

      tot = .false.
      if( present( totonly)) tot = totonly
      
      allocate( logsum(  wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_n_ntot))
      allocate( log2sum( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_n_ntot))
      allocate( abssum(  wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_n_ntot))
      allocate( abs2sum( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf, wf_n_ntot))

      !if( .not. wf_initialized) call wannier_init
      if( .not. allocated( wf_m0)) then
        write(*,*)
        write(*, '("Error (mlwf_loc): Matrix elements not available.")')
        stop
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

      call wfomega_m

      logsum = 0.d0
      log2sum = 0.d0
      abssum = 0.d0
      abs2sum = 0.d0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( j, idxn, ik, tmp, k, m, p1, p2)
!$OMP DO COLLAPSE(2)
#endif
      do idxn = 1, wf_n_ntot
        do j = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
          do ik = 1, wf_kset%nkpt
            tmp = 0.d0
            if( abs( wf_m( j, j, ik, idxn)) .gt. 1.d-14) then
              m = wf_m( j, j, ik, idxn)
              p1 = wf_phases( j, ik)
              p2 = wf_phases( j, wf_n_ik( idxn, ik))
              tmp = atan2( aimag( m), dble( m)) - p1 + p2
              logsum( j, idxn) = logsum( j, idxn) + tmp
              log2sum( j, idxn) = log2sum( j, idxn) + tmp*tmp
              abssum( j, idxn) = abssum( j, idxn) + dble( m*conjg( m))
            end if
            if( .not. tot) then
              do k = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
                abs2sum( j, idxn) = abs2sum( j, idxn) + 0.5d0*( dble( wf_m( j, k, ik, idxn)*conjg( wf_m( j, k, ik, idxn))) + &
                                                                dble( wf_m( k, j, ik, idxn)*conjg( wf_m( k, j, ik, idxn))))
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

      if( tot) then
        wf_omega( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf) = 0.d0

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
        wf_omega(    wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf) = 0.d0
        wf_omega_i(  wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf) = 0.d0
        wf_omega_d(  wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf) = 0.d0
        wf_omega_od( wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf) = 0.d0
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

      deallocate( logsum, log2sum, abssum, abs2sum)
      return
    end subroutine wfomega_gen

    !=====================================================================================
    ! calculates the gradient of the spread w.r.t. the matrices U
    subroutine wfomega_gradu( ik, g, ldg)
      use constants, only: zi 
      integer, intent( in)     :: ik, ldg
      complex(8), intent( out) :: g(ldg,*)

      integer :: idxn, ist, jst
      real(8) :: p1, p2
      complex(8) :: z, m0, m
      complex(8), allocatable :: auxmat(:,:)

      allocate( auxmat( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf))
      g(:,1:wf_groups( wf_group)%nwf) = zzero
      do idxn = 1, wf_n_ntot
        ! positive neighbor
        call zgemm( 'n', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nst, zone, &
               wf_m0( wf_groups( wf_group)%fst, wf_groups( wf_group)%fst, ik, idxn), wf_nst, &
               wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, wf_n_ik( idxn, ik)), wf_nst, zzero, &
               auxmat, wf_groups( wf_group)%nst)
        do ist = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
          jst = ist - wf_groups( wf_group)%fwf + 1
          p1 = wf_phases( ist, ik)
          p2 = wf_phases( ist, wf_n_ik( idxn, ik))
          m = wf_m( ist, ist, ik, idxn)
          m0 = m*exp( zi*(p1-p2))
          z = conjg(m0) + zi*(atan2( aimag(m), dble(m)) - p1 + p2 + dot_product( wf_centers(:,ist), wf_n_vc(:,idxn)))/m0
          z = -(4.d0/dble( wf_kset%nkpt))*wf_n_wgt( idxn)*z
          g(:,jst) = g(:,jst) + z*auxmat(:,jst)
        end do
        ! negative neighbor
        call zgemm( 'c', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nst, zone, &
               wf_m0( wf_groups( wf_group)%fst, wf_groups( wf_group)%fst, wf_n_ik2( idxn, ik), idxn), wf_nst, &
               wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, wf_n_ik2( idxn, ik)), wf_nst, zzero, &
               auxmat, wf_groups( wf_group)%nst)
        do ist = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
          jst = ist - wf_groups( wf_group)%fwf + 1
          p1 = wf_phases( ist, wf_n_ik2( idxn, ik))
          p2 = wf_phases( ist, ik)
          m = conjg( wf_m( ist, ist, wf_n_ik2( idxn, ik), idxn))
          m0 = m*exp( zi*(p1-p2))
          z = conjg(m0) + zi*(atan2( aimag(m), dble(m)) - p1 + p2 - dot_product( wf_centers(:,ist), wf_n_vc(:,idxn)))/m0
          z = -(4.d0/dble( wf_kset%nkpt))*wf_n_wgt( idxn)*z
          g(:,jst) = g(:,jst) + z*auxmat(:,jst)
        end do
      end do
      deallocate( auxmat)
      return      
    end subroutine wfomega_gradu

    !=====================================================================================
    ! calculates the gradient of the gauge independent part of the spread w.r.t. the matrices U
    subroutine wfomega_gradiu( ik, g, ldg)
      integer, intent( in)     :: ik, ldg
      complex(8), intent( out) :: g(ldg,*)

      integer :: idxn
      complex(8) :: z
      complex(8), allocatable :: auxmat(:,:)

      allocate( auxmat( wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf))
      g(:,1:wf_groups( wf_group)%nwf) = zzero
      do idxn = 1, wf_n_ntot
        z = cmplx( -4.d0*wf_n_wgt( idxn)/dble( wf_kset%nkpt), 0.d0, 8)
        ! positive neighbor
        call zgemm( 'n', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nst, zone, &
               wf_m0( wf_groups( wf_group)%fst, wf_groups( wf_group)%fst, ik, idxn), wf_nst, &
               wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, wf_n_ik( idxn, ik)), wf_nst, zzero, &
               auxmat, wf_groups( wf_group)%nst)
        call zgemm( 'n', 'c', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, z, &
               auxmat, wf_groups( wf_group)%nst, &
               wf_m( wf_groups( wf_group)%fwf, wf_groups( wf_group)%fwf, ik, idxn), wf_nwf, zone, &
               g(1,1), ldg)
        ! negative neighbor
        call zgemm( 'c', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nst, zone, &
               wf_m0( wf_groups( wf_group)%fst, wf_groups( wf_group)%fst, wf_n_ik2( idxn, ik), idxn), wf_nst, &
               wf_transform( wf_groups( wf_group)%fst, wf_groups( wf_group)%fwf, wf_n_ik2( idxn, ik)), wf_nst, zzero, &
               auxmat, wf_groups( wf_group)%nst)
        call zgemm( 'n', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nwf, wf_groups( wf_group)%nwf, z, &
               auxmat, wf_groups( wf_group)%nst, &
               wf_m( wf_groups( wf_group)%fwf, wf_groups( wf_group)%fwf, wf_n_ik2( idxn, ik), idxn), wf_nwf, zone, &
               g(1,1), ldg)
      end do
      deallocate( auxmat)
      return      
    end subroutine wfomega_gradiu

    !=====================================================================================
    ! finds a unitary diagonal matrix (phase factors) that minimize
    ! the diagonal part of the spread
    ! (i.e. finds also the optimal sheet of the logarithm)
    subroutine wfomega_diagphases( u, ldu1, ldu2, nst)
      integer, intent( in)       :: ldu1, ldu2, nst
      complex(8), intent( inout) :: u(ldu1,ldu2,*)
      integer :: ik, ist, jst, idxn
      real(8) :: swgt
      complex(8) :: z
      integer, allocatable :: ipiv(:)
      real(8), allocatable :: m(:,:), p(:,:)

      allocate( m( wf_kset%nkpt, wf_kset%nkpt))
      allocate( p( wf_kset%nkpt, wf_groups( wf_group)%fwf:wf_groups( wf_group)%lwf))
      allocate( ipiv( wf_kset%nkpt))

      swgt = 2.d0*sum( wf_n_wgt( 1:wf_n_ntot))
      m = 0.d0
      do ik = 1, wf_kset%nkpt
        do idxn = 1, wf_n_ntot
          m( ik, ik) = m( ik, ik) + 2.d0*wf_n_wgt( idxn)
          m( ik, wf_n_ik(  idxn, ik)) = m( ik, wf_n_ik(  idxn, ik)) - wf_n_wgt( idxn)
          m( ik, wf_n_ik2( idxn, ik)) = m( ik, wf_n_ik2( idxn, ik)) - wf_n_wgt( idxn)
        end do
      end do
      p = 0.d0
      do ist = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
        do ik = 1, wf_kset%nkpt
          do idxn = 1, wf_n_ntot
            z = wf_m( ist, ist, ik, idxn)
            p(ik,ist) = p(ik,ist) + wf_n_wgt( idxn)*atan2(  aimag(z), dble(z))
            z = wf_m( ist, ist, wf_n_ik2( idxn, ik), idxn)
            p(ik,ist) = p(ik,ist) + wf_n_wgt( idxn)*atan2( -aimag(z), dble(z))
          end do
        end do
      end do
      call dgesv( wf_kset%nkpt, wf_groups( wf_group)%nwf, &
             m, wf_kset%nkpt, ipiv, &
             p, wf_kset%nkpt, ist)
      if( ist .ne. 0) then
        write(*,*)
        write(*,'("Error (wfomega_fixphases): DGESV returned non-zero info:",i4)') ist
        stop
      end if
            
      do ik = 1, wf_kset%nkpt
        do ist = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
          jst = ist - wf_groups( wf_group)%fwf + 1
          z = cmplx( cos( p(ik,ist)), sin( p(ik,ist)), 8)
          u(1:nst,jst,ik) = z*u(1:nst,jst,ik)
          wf_phases( ist, ik) = p(ik,ist)
        end do
      end do
      deallocate( m, p, ipiv)
      return      
    end subroutine wfomega_diagphases

    !=====================================================================================
    ! reorders the U and M0 matrices such that the rows are in the order:
    ! outer window bands, inner window bands, bands not used
    subroutine wfomega_shuffle( dir)
      integer, intent( in) :: dir

      integer :: ik, idxn, ist
      integer, allocatable :: map(:)

      if( .not.(wf_groups( wf_group)%method .eq. 'disSMV' .or. wf_groups( wf_group)%method .eq. 'disFull') .or. (dir .eq. 0)) return

      allocate( map( wf_groups( wf_group)%nst))
      do ik = 1, wf_kset%nkpt
        ! build map
        call getmap( ik)
        ! shuffle rows of transformation matrices
        do ist = wf_groups( wf_group)%fwf, wf_groups( wf_group)%lwf
          wf_transform( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, ist, ik) = wf_transform(map,ist,ik)
        end do
        ! shuffle rows of subspace
        if( allocated( wf_subspace)) then
          do ist = 1, wf_groups( wf_group)%nwf
            wf_subspace(:,ist,ik) = wf_subspace(map-wf_groups( wf_group)%fst+1,ist,ik)
          end do
        end if
        ! shuffle rows of matrix elements
        do idxn = 1, wf_n_ntot
          do ist = wf_groups( wf_group)%fst, wf_groups( wf_group)%lst
            wf_m0( wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, ist, ik, idxn) = wf_m0(map,ist,ik,idxn)
          end do
        end do
        ! shuffle columns of matrix elements
        do idxn = 1, wf_n_ntot
          call getmap( wf_n_ik( idxn, ik))
          do ist = wf_groups( wf_group)%fst, wf_groups( wf_group)%lst
            wf_m0( ist, wf_groups( wf_group)%fst:wf_groups( wf_group)%lst, ik, idxn) = wf_m0(ist,map,ik,idxn)
          end do
        end do
      end do
      deallocate( map)
      return

      contains
        subroutine getmap( k)
          use sorting, only: sort_index_1d
          integer, intent( in) :: k
          integer :: nik, nok, i, j
          integer, allocatable :: idx(:)
          map = 0
          nik = wf_groups( wf_group)%win_ni( k)
          nok = wf_groups( wf_group)%win_no( k)
          do i = 1, nok
            map(i) = wf_groups( wf_group)%win_io(i,k)
          end do
          do i = 1, nik
            map( nok+i) = wf_groups( wf_group)%win_ii(i,k)
          end do
          i = nik + nok
          do j = wf_groups( wf_group)%fst, wf_groups( wf_group)%lst
            if( any( map .eq. j)) cycle
            i = i + 1
            map(i) = j
          end do
          ! invert map if necessary
          if( dir .lt. 0) then
            allocate( idx( wf_groups( wf_group)%nst))
            idx = sort_index_1d( wf_groups( wf_group)%nst, map)
            map = idx + wf_groups( wf_group)%fst - 1
            deallocate( idx)
          end if
        end subroutine getmap
    end subroutine wfomega_shuffle

end module mod_wannier_omega
