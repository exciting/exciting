module mod_wannier_interpolate
  use modmain
  use mod_wannier
  use m_linalg
  use m_plotmat
  use mod_opt_tetra
  implicit none

! module variables
  type( k_set)  :: wfint_kset                       ! k-point set on which the interpolation is performed
  type( Gk_set) :: wfint_Gkset                      ! G+k-point set on which the interpolation is performed
  logical	:: wfint_initialized = .false.
  logical	:: wfint_mindist                    ! use minimum distances interpolation method
  real(8)	:: wfint_efermi                     ! interpolated fermi energy
  type( t_set)  :: wfint_tetra                      ! tetrahedra for tetrahedron integration
  real(8)       :: wfint_vvbm(3), wfint_evbm
  real(8)       :: wfint_vcbm(3), wfint_ecbm

  real(8), allocatable    :: wfint_eval(:,:)        ! interpolated eigenenergies
  complex(8), allocatable :: wfint_transform(:,:,:) ! corresponding expansion coefficients
  integer, allocatable    :: wfint_bandmap(:)       ! map from interpolated bands to original bands

  real(8), allocatable    :: wfint_occ(:,:)         ! interpolated occupation numbers
  real(8), allocatable    :: wfint_phase(:,:)       ! summed phase factors in interpolation
  complex(8), allocatable :: wfint_pqr(:,:)

  ! often used matrices
  complex(8), allocatable :: wfint_hwr(:,:,:)       ! R-dependent Hamiltonian in Wannier gauge

! methods
  contains
    !BOP
    ! !ROUTINE: wfint_init
    ! !INTERFACE:
    !
    subroutine wfint_init( int_kset, evalin, serial)
      ! !USES:
        use m_getunit 
      ! !INPUT PARAMETERS:
      !   int_kset : k-point set on which the interpolation is performed on (in, type k_set)
      ! !DESCRIPTION:
      !   Sets up the interpolation grid and calculates the interpolated eigenenergies as well as 
      !   the corresponding expansion coefficients and phasefactors for the interpolated wavefunctions.
      !
      ! !REVISION HISTORY:
      !   Created July 2017 (SeTi)
      !EOP
      !BOC
      type( k_set), intent( in)      :: int_kset
      real(8), optional, intent( in) :: evalin( wf_fst:wf_lst, wf_kset%nkpt)
      logical, optional, intent( in) :: serial
    
      integer :: fst, lst
      real(8), allocatable :: evalfv(:,:), evalin_(:,:)

      if( wfint_initialized) call wfint_destroy

      wfint_mindist = input%properties%wannier%mindist
      wfint_kset = int_kset
    
      allocate( evalin_( wf_fst:wf_lst, wf_kset%nkpt))
    
      if( present( evalin)) then
        evalin_ = evalin
      else
        call wfhelp_geteval( evalfv, fst, lst)
        if( (fst .gt. wf_fst) .or. (lst .lt. wf_lst)) then
          write(*,*)
          write(*, '("Error (wfint_init): Eigenenergies out of range of Wannier functions.")')
          stop
        end if
        evalin_ = evalfv( wf_fst:wf_lst, :)
        deallocate( evalfv)
      end if
    
      if( wfint_mindist) call wannier_mindist( input%structure%crystal%basevect, wf_kset, wf_nrpt, wf_rvec, wf_nwf, wf_centers, wf_wdistvec, wf_wdistmul)
      call wfint_fourierphases( wfint_kset, wf_nrpt, wf_rvec, wf_rmul, wfint_pqr, wfint_phase)

      allocate( wfint_eval( wf_nwf, wfint_kset%nkpt))
      allocate( wfint_transform( wf_nwf, wf_nwf, wfint_kset%nkpt))
      
      if( present( serial)) then
        call wfint_interpolate_eigsys( evalin_, serial=serial)
      else
        call wfint_interpolate_eigsys( evalin_)
      end if
    
      wfint_initialized = .true.
      deallocate( evalin_)

      return
    end subroutine wfint_init
    !EOC

!--------------------------------------------------------------------------------------
    
    subroutine wfint_destroy
      if( allocated( wfint_phase))     deallocate( wfint_phase)
      if( allocated( wfint_pqr))       deallocate( wfint_pqr)
      if( allocated( wfint_eval))      deallocate( wfint_eval)
      if( allocated( wfint_transform)) deallocate( wfint_transform)
      if( allocated( wfint_occ))       deallocate( wfint_occ)
      if( allocated( wfint_hwr))       deallocate( wfint_hwr)
      !if( allocated( wfint_bandmap))   deallocate( wfint_bandmap)
      !if( allocated( wfint_rvec))      deallocate( wfint_rvec)
      !if( allocated( wfint_rmul))      deallocate( wfint_rmul)
      !if( allocated( wfint_wdistvec))  deallocate( wfint_wdistvec)
      !if( allocated( wfint_wdistmul))  deallocate( wfint_wdistmul)
      call delete_k_vectors( wfint_kset)
      call delete_Gk_vectors( wfint_Gkset)
      wfint_initialized = .false.
      return
    end subroutine wfint_destroy

!--------------------------------------------------------------------------------------
!   BASIC FUNCTIONS OFTEN NEEDED IN GENERIC INTERPOLATION
!--------------------------------------------------------------------------------------
    
    ! precomputes Fourier factors appearing in generic interpolation
    subroutine wfint_fourierphases( qset, nrpt, rvec, rmul, pqr, phase)
      type( k_set), intent( in)                    :: qset
      integer, intent( in)                         :: nrpt
      integer, intent( in)                         :: rvec( 3, nrpt)
      integer, intent( in)                         :: rmul( nrpt)
      complex(8), allocatable, intent( out)        :: pqr(:,:)
      real(8), allocatable, optional, intent( out) :: phase(:,:)

      integer :: ir, iq
      real(8) :: d
      complex(8), allocatable :: auxmat(:,:)

      if( allocated( pqr)) deallocate( pqr)
      allocate( pqr( qset%nkpt, nrpt))

      do ir = 1, nrpt
        do iq = 1, qset%nkpt
          d = twopi*dot_product( qset%vkl( :, iq), dble( rvec( :, ir)))
          pqr( iq, ir) = cmplx( cos( d), sin( d), 8)/dble( rmul( ir))
        end do
      end do

      if( present( phase)) then
        if( allocated( phase)) deallocate( phase)
        allocate( phase( wf_kset%nkpt, qset%nkpt))
        allocate( auxmat( qset%nkpt, wf_kset%nkpt))
        call zgemm( 'n', 't', qset%nkpt, wf_kset%nkpt, nrpt, zone, &
               pqr, qset%nkpt, &
               wf_pkr, wf_kset%nkpt, zzero, &
               auxmat, qset%nkpt)
        phase = dble( transpose( auxmat))/wf_kset%nkpt
        deallocate( auxmat)
      end if

      return
    end subroutine wfint_fourierphases

!--------------------------------------------------------------------------------------

    ! generates real-space Hamiltonian in Wannier gauge H_{mn}(R)
    subroutine wfint_genhwr( evalin)
      use m_plotmat
      real(8), intent( in) :: evalin( wf_fst:wf_lst, wf_kset%nkpt)

      integer :: ik, ist, jst, igroup
      complex(8), allocatable :: hwk(:,:,:), auxmat(:,:)

      if( allocated( wfint_hwr)) deallocate( wfint_hwr)
      allocate( wfint_hwr( wf_nwf, wf_nwf, wf_nrpt))

      allocate( hwk( wf_nwf, wf_nwf, wf_kset%nkpt), auxmat( wf_fst:wf_lst, wf_nwf))
      hwk = zzero
      do ik = 1, wf_kset%nkpt
        do ist = wf_fst, wf_lst
          do jst = 1, wf_nwf
            auxmat( ist, jst) = wf_transform( ist, jst, ik)*evalin( ist, ik)
          end do
        end do
        do igroup = 1, wf_ngroups
          call zgemm( 'c', 'n', wf_groups( igroup)%nwf, wf_groups( igroup)%nwf, wf_nst, zone, &
                 wf_transform( :, wf_groups( igroup)%fwf, ik), wf_nst, &
                 auxmat( :, wf_groups( igroup)%fwf), wf_nst, zzero, &
                 hwk( wf_groups( igroup)%fwf, wf_groups( igroup)%fwf, ik), wf_nwf)
        end do
      end do
      call wfint_ftk2r( hwk, wfint_hwr, 1)
      deallocate( hwk, auxmat)

      return
    end subroutine wfint_genhwr

!--------------------------------------------------------------------------------------

    ! transform matrix from coarse reciprocal-space grid to real-space grid
    subroutine wfint_ftk2r( mwk, mwr, ld)
      complex(8), intent( in)  :: mwk( wf_nwf, wf_nwf, ld, wf_kset%nkpt)
      complex(8), intent( out) :: mwr( wf_nwf, wf_nwf, ld, wf_nrpt)
      integer, intent( in)     :: ld

      call zgemm( 'n', 'n', wf_nwf*wf_nwf*ld, wf_nrpt, wf_kset%nkpt, cmplx( 1.d0/wf_kset%nkpt, 0.d0, 8), &
             mwk, wf_nwf*wf_nwf*ld, &
             wf_pkr, wf_kset%nkpt, zzero, &
             mwr, wf_nwf*wf_nwf*ld)
      return
    end subroutine wfint_ftk2r

!--------------------------------------------------------------------------------------

    ! transform matrix from double coarse reciprocal-space grid to double real-space grid
    subroutine wfint_ftkk2rr( mwkk, mwrr, ld)
      complex(8), intent( in)  :: mwkk( wf_nwf, wf_nwf, ld, wf_kset%nkpt, wf_kset%nkpt)
      complex(8), intent( out) :: mwrr( wf_nwf, wf_nwf, ld, wf_nrpt, wf_nrpt)
      integer, intent( in)     :: ld

      integer :: ik1, ik2, ir1, ir2
      complex(8) :: phase( wf_kset%nkpt, wf_kset%nkpt)

#ifdef USEOMP
!$omp parallel default( shared) private( ir1, ir2, ik1, ik2, phase)
!$omp do collapse(2)
#endif
      do ir2 = 1, wf_nrpt
        do ir1 = 1, wf_nrpt
          do ik2 = 1, wf_kset%nkpt
            do ik1 = 1, wf_kset%nkpt
              phase(ik1,ik2) = wf_pkr( ik1, ir1)*wf_pkr( ik2, ir2)
            end do
          end do
          call zgemv( 'n', wf_nwf*wf_nwf*ld, wf_kset%nkpt*wf_kset%nkpt, cmplx( 1.d0/wf_kset%nkpt**2, 0.d0, 8), &
                 mwkk, wf_nwf*wf_nwf*ld, &
                 phase, 1, zzero, &
                 mwrr(1,1,1,ir1,ir2), 1)
        end do
      end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif

      return
    end subroutine wfint_ftkk2rr

!--------------------------------------------------------------------------------------

    ! transform matrix from real-space grid to fine reciprocal-space grid
    ! R-dependent summation weight can be provided
    subroutine wfint_ftr2q( mwr, mwq, iq, ld, rwgt)
      complex(8), intent( in)  :: mwr( wf_nwf, wf_nwf, ld, wf_nrpt)
      complex(8), intent( out) :: mwq( wf_nwf, wf_nwf, ld)
      integer, intent( in)     :: iq, ld
      complex(8), optional, intent( in) :: rwgt( wf_nrpt)

      integer :: ir, ist, jst
      complex(8), allocatable :: wgts(:,:)

      if( wfint_mindist) then
        mwq = zzero
        allocate( wgts( wf_nwf, wf_nwf))
        do ir = 1, wf_nrpt
          call wfint_getmindistwgts( ir, iq, wgts)
          do ist = 1, wf_nwf
            do jst = 1, wf_nwf
              if( present( rwgt)) then
                mwq( ist, jst, :) = mwq( ist, jst, :) + rwgt( ir)*wgts( ist, jst)*mwr( ist, jst, :,  ir)
              else
                mwq( ist, jst, :) = mwq( ist, jst, :) + wgts( ist, jst)*mwr( ist, jst, :, ir)
              end if
            end do
          end do
        end do
        deallocate( wgts)
      else
        if( present( rwgt)) then
          call zgemv( 'n', wf_nwf*wf_nwf*ld, wf_nrpt, zone, &
                 mwr, wf_nwf*wf_nwf*ld, &
                 wfint_pqr( iq, :)*rwgt, 1, zzero, &
                 mwq, 1)
        else
          call zgemv( 'n', wf_nwf*wf_nwf*ld, wf_nrpt, zone, &
                 mwr, wf_nwf*wf_nwf*ld, &
                 wfint_pqr( iq, :), 1, zzero, &
                 mwq, 1)
        end if
      end if
      return
    end subroutine wfint_ftr2q

!--------------------------------------------------------------------------------------

    ! transform matrix from double real-space grid to double fine reciprocal-space grid
    ! R-dependent summation weight can be provided
    subroutine wfint_ftrr2qq( mwrr, mwqq, iq1, iq2, ld, rwgt)
      complex(8), intent( in)  :: mwrr( wf_nwf, wf_nwf, ld, wf_nrpt, wf_nrpt)
      complex(8), intent( out) :: mwqq( wf_nwf, wf_nwf, ld)
      integer, intent( in)     :: iq1, iq2, ld
      complex(8), optional, intent( in) :: rwgt( wf_nrpt, wf_nrpt)

      integer :: ir1, ir2
      complex(8) :: phase( wf_nrpt, wf_nrpt)

      do ir1 = 1, wf_nrpt
        do ir2 = 1, wf_nrpt
          phase( ir1, ir2) = wfint_pqr( iq1, ir1)*wfint_pqr( iq2, ir2)
        end do
      end do

      if( present( rwgt)) then
        call zgemv( 'n', wf_nwf*wf_nwf*ld, wf_nrpt*wf_nrpt, zone, &
               mwrr, wf_nwf*wf_nwf*ld, &
               phase*rwgt, 1, zzero, &
               mwqq, 1)
      else
        call zgemv( 'n', wf_nwf*wf_nwf*ld, wf_nrpt*wf_nrpt, zone, &
               mwrr, wf_nwf*wf_nwf*ld, &
               phase, 1, zzero, &
               mwqq, 1)
      end if
      return
    end subroutine wfint_ftrr2qq

!--------------------------------------------------------------------------------------

    ! transform matrix at k from Hamiltonian to Wannier gauge
    subroutine wfint_h2wk( mhk, mwk, ik, ld, ik2)
      complex(8), intent( in)        :: mhk( wf_fst:wf_lst, wf_fst:wf_lst, ld)
      complex(8), intent( out)       :: mwk( wf_nwf, wf_nwf, ld)
      integer, intent( in)           :: ik, ld
      integer, optional, intent( in) :: ik2

      integer :: d, ip
      complex(8) :: auxmat( wf_nwf, wf_fst:wf_lst, ld)

      ip = ik
      if( present( ik2)) ip = ik2

      call zgemm( 'c', 'n', wf_nwf, wf_nst*ld, wf_nst, zone, &
             wf_transform( :, :, ik), wf_nst, &
             mhk, wf_nst, zzero, &
             auxmat, wf_nwf)
      do d = 1, ld
        call zgemm( 'n', 'n', wf_nwf, wf_nwf, wf_nst, zone, &
               auxmat( :, :, d), wf_nwf, &
               wf_transform( :, :, ip), wf_nst, zzero, &
               mwk( :, :, d), wf_nwf)
      end do
      return
    end subroutine wfint_h2wk

!--------------------------------------------------------------------------------------

    ! transform matrix at q from Wannier to Hamiltonian gauge
    subroutine wfint_w2hq( mwq, mhq, iq, ld, iq2)
      complex(8), intent( in)        :: mwq( wf_nwf, wf_nwf, ld)
      complex(8), intent( out)       :: mhq( wf_nwf, wf_nwf, ld)
      integer, intent( in)           :: iq, ld
      integer, optional, intent( in) :: iq2

      integer :: d, ip
      complex(8) :: auxmat( wf_nwf, wf_nwf, ld)

      ip = iq
      if( present( iq2)) ip = iq2

      call zgemm( 'c', 'n', wf_nwf, wf_nwf*ld, wf_nwf, zone, &
             wfint_transform( :, :, iq), wf_nwf, &
             mwq, wf_nwf, zzero, &
             auxmat, wf_nwf)
      do d = 1, ld
        call zgemm( 'n', 'n', wf_nwf, wf_nwf, wf_nwf, zone, &
               auxmat( :, :, d), wf_nwf, &
               wfint_transform( :, :, ip), wf_nwf, zzero, &
               mhq( :, :, d), wf_nwf)
      end do
      return
    end subroutine wfint_w2hq

!--------------------------------------------------------------------------------------

    ! get weights needed for minimal distance mode interpolation
    subroutine wfint_getmindistwgts( ir, iq, wgts)
      integer, intent( in)     :: ir, iq
      complex(8), intent( out) :: wgts( wf_nwf, wf_nwf)

      integer :: ist, jst, m
      real(8) :: dotp

      wgts = zzero
      do ist = 1, wf_nwf
        do jst = 1, wf_nwf
          do m = 1, wf_wdistmul( ist, jst, ir)
            dotp = twopi*dot_product( wfint_kset%vkl( :, iq), dble( wf_wdistvec( :, m, ist, jst, ir)))
            dotp = dotp + dot_product( wfint_kset%vkc( :, iq), wf_centers( :, jst) - wf_centers( :, ist))
            wgts( ist, jst) = wgts( ist, jst) + cmplx( cos( dotp), sin( dotp), 8)
          end do
          wgts( ist, jst) = wgts( ist, jst)/wf_rmul( ir)/wf_wdistmul( ist, jst, ir)
        end do
      end do
      return
    end subroutine wfint_getmindistwgts
    
!--------------------------------------------------------------------------------------
! INTERPOLATION FUNCTIONS FOR VARIOUS QUANTITIES
!--------------------------------------------------------------------------------------

    !BOP
    ! !ROUTINE: wfint_interpolate_eigsys
    ! !INTERFACE:
    !
    subroutine wfint_interpolate_eigsys( evalin, serial)
      ! !USES:
      ! !INPUT PARAMETERS:
      !   evalin : eigenenergies on original grid (in, real( wf_fst:wf_lst, wf_kset%nkpt))
      ! !DESCRIPTION:
      !   Calculates the interpolated eigenenergies as well as the corresponding expansion 
      !   coefficients and phasefactors for the interpolated wavefunctions.
      !
      ! !REVISION HISTORY:
      !   Created July 2017 (SeTi)
      !EOP
      !BOC

      real(8), intent( in)           :: evalin( wf_fst:wf_lst, wf_kset%nkpt)
      logical, optional, intent( in) :: serial
      
      integer :: iq, q1, q2
      logical :: parallel

      complex(8), allocatable :: hwq(:,:)

      parallel = .true.
      if( present( serial)) parallel = .not. serial

      !**********************************************
      ! interpolated eigenenergies and corresponding 
      ! eigenvectors in Wannier basis
      !**********************************************
    
      ! get R-dependent Hamiltonian in Wannier gauge
      if( .not. allocated( wfint_hwr)) call wfint_genhwr( evalin)

      allocate( hwq( wf_nwf, wf_nwf))
      if( parallel) then
        q1 = firstofset( mpiglobal%rank, wfint_kset%nkpt)
        q2 = lastofset( mpiglobal%rank, wfint_kset%nkpt)
      else
        q1 = 1
        q2 = wfint_kset%nkpt
      end if
#ifdef USEOMP
!$omp parallel default( shared) private( iq, hwq)
!$omp do
#endif
      do iq = q1, q2
        ! transform it to the q-dependent Hamiltonian in Wannier gauge
        call wfint_ftr2q( wfint_hwr, hwq, iq, 1)
        ! find the matrices that transform it to the Hamilton gauge
        call zhediag( hwq, wfint_eval( :, iq), evec=wfint_transform( :, :, iq))
      end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
      if( parallel) then
        call mpi_allgatherv_ifc( set=wfint_kset%nkpt, rlen=wf_nwf, rbuf=wfint_eval)
        call mpi_allgatherv_ifc( set=wfint_kset%nkpt, rlen=wf_nwf*wf_nwf, zbuf=wfint_transform)
        call barrier
      end if
      deallocate( hwq)
    
      return
    end subroutine wfint_interpolate_eigsys
    !EOC

!--------------------------------------------------------------------------------------
    
    !BOP
    ! !ROUTINE: wfint_interpolate_occupancy
    ! !INTERFACE:
    !
    subroutine wfint_interpolate_occupancy( usetetra)
      ! !USES:
      use mod_opt_tetra
      ! !INPUT PARAMETERS:
      !   usetetra : to use tetrahedron integration or not (in, optional, logical)
      ! !DESCRIPTION:
      !   Calclulates the interpolated occupation numbers for the wannierized bands and
      !   interpolated Fermi energy.
      !
      ! !REVISION HISTORY:
      !   Created July 2017 (SeTi)
      !EOP
      logical, optional, intent( in) :: usetetra
      !BOC
      integer :: fst, lst
      logical :: usetetra_ = .false.

      if( present( usetetra)) usetetra_ = usetetra
      if( allocated( wfint_occ)) deallocate( wfint_occ)
      allocate( wfint_occ( wf_nwf, wfint_kset%nkpt))

      fst = wf_fst
      lst = wf_fst + wf_nwf - 1
      if( allocated( wfint_bandmap)) then
        fst = wfint_bandmap(1)
        lst = wfint_bandmap( wf_nwf)
      end if

      if( usetetra_) then
        call wfhelp_init_tetra( wfint_tetra, wfint_kset, ttype=2, reduce=.true.)
        call wfhelp_occupy( wfint_kset, wfint_eval, fst, lst, wfint_efermi, wfint_occ, tetra=wfint_tetra)
      else
        call wfhelp_occupy( wfint_kset, wfint_eval, fst, lst, wfint_efermi, wfint_occ)
      end if

      return
    end subroutine wfint_interpolate_occupancy
    !EOC

!--------------------------------------------------------------------------------------

    !BOP
    ! !ROUTINE: wfint_interpolate_ederiv
    ! !INTERFACE:
    !
    subroutine wfint_interpolate_ederiv( velo_, mass_)
      ! !INPUT PARAMETERS:
      !   velo : group-velocity vector (out, real( 3, wf_nwf, wfint_kset%nkpt))
      !   mass : effective-mass tensor (out, real( 3, 3, wf_nwf, wfint_kset%nkpt))
      ! !DESCRIPTION:
      ! Calculates the group-velocity vector (first derivative) and 
      ! effective-mass tensors (inverse of second derivative) for the Wannier-interpolated bands.
      !
      ! !REVISION HISTORY:
      !   Created July 2017 (SeTi)
      !EOP
      real(8), intent( out) :: velo_( 3, wf_nwf, wfint_kset%nkpt)
      real(8), intent( out) :: mass_( 3, 3, wf_nwf, wfint_kset%nkpt)
      !BOC

      integer :: iq, ir, ist, jst, im, d1, d2, ndeg, sdeg, ddeg
      real(8) :: dotp, eps1, eps2, vr(3)
      complex(8) :: ftweight, hamwk( wf_nwf, wf_nwf)
      complex(8) :: velo( wf_nwf, wf_nwf, 3, wfint_kset%nkpt)
      complex(8) :: dmat( wf_nwf, wf_nwf, 3, wfint_kset%nkpt)
      complex(8) :: mass( wf_nwf, wf_nwf, 3, 3, wfint_kset%nkpt)
      real(8) :: degeval( wf_nwf)
      complex(8) :: degmat( wf_nwf, wf_nwf), degevec( wf_nwf, wf_nwf)
      
      real(8), allocatable :: evalin(:,:)
      complex(8), allocatable :: auxmat(:,:)

      eps1 = 1.d-4
      eps2 = 1.d-2
      ddeg = 0

      velo_ = 0.d0
      mass_ = 0.d0

      call wfhelp_geteval( evalin, ist, jst)
      if( .not. allocated( wfint_hwr)) call wfint_genhwr( evalin( wf_fst:wf_lst, :))

      ! get R-dependent Hamiltonian in Wannier gauge
      allocate( auxmat( wf_nwf, wf_nwf))

      velo = zzero
      mass = zzero
      dmat = zzero

      ! get first band-derivative (velocity)
      do d1 = 1, 3
        do iq = 1, wfint_kset%nkpt

          hamwk = zzero
          do ir = 1, wf_nrpt
            call r3mv( input%structure%crystal%basevect, dble( wf_rvec( :, ir)), vr)
            if( wfint_mindist) then
              do ist = 1, wf_nwf
                do jst = 1, wf_nwf
                  ftweight = zzero
                  do im = 1, wf_wdistmul( ist, jst, ir)
                    call r3mv( input%structure%crystal%basevect, dble( wf_wdistvec( :, im, ist, jst, ir)), vr)
                    dotp = dot_product( wfint_kset%vkc( :, iq), vr + wf_centers( :, jst) - wf_centers( :, ist))
                    ftweight = ftweight + cmplx( cos( dotp), sin( dotp), 8)
                  end do
                  hamwk( ist, jst) = hamwk( ist, jst) + zi*vr( d1)*ftweight/wf_rmul( ir)/wf_wdistmul( ist, jst, ir)*wfint_hwr( ist, jst, ir)
                end do
              end do
            else
              hamwk = hamwk + zi*vr( d1)*wfint_pqr( iq, ir)*wfint_hwr( :, :, ir)
            end if
          end do
          ! force hermiticity (not guaranteed if wfint_mindist)
          hamwk = cmplx( 0.5d0, 0.d0, 8)*(hamwk + conjg( transpose( hamwk)))

          call wfint_w2hq( hamwk, velo( :, :, d1, iq), iq, 1)
          ! handle degeneracies of first order
          sdeg = 1
          do while( (sdeg .lt. wf_nwf) .and. (d1 .eq. ddeg))
            ndeg = 1
            do while( (sdeg .lt. wf_nwf))
              if( abs( wfint_eval( sdeg, iq) - wfint_eval( sdeg+1, iq)) .ge. eps1) exit
              ndeg = ndeg + 1
              sdeg = sdeg + 1
            end do
            if( ndeg .gt. 1) then
              degmat( 1:ndeg, 1:ndeg) = velo( (sdeg-ndeg+1):sdeg, (sdeg-ndeg+1):sdeg, d1, iq)
              call zhediag( degmat( 1:ndeg, 1:ndeg), degeval( 1:ndeg), evec=degevec( 1:ndeg, 1:ndeg))
              call zgemm( 'n', 'n', wf_nwf, ndeg, ndeg, zone, &
                     wfint_transform( :, (sdeg-ndeg+1):sdeg, iq), wf_nwf, &
                     degevec, wf_nwf, zzero, &
                     auxmat( :, 1:ndeg), wf_nwf)
              wfint_transform( :, (sdeg-ndeg+1):sdeg, iq) = auxmat( :, 1:ndeg)
              call wfint_w2hq( hamwk, velo( :, :, d1, iq), iq, 1)
            else
              sdeg = sdeg + 1
            end if
          end do

          ! set up D-matrix
          do ist = 1, wf_nwf
            do jst = 1, wf_nwf
              if( abs( wfint_eval( ist, iq) - wfint_eval( jst, iq)) .gt. eps1) then
                dmat( ist, jst, d1, iq) = velo( ist, jst, d1, iq)/(wfint_eval( jst, iq) - wfint_eval( ist, iq))
              end if
            end do
          end do

        end do
      end do

      do d1 = 1, 3
        ! get second band-derivative (inverse effective mass)
        do iq = 1, wfint_kset%nkpt
          do d2 = 1, 3
            hamwk = zzero
            do ir = 1, wf_nrpt
              if( wfint_mindist) then
                do ist = 1, wf_nwf
                  do jst = 1, wf_nwf
                    ftweight = zzero
                    do im = 1, wf_wdistmul( ist, jst, ir)
                      call r3mv( input%structure%crystal%basevect, dble( wf_wdistvec( :, im, ist, jst, ir)), vr)
                      dotp = dot_product( wfint_kset%vkc( :, iq), vr + wf_centers( :, jst) - wf_centers( :, ist))
                      ftweight = ftweight + cmplx( cos( dotp), sin( dotp), 8)
                    end do
                    hamwk( ist, jst) = hamwk( ist, jst) - vr( d1)*vr( d2)*ftweight/wf_rmul( ir)/wf_wdistmul( ist, jst, ir)*wfint_hwr( ist, jst, ir)
                  end do
                end do
              else
                call r3mv( input%structure%crystal%basevect, dble( wf_rvec( :, ir)), vr)
                hamwk = hamwk - vr( d1)*vr( d2)*wfint_pqr( iq, ir)*wfint_hwr( :, :, ir)
              end if
            end do
            ! force hermiticity (not guaranteed if wfint_mindist)
            hamwk = cmplx( 0.5d0, 0.d0, 8)*(hamwk + conjg( transpose( hamwk)))
            call wfint_w2hq( hamwk, mass( :, :, d1, d2, iq), iq, 1)
            call zgemm( 'n', 'n', wf_nwf, wf_nwf, wf_nwf, zone, &
                   velo( :, :, d1, iq), wf_nwf, &
                   dmat( :, :, d2, iq), wf_nwf, zzero, &
                   auxmat, wf_nwf)
            mass( :, :, d1, d2, iq) = mass( :, :, d1, d2, iq) + auxmat + conjg( transpose( auxmat))
            ! handle degeneracies of second order
            sdeg = 1
            do while( (sdeg .lt. wf_nwf) .and. (d1 .eq. ddeg) .and. (d2 .eq. ddeg))
              ndeg = 1
              do while( (sdeg .lt. wf_nwf))
                if( (abs( velo( sdeg, sdeg, d1, iq) - velo( sdeg+1, sdeg+1, d1, iq)) .ge. eps2) .and. &
                    (abs( wfint_eval( sdeg, iq) - wfint_eval( sdeg+1, iq)) .ge. eps1)) exit
                ndeg = ndeg + 1
                sdeg = sdeg + 1
              end do
              if( ndeg .gt. 1) then
                degmat( 1:ndeg, 1:ndeg) = mass( (sdeg-ndeg+1):sdeg, (sdeg-ndeg+1):sdeg, d1, d2, iq)
                call zhediag( degmat( 1:ndeg, 1:ndeg), degeval( 1:ndeg), evec=degevec( 1:ndeg, 1:ndeg))
                do ist = 1, ndeg
                  mass( sdeg-ndeg+ist, sdeg-ndeg+ist, d1, d2, iq) = cmplx( degeval( ist), 0, 8)
                end do
              else
                sdeg = sdeg + 1
              end if
            end do
          end do
        end do
      end do

      velo_ = 0.d0
      mass_ = 0.d0
      do iq = 1, wfint_kset%nkpt
        do ist = 1, wf_nwf
          velo_( :, ist, iq) = dble( velo( ist, ist, :, iq))
          mass_( :, :, ist, iq) = dble( mass( ist, ist, :, :, iq))
        end do
      end do

      deallocate( evalin, auxmat)

      return
    end subroutine wfint_interpolate_ederiv
    !EOC

!--------------------------------------------------------------------------------------

    !BOP
    ! !ROUTINE: wfint_interpolate_bandchar
    ! !INTERFACE:
    !
    subroutine wfint_interpolate_bandchar( lmax, bc)
      ! !INPUT PARAMETERS:
      !   lmax : maximum angular quantum number (in, integer)
      !   bc   : band character (out, real( 0:lmax, natmtot, wf_nwf, wfint_kset%nkpt))
      ! !DESCRIPTION:
      ! Calculates the band character for the Wannier-interpolated bands.
      !
      ! !REVISION HISTORY:
      !   Created July 2017 (SeTi)
      !EOP
      integer, intent( in) :: lmax
      real(8), intent( out) :: bc( 0:lmax, natmtot, wf_nwf, wfint_kset%nkpt)
      !BOC

      integer :: iq, l, m, lm, lmmax, lammax, ist, ias
      real(8), allocatable :: radolp(:,:,:,:), elm(:,:)
      complex(8), allocatable :: dmat(:,:,:,:), radcoeffr(:,:,:,:,:), ulm(:,:,:), auxmat(:,:)

      lmmax = (lmax + 1)**2

      call wfint_gen_radcoeffr( lmax, lammax, radcoeffr)
      call wfint_gen_radolp( lmax, lammax, radolp)

      allocate( elm( lmmax, natmtot))
      allocate( ulm( lmmax, lmmax, natmtot))
      call genlmirep( lmax, lmmax, elm, ulm)

      bc = 0.d0
#ifdef USEOMP
!$omp parallel default( shared) private( iq, dmat, auxmat, l, m, lm, ist, ias)
#endif
      allocate( dmat( lmmax, lmmax, wf_nwf, natmtot))
      allocate( auxmat( lmmax, lmmax))
#ifdef USEOMP
!$omp do
#endif
      do iq = firstofset( mpiglobal%rank, wfint_kset%nkpt), lastofset( mpiglobal%rank, wfint_kset%nkpt)
        call wfint_interpolate_dmat( lmax, lammax, iq, radcoeffr, radolp, dmat, diagonly=.false.)
        do ias = 1, natmtot
          do ist = 1, wf_nwf
            call zgemm( 'n', 'n', lmmax, lmmax, lmmax, zone, &
                   ulm( :, :, ias), lmmax, &
                   dmat( :, :, ist, ias), lmmax, zzero, &
                   auxmat, lmmax)
            call zgemm( 'n', 'c', lmmax, lmmax, lmmax, zone, &
                   auxmat, lmmax, &
                   ulm( :, :, ias), lmmax, zzero, &
                   dmat( :, :, ist, ias), lmmax)
            do l = 0, lmax
              do m = -l, l
                lm = idxlm( l, m)
                bc( l, ias, ist, iq) = bc( l, ias, ist, iq) + dble( dmat( lm, lm, ist, ias))
              end do
            end do
          end do
        end do
      end do
#ifdef USEOMP
!$omp end do
#endif
      deallocate( dmat, auxmat)
#ifdef USEOMP
!$omp end parallel
#endif
      call mpi_allgatherv_ifc( set=wfint_kset%nkpt, rlen=(lmax+1)*natmtot*wf_nwf, rbuf=bc)
      call barrier

      if( allocated( radcoeffr)) deallocate( radcoeffr)
      if( allocated( radolp)) deallocate( radolp)

      return
      
    end subroutine wfint_interpolate_bandchar
    !EOC

!--------------------------------------------------------------------------------------

    !BOP
    ! !ROUTINE: wfint_interpolate_dos
    ! !INTERFACE:
    !
    subroutine wfint_interpolate_dos( lmax, nsmooth, intgrid, neffk, nsube, ewin, tdos, scissor, inttype, pdos, lonly, jdos, ntrans, mtrans)
      ! !USES
      use mod_eigenvalue_occupancy, only: occmax
      use mod_opt_tetra
      ! !INPUT PARAMETERS:
      !   lmax    : maximum angular quantum number (in, integer)
      !   nsmooth : smoothing of output result (in, integer)
      !   intgrid : interpolation/integration grid (in, integer(3))
      !   neffk   : effective number of k-points in one direction (in, integer)
      !   nsube   : number of energy subdivisions in energy window (in, integer)
      !   ewin    : energy window (in, real(2))
      !   tdos    : total DOS (out, real( nsube))
      !   scissor : energy scissor to apply (in, optional, real)
      !   inttype : integration method (in, optional, character(64))
      !   pdos    : partial DOS (out, optional, real( nsube, (lmax+1)**2, natmtot))
      !   lonly   : partial DOS only l-resolved (in, optional, logical)
      !   jdos    : joint DOS (out, optional, real( nsube, 0:wf_nwf))
      !   ntrans  : number of transistion in JDOS (out, optional, integer)
      !   mtrans  : number of transistion in JDOS (out, optional, integer)
      ! !DESCRIPTION:
      ! Calculates the total density of states (DOS) and optionally also the partial or joint DOS.
      !
      ! !REVISION HISTORY:
      !   Created July 2017 (SeTi)
      !EOP
      integer, intent( in) :: lmax, nsmooth, intgrid(3), neffk, nsube
      real(8), intent( in) :: ewin(2)
      real(8), intent( out) :: tdos( nsube)
      ! optional arguments
      real(8), optional, intent( in)       :: scissor
      logical, optional, intent( in)       :: lonly
      character(64), optional, intent( in) :: inttype
      real(8), optional, intent( out)      :: pdos( nsube, (lmax+1)**2, natmtot)
      real(8), optional, intent( out)      :: jdos( nsube, 0:wf_nwf)
      integer, optional, intent( out)      :: ntrans, mtrans
      !BOC

      logical :: genpdos, genjdos, pdoslonly
      integer :: lmmax, ias, l, m, lm, ist, jst, iq, q1, q2, ie, nk(3), lammax, n
      real(8) :: dosscissor, tmpfermi
      character(64) :: integraltype
      type( k_set) :: tmp_kset
      type( t_set) :: tset

      real(8), allocatable :: energies(:,:), radolp(:,:,:,:), elm(:,:), e(:), ftdos(:,:,:), fjdos(:,:), fpdos(:,:,:,:), edif(:,:,:), ejdos(:,:), ijdos(:,:,:), wjdos(:,:,:)
      complex(8), allocatable :: ulm(:,:,:), radcoeffr(:,:,:,:,:), dmat(:,:,:,:), auxmat(:,:)

      dosscissor = 0.d0
      if( present( scissor)) dosscissor = scissor

      genpdos = .false.
      if( present( pdos)) genpdos = .true.

      genjdos = .false.
      if( present( jdos)) genjdos = .true.

      integraltype = 'tetra'
      if( present( inttype)) integraltype = inttype

      pdoslonly = .true.
      if( present( lonly)) pdoslonly = lonly

      lmmax = (lmax+1)**2

      call wfint_findbandmap

      ! we don't want to use the libzint routine in case of Hybrids (too slow for dense grids)
      ! in case of Hybrids, stypenumber < 0 --> libzint k-point generation is invoked inside generate_k_vectors
      ! temporarilly set stypenumber = 1 to invoke exciting internal k-point generation
      l = input%groundstate%stypenumber
      input%groundstate%stypenumber = 1
      call generate_k_vectors( tmp_kset, bvec, intgrid, wf_kset%vkloff, .true., uselibzint=.false.)
      input%groundstate%stypenumber = l
      nk(:) = max( neffk/tmp_kset%ngridk, 1)

      allocate( e( nsube))
      do ie = 1, nsube
        e( ie) = ewin(1) + dble( ie-1)*(ewin(2)-ewin(1))/(nsube-1)
      end do

      call wfint_init( tmp_kset)

      call wfint_interpolate_occupancy( usetetra=(integraltype == 'tetra'))
      
      allocate( energies( wf_nwf, wfint_kset%nkpt))
      energies = wfint_eval
      tmpfermi = wfint_efermi
      if( wf_fermizero) tmpfermi = 0.d0
      if( wf_fermizero) energies = energies - wfint_efermi

      ! apply scissor
      do iq = 1, wfint_kset%nkpt
        do ist = 1, wf_nwf
          if( energies( ist, iq) .gt. tmpfermi) energies( ist, iq) = energies( ist, iq) + dosscissor
        end do
      end do

      !--------------------------------------!
      !              total DOS               !
      !--------------------------------------!
      allocate( ftdos( wf_nwf, wfint_kset%nkpt, nsube))
      ftdos = 1.d0
      if( integraltype == 'trilin') then
        call brzint( nsmooth, wfint_kset%ngridk, nk, wfint_kset%ikmap, nsube, ewin, wf_nwf, wf_nwf, &
             energies, &
             ftdos(:,:,1), &
             tdos)
      else if( integraltype == 'trilin+') then
        call brzint_new( nsmooth, wfint_kset%ngridk, nk, wfint_kset%ikmap, nsube, ewin, wf_nwf, wf_nwf, &
             energies, &
             ftdos(:,:,1), &
             tdos)
      else
        call opt_tetra_wgt_delta( wfint_tetra, wfint_kset%nkpt, wf_nwf, energies, nsube, e, ftdos)
        tdos = 0.d0
#ifdef USEOMP
!$omp parallel default( shared) private( ie)
!$omp do
#endif
        do ie = 1, nsube
          tdos( ie) = sum( ftdos( :, :, ie))
        end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
      end if

      !--------------------------------------!
      !             partial DOS              !
      !--------------------------------------!
      if( genpdos) then
        allocate( elm( lmmax, natmtot))
        allocate( ulm( lmmax, lmmax, natmtot))
        call genlmirep( lmax, lmmax, elm, ulm)
        deallocate( elm)
        allocate( fpdos( wf_nwf, lmmax, natmtot, wfint_kset%nkpt))
        fpdos(:,:,:,:) = 0.d0
        call wfint_gen_radcoeffr( lmax, lammax, radcoeffr)
        call wfint_gen_radolp( lmax, lammax, radolp)

        q1 = firstofset( mpiglobal%rank, wfint_kset%nkpt)
        q2 = lastofset( mpiglobal%rank, wfint_kset%nkpt)
        allocate( dmat( lmmax, lmmax, wf_nwf, natmtot))
        allocate( auxmat( lmmax, lmmax))
        do iq = q1, q2
          call wfint_interpolate_dmat( lmax, lammax, iq, radcoeffr, radolp, dmat)
          do ias = 1, natmtot
            do ist = 1, wf_nwf
              call zgemm( 'n', 'n', lmmax, lmmax, lmmax, zone, &
                     ulm( :, :, ias), lmmax, &
                     dmat( :, :, ist, ias), lmmax, zzero, &
                     auxmat, lmmax)
              call zgemm( 'n', 'c', lmmax, lmmax, lmmax, zone, &
                     auxmat, lmmax, &
                     ulm( :, :, ias), lmmax, zzero, &
                     dmat( :, :, ist, ias), lmmax)
              do l = 0, lmax
                do m = -l, l
                  lm = idxlm( l, m)
                  if( pdoslonly) then
                    fpdos( ist, l+1, ias, iq) = fpdos( ist, l+1, ias, iq) + dble( dmat( lm, lm, ist, ias))
                  else
                    fpdos( ist, lm, ias, iq) = dble( dmat( lm, lm, ist, ias))
                  end if
                end do
              end do
            end do
          end do
        end do
        deallocate( dmat, auxmat, ulm)
        call mpi_allgatherv_ifc( set=wfint_kset%nkpt, rlen=lmmax*natmtot*wf_nwf, rbuf=fpdos)
        call barrier

        do ias = 1, natmtot
          do l = 0, lmax
            if( pdoslonly) then
              if( integraltype == 'trilin') then
                call brzint( nsmooth, wfint_kset%ngridk, nk, wfint_kset%ikmap, nsube, ewin, wf_nwf, wf_nwf, &
                       energies, &
                       fpdos( :, l+1, ias, :), &
                       pdos( :, l+1, ias))
              else if( integraltype == 'trilin+') then
                call brzint_new( nsmooth, wfint_kset%ngridk, nk, wfint_kset%ikmap, nsube, ewin, wf_nwf, wf_nwf, &
                       energies, &
                       fpdos( :, l+1, ias, :), &
                       pdos( :, l+1, ias))
              else
#ifdef USEOMP
!$omp parallel default( shared) private( ie)
!$omp do
#endif
                do ie = 1, nsube
                  pdos( ie, l+1, ias) = sum( fpdos( :, l+1, ias, :)*ftdos( :, :, ie))
                end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
              end if
            else
              do m = -l, l
                lm = idxlm( l, m)
                if( integraltype == 'trilin') then
                  call brzint( nsmooth, wfint_kset%ngridk, nk, wfint_kset%ikmap, nsube, ewin, wf_nwf, wf_nwf, &
                         energies, &
                         fpdos( :, lm, ias, :), &
                         pdos( :, lm, ias))
                else if( integraltype == 'trilin+') then
                  call brzint_new( nsmooth, wfint_kset%ngridk, nk, wfint_kset%ikmap, nsube, ewin, wf_nwf, wf_nwf, &
                         energies, &
                         fpdos( :, lm, ias, :), &
                         pdos( :, lm, ias))
                else
#ifdef USEOMP
!$omp parallel default( shared) private( ie)
!$omp do
#endif
                  do ie = 1, nsube
                    pdos( ie, lm, ias) = sum( fpdos( :, lm, ias, :)*ftdos( :, :, ie))
                  end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
                end if
              end do
            end if
          end do
        end do

        deallocate( fpdos)
      end if
      deallocate( ftdos)

      !--------------------------------------!
      !              joint DOS               !
      !--------------------------------------!
      if( genjdos) then
        allocate( edif( wf_nwf, wfint_kset%nkpt, wf_nwf))
        edif = -1.d100
        mtrans = 0
        ntrans = 0
        do iq = 1, wfint_kset%nkpt
          n = 0
          do ist = wf_nwf, 1, -1
            if( energies( ist, iq) .le. tmpfermi) then
              n = n + 1
              m = 0
              do jst = 1, wf_nwf
                if( energies( jst, iq) .gt. tmpfermi) then
                  m = m + 1
                  edif( m, iq, n) = energies( jst, iq) - energies( ist, iq)
                end if
              end do
              mtrans = max( mtrans, m)
            end if
          end do
          ntrans = max( ntrans, n)
        end do
        ! state dependent JDOS
        if( integraltype == 'trilin' .or. integraltype == 'trilin+') then
          allocate( fjdos( mtrans*ntrans, wfint_kset%nkpt))
          allocate( ejdos( mtrans*ntrans, wfint_kset%nkpt))
          fjdos(:,:) = 1.d0
          ejdos(:,:) = -1.d100
#ifdef USEOMP
!$omp parallel default( shared) private( ist)
!$omp do
#endif
          do ist = 1, ntrans
            if( integraltype == 'trilin') then
              call brzint( nsmooth, wfint_kset%ngridk, nk, wfint_kset%ikmap, nsube, ewin, mtrans, mtrans, &
                     edif( 1:mtrans, :, ist), &
                     fjdos( 1:mtrans, :), &
                     jdos( :, ist))
            else
              call brzint_new( nsmooth, wfint_kset%ngridk, nk, wfint_kset%ikmap, nsube, ewin, mtrans, mtrans, &
                     edif( 1:mtrans, :, ist), &
                     fjdos( 1:mtrans, :), &
                     jdos( :, ist))
            end if
          end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
          ! total JDOS
          iq = 0
          do ist = 1, ntrans
            do jst = 1, mtrans
              if( (.not. (maxval( edif( jst, :, ist)) .lt. ewin(1))) .and. (.not. (minval( edif( jst, :, ist)) .gt. ewin(2)))) then
                iq = iq + 1
                ejdos( iq, :) = edif( jst, :, ist)
              end if
            end do
          end do
  
          if( integraltype == 'trilin') then
            call brzint( nsmooth, wfint_kset%ngridk, nk, wfint_kset%ikmap, nsube, ewin, iq, iq, &
                   ejdos( 1:iq, :), &
                   fjdos( 1:iq, :), &
                   jdos( :, 0))
          else
            call brzint_new( nsmooth, wfint_kset%ngridk, nk, wfint_kset%ikmap, nsube, ewin, iq, iq, &
                   ejdos( 1:iq, :), &
                   fjdos( 1:iq, :), &
                   jdos( :, 0))
          end if
          jdos = jdos*occmax
          deallocate( fjdos, ejdos)            
        else
          do l = 1, wf_nwf      ! lowest (partially) unoccupied band
            if( occmax-minval( wfint_occ(l,:)) .gt. 1.d-10) exit
          end do
          do m = wf_nwf, 1, -1  ! highest (partially) occupied band
            if( maxval( wfint_occ(m,:)) .gt. 1.d-10) exit
          end do
            
          allocate( ijdos( nsube, l:wf_nwf, m))
          allocate( wjdos( l:wf_nwf, m, wfint_kset%nkpt))
          do iq = 1, wfint_kset%nkpt
            do jst = 1, m
              do ist = l, wf_nwf
                wjdos( ist, jst, iq) = min( wfint_occ( jst, iq), occmax-wfint_occ( ist, iq))
              end do
            end do
          end do

          call opt_tetra_init( tset, wfint_kset, 2, reduce=.true.)
          call opt_tetra_int_deltadiff( tset, wfint_kset%nkpt, wf_nwf-l+1, energies( l:wf_nwf, :), m, energies( 1:m, :), &
                 nsube, e, 1, resr=ijdos, matr=wjdos)

          jdos = 0.d0
          do jst = 1, m
            jdos(:,jst) = sum( ijdos(:,:,jst), 2)
          end do
          jdos(:,0) = sum( jdos(:,1:m), 2)
                 
          deallocate( ijdos, wjdos)
        end if
        deallocate( edif)
      end if

      deallocate( e, energies)
      if( allocated( radcoeffr)) deallocate( radcoeffr)
      if( allocated( radolp)) deallocate( radolp)

      return
      
    end subroutine wfint_interpolate_dos
    !EOC

!--------------------------------------------------------------------------------------
      
    ! interpolates the density matrix
    ! needed for band character and PDOS
    subroutine wfint_interpolate_dmat( lmax, lammax, iq, radcoeffr, radolp, dmat, diagonly)
      integer, intent( in)           :: lmax, lammax, iq
      complex(8), intent( in)        :: radcoeffr( wf_nwf, lammax, (lmax+1)**2, natmtot, wf_nrpt)
      real(8), intent( in)           :: radolp( lammax, lammax, 0:lmax, natmtot)
      complex(8), intent( out)       :: dmat( (lmax+1)**2, (lmax+1)**2, wf_nwf, natmtot)
      logical, optional, intent( in) :: diagonly

      integer :: ir, is, ia, ias, o, l1, m1, lm1, m2, lm2, lmmax, ilo1, maxdim, ist
      integer :: lamcnt( 0:lmax, nspecies), o2idx( apwordmax, 0:lmax, nspecies), lo2idx( nlomax, 0:lmax, nspecies), lm2l( (lmax+1)**2)

      complex(8), allocatable :: radcoeffq1(:,:), radcoeffq2(:,:), U(:,:), auxmat(:,:), wgts(:,:,:)
      logical :: diag

      complex(8) :: zdotc

      diag = .false.
      if( present( diagonly)) diag = diagonly

      lm2l = 0
      do l1 = 0, lmax
        do m1 = -l1, l1
          lm2l( idxlm( l1, m1)) = l1
        end do
      end do
      lmmax = (lmax+1)**2
      maxdim = 0
      do is = 1, nspecies
        o2idx( :, :, is) = 0
        lo2idx( :, :, is) = 0
        lamcnt( :, is) = 0
        do l1 = 0, lmax
          do o = 1, apword( l1, is)
            lamcnt( l1, is) = lamcnt( l1, is) + 1
            o2idx( o, l1, is) = lamcnt( l1, is)
          end do
        end do
        do ilo1 = 1, nlorb( is)
          l1 = lorbl( ilo1, is)
          if( l1 .le. lmax) then
            lamcnt( l1, is) = lamcnt( l1, is) + 1
            lo2idx( ilo1, l1, is) = lamcnt( l1, is)
          end if
        end do
        maxdim = max( maxdim, maxval( lamcnt( :, is)))
      end do
      if( maxdim .ne. lammax) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Error (wfint_interpolate_dmat): Inconsistent input. Check lmax and lammax.")')
        end if
        stop
      end if

      ! build q-point density coefficients
      if( wfint_mindist) then
        allocate( wgts( wf_nwf, wf_nwf, wf_nrpt))
        do ir = 1, wf_nrpt
          call wfint_getmindistwgts( ir, iq, wgts( :, :, ir))
          wgts( :, :, ir) = transpose( wgts( :, :, ir))
        end do
      end if
          
      dmat = zzero

#ifdef USEOMP
!$omp parallel default( shared) private( is, ia, ias, l1, lm1, m2, lm2, radcoeffq1, radcoeffq2, ir, ist, U, auxmat)
#endif
      allocate( radcoeffq1( maxdim, wf_nwf))
      allocate( radcoeffq2( maxdim, wf_nwf))
      allocate( auxmat( maxdim, wf_nwf))
      allocate( U( wf_nwf, wf_nwf))
      U = wfint_transform(:,:,iq)
#ifdef USEOMP
!$omp do
#endif
      do lm1 = 1, lmmax
        l1 = lm2l( lm1)

        do is = 1, nspecies
          do ia = 1, natoms( is)

            ias = idxas( ia, is)
            radcoeffq1 = zzero
            do ir = 1, wf_nrpt
              if( wfint_mindist) U = wfint_transform(:,:,iq)*wgts(:,:,ir)
              call zgemm( 't', 'n', maxdim, wf_nwf, wf_nwf, wfint_pqr( iq, ir), &
                     radcoeffr( :, :, lm1, ias, ir), wf_nwf, &
                     U, wf_nwf, zone, &
                     radcoeffq1, maxdim)
            end do
            if( diag) then
              call zgemm( 't', 'n', maxdim, wf_nwf, maxdim, zone, &
                     cmplx( radolp( :, :, l1, ias), 0, 8), maxdim, &
                     radcoeffq1, maxdim, zzero, &
                     auxmat, maxdim)
              do ist = 1, wf_nwf
                dmat( lm1, lm1, ist, ias) = zdotc( lamcnt( l1, is), radcoeffq1( :, ist), 1, auxmat( :, ist), 1) 
              end do
            else
              do m2 = -l1, l1
                lm2 = idxlm( l1, m2)
                ! calculate q-point density coefficient
                radcoeffq2 = zzero
                do ir = 1, wf_nrpt
                  if( wfint_mindist) U = wfint_transform(:,:,iq)*wgts(:,:,ir)
                  call zgemm( 't', 'n', maxdim, wf_nwf, wf_nwf, wfint_pqr( iq, ir), &
                         radcoeffr( :, :, lm2, ias, ir), wf_nwf, &
                         U, wf_nwf, zone, &
                         radcoeffq2, maxdim)
                end do
                call zgemm( 't', 'n', maxdim, wf_nwf, maxdim, zone, &
                       cmplx( radolp( :, :, l1, ias), 0, 8), maxdim, &
                       radcoeffq1, maxdim, zzero, &
                       auxmat, maxdim)
                do ist = 1, wf_nwf
                  dmat( lm1, lm2, ist, ias) = zdotc( lamcnt( l1, is), radcoeffq2( :, ist), 1, auxmat( :, ist), 1) 
                end do
              end do
            end if

          end do
        end do

      end do
#ifdef USEOMP
!$omp end do
#endif
      deallocate( radcoeffq1, radcoeffq2, auxmat, U)
#ifdef USEOMP
!$omp end parallel
#endif
      if( wfint_mindist) deallocate( wgts)
      return
      
    end subroutine wfint_interpolate_dmat

!--------------------------------------------------------------------------------------
      
    ! generates radial overlaps of basis functions
    subroutine wfint_gen_radolp( lmax, lammax, radolp)
      integer, intent( in) :: lmax
      integer, intent( out) :: lammax
      real(8), allocatable, intent( out) :: radolp(:,:,:,:)

      integer :: is, ia, ias, o, l1, lmmax, ilo1, ilo2, maxdim
      integer :: lamcnt( 0:lmax, nspecies), o2idx( apwordmax, 0:lmax, nspecies), lo2idx( nlomax, 0:lmax, nspecies)

      lmmax = (lmax + 1)**2

      call wfhelp_genradfun

      maxdim = 0
      do is = 1, nspecies
        o2idx( :, :, is) = 0
        lo2idx( :, :, is) = 0
        lamcnt( :, is) = 0
        do l1 = 0, lmax
          do o = 1, apword( l1, is)
            lamcnt( l1, is) = lamcnt( l1, is) + 1
            o2idx( o, l1, is) = lamcnt( l1, is)
          end do
        end do
        do ilo1 = 1, nlorb( is)
          l1 = lorbl( ilo1, is)
          if( l1 .le. lmax) then
            lamcnt( l1, is) = lamcnt( l1, is) + 1
            lo2idx( ilo1, l1, is) = lamcnt( l1, is)
          end if
        end do
        maxdim = max( maxdim, maxval( lamcnt( :, is)))
      end do
      lammax = maxdim

      ! radial overlap-intergrals
      if( allocated( radolp)) deallocate( radolp)
      allocate( radolp( maxdim, maxdim, 0:lmax, natmtot))
      radolp(:,:,:,:) = 0.d0
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          ! APW-APW integral
          do l1 = 0, lmax
            do o = 1, apword( l1, is)
              radolp( o2idx( o, l1, is), o2idx( o, l1, is), l1, ias) = 1.d0
            end do
          end do
          do ilo1 = 1, nlorb( is)
            l1 = lorbl( ilo1, is)
            if( l1 .le. lmax) then
              ! APW-LO integral
              do o = 1, apword( l1, is)
                if( (o2idx( o, l1, is) .gt. 0) .and. (lo2idx( ilo1, l1, is) .gt. 0)) then
                  radolp( o2idx( o, l1, is), lo2idx( ilo1, l1, is), l1, ias) = oalo( o, ilo1, ias)
                  radolp( lo2idx( ilo1, l1, is), o2idx( o, l1, is), l1, ias) = oalo( o, ilo1, ias)
                end if
              end do
              ! LO-LO integral
              do ilo2 = 1, nlorb( is)
                if( lorbl( ilo2, is) .eq. l1) then
                  if( (lo2idx( ilo1, l1, is) .gt. 0) .and. (lo2idx( ilo2, l1, is) .gt. 0)) then
                    radolp( lo2idx( ilo1, l1, is), lo2idx( ilo2, l1, is), l1, ias) = ololo( ilo1, ilo2, ias)
                    radolp( lo2idx( ilo2, l1, is), lo2idx( ilo1, l1, is), l1, ias) = ololo( ilo2, ilo1, ias)
                  end if
                end if
              end do
            end if
          end do
        end do
      end do
   
      return
      
    end subroutine wfint_gen_radolp

!--------------------------------------------------------------------------------------
      
    ! generates R-dependent muffin-tin density coefficients
    subroutine wfint_gen_radcoeffr( lmax, lammax, radcoeffr)
      integer, intent( in) :: lmax
      integer, intent( out) :: lammax
      complex(8), allocatable, intent( out) :: radcoeffr(:,:,:,:,:)

      integer :: ik, ir, is, ia, ias, l1, m1, lm1, o, ilo1, lmmax, ngknr, maxdim
      integer :: lamcnt( 0:lmax, nspecies), o2idx( apwordmax, 0:lmax, nspecies), lo2idx( nlomax, 0:lmax, nspecies)

      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), radcoeffk(:,:,:,:,:)

      lmmax = (lmax + 1)**2

      call wfhelp_genradfun

      maxdim = 0
      do is = 1, nspecies
        o2idx( :, :, is) = 0
        lo2idx( :, :, is) = 0
        lamcnt( :, is) = 0
        do l1 = 0, lmax
          do o = 1, apword( l1, is)
            lamcnt( l1, is) = lamcnt( l1, is) + 1
            o2idx( o, l1, is) = lamcnt( l1, is)
          end do
        end do
        do ilo1 = 1, nlorb( is)
          l1 = lorbl( ilo1, is)
          if( l1 .le. lmax) then
            lamcnt( l1, is) = lamcnt( l1, is) + 1
            lo2idx( ilo1, l1, is) = lamcnt( l1, is)
          end if 
        end do
        maxdim = max( maxdim, maxval( lamcnt( :, is)))
      end do
      lammax = maxdim
   
      ! build k-point density coefficients
      allocate( radcoeffk( wf_fst:wf_lst, maxdim, lmmax, natmtot, wf_kset%nkpt))
      allocate( evecfv( nmatmax_ptr, nstfv, nspnfv))
      allocate( apwalm( ngkmax_ptr, apwordmax, lmmaxapw, natmtot, nspnfv))
      radcoeffk(:,:,:,:,:) = zzero

      do ik = 1, wf_kset%nkpt
        ngknr = wf_Gkset%ngk( 1, ik)

        ! get matching coefficients
        call match( ngknr, wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm( :, :, :, :, 1))
          
        ! read eigenvector      
        call wfhelp_getevec( ik, evecfv)

        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            ! APW contribution
            do l1 = 0, lmax
              do m1 = -l1, l1
                lm1 = idxlm( l1, m1)
                do o = 1, apword( l1, is)
                  call zgemv( 't', ngknr, wf_nwf, zone, &
                         evecfv( 1, wf_fst, 1), nmatmax_ptr, &
                         apwalm( 1, o, lm1, ias, 1), 1, zzero, &
                         radcoeffk( wf_fst, o2idx( o, l1, is), lm1, ias, ik), 1)
                end do
              end do
            end do
            ! LO contribution
            do ilo1 = 1, nlorb( is)
              l1 = lorbl( ilo1, is)
              if( l1 .le. lmax) then
                do m1 = -l1, l1
                  lm1 = idxlm( l1, m1)
                  radcoeffk( :, lo2idx( ilo1, l1, is), lm1, ias, ik) = evecfv( ngknr+idxlo( lm1, ilo1, ias), wf_fst:wf_lst, 1)
                end do
              end if
            end do
          end do
        end do
      end do

      deallocate( evecfv, apwalm)

      ! build R-point density coefficients
      if( allocated( radcoeffr)) deallocate( radcoeffr)
      allocate( radcoeffr( wf_nwf, maxdim, lmmax, natmtot, wf_nrpt))
      radcoeffr(:,:,:,:,:) = zzero
#ifdef USEOMP                
!$omp parallel default( shared) private( ir, ik, is, ia, ias, lm1)
!$omp do
#endif
      do ir = 1, wf_nrpt
        do ik = 1, wf_kset%nkpt
          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)
              do lm1 = 1, lmmax
                call zgemm( 't', 'n', wf_nwf, maxdim, wf_nst, wf_pkr( ik, ir)/wf_kset%nkpt, &
                       wf_transform( :, :, ik), wf_nst, &
                       radcoeffk( :, :, lm1, ias, ik), wf_nst, zone, &
                       radcoeffr( :, :, lm1, ias, ir), wf_nwf)
              end do
            end do
          end do
        end do
      end do
#ifdef USEOMP                
!$omp end do
!$omp end parallel
#endif
      
      deallocate( radcoeffk)

      return
      
    end subroutine wfint_gen_radcoeffr

!--------------------------------------------------------------------------------------

    subroutine wfint_findbandmap
      integer :: ik, ist, jst, fst, lst
      
      integer, allocatable :: map(:,:)
      real(8), allocatable :: eval(:,:)

      if( allocated( wfint_bandmap)) deallocate( wfint_bandmap)
      allocate( wfint_bandmap( wf_nwf))
      wfint_bandmap = 0

      if( input%properties%wannier%input .eq. 'gw') then
        call wfhelp_geteval( eval, fst, lst, mode='gwks')
      else
        call wfhelp_geteval( eval, fst, lst)
      end if
      call wfint_init( wf_kset, evalin=eval( wf_fst:wf_lst, :))

      allocate( map( wf_nwf, wf_kset%nkpt))
      map = 0
      do ik = 1, wf_kset%nkpt
        do ist = 1, wf_nwf
          do jst = fst, lst
            if( abs( eval( jst, ik) - wfint_eval( ist, ik)) .lt. 1.d-6) then
              map( ist, ik) = jst
              exit
            end if
          end do
        end do
      end do
      fst = 0
      do ist = 1, wf_nwf
        jst = minval( map( ist, :))
        if( (jst .eq. maxval( map( ist, :))) .and. (jst .ne. 0)) then
          wfint_bandmap( ist) = jst
          if( fst .eq. 0) fst = ist
        end if
      end do
      if( fst .eq. 0) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Error (wfint_findbandmap) Bands cannot be maped!")')
        end if
        stop
      end if

      if( fst .gt. 1) then
        do ist = fst - 1, 1, -1
          wfint_bandmap( ist) = wfint_bandmap( ist + 1) - 1
        end do
      end if
      do ist = fst + 1, wf_nwf
        if( wfint_bandmap( ist) .eq. 0) wfint_bandmap( ist) = wfint_bandmap( ist - 1) + 1
      end do

      deallocate( eval, map)
      return
    end subroutine wfint_findbandmap
    
!--------------------------------------------------------------------------------------

    !BOP
    ! !ROUTINE: wfint_matchksgw_linreal
    ! !INTERFACE:
    !
    subroutine wfint_matchksgw_linreal( kset, fst, lst, rin, rout, fill)
      ! !INPUT PARAMETERS:
      !   kset     : kset on which the quantity is given (in, type( k_set))
      !   fst, lst : first and last band index for which the input is provided (in, integer)
      !   rin      : real scalar quantity corresponding to KS bands (in, real( fst:lst, kset%nkpt))
      !   rout     : quantity matched onto GW Wannier bands (out, real( wf_nwf, kset%nkpt))
      !   fill     : fill value for bands outside the wannierized bands (default: 0.d0) (in, optional, real)
      ! !DESCRIPTION:
      ! Matches a scalar real quatity corresponding to the KS bands at the given kset
      ! onto the Wannier interpolated GW bands at the given kset.
      !
      ! !REVISION HISTORY:
      !   Created June 2019 (SeTi)
      !EOP
      type( k_set), intent( in)      :: kset
      integer, intent( in)           :: fst, lst
      real(8), intent( in)           :: rin( fst:lst, kset%nkpt)
      real(8), intent( out)          :: rout( wf_nwf, kset%nkpt)
      real(8), optional, intent( in) :: fill
      !BOC

      integer :: iq, ist, jst
      real(8) :: f
      real(8), allocatable :: r(:), evalinks(:,:), evalinqp(:,:)
      complex(8), allocatable :: olp(:,:), evecintks(:,:,:)

      f = 0.d0
      if( present( fill)) f = fill

      call wfint_findbandmap

      allocate( olp( wf_nwf, wf_nwf), evecintks( wf_nwf, wf_nwf, kset%nkpt), r( wf_nwf))
      call wfhelp_geteval( evalinks, ist, jst, mode='gwks')
      call wfhelp_geteval( evalinqp, ist, jst, mode='gw')
      call wfint_init( kset, evalin=evalinks( wf_fst:wf_lst, :))
      evecintks = wfint_transform
      call wfint_init( kset, evalin=evalinqp( wf_fst:wf_lst, :))

      do iq = 1, kset%nkpt
        r = f
        do ist = 1, wf_nwf
          jst = wfint_bandmap( ist)
          if( (jst .ge. fst) .and. (jst .le. lst)) r( ist) = rin( jst, iq)
        end do
        call zgemm( 'c', 'n', wf_nwf, wf_nwf, wf_nwf, zone, &
               evecintks( :, :, iq), wf_nwf, &
               wfint_transform( :, :, iq), wf_nwf, zzero, &
               olp, wf_nwf)
        call dgemv( 't', wf_nwf, wf_nwf, 1.d0, &
               abs( olp)**2, wf_nwf, &
               r, 1, 0.d0, &
               rout( :, iq), 1)
      end do

      deallocate( r, olp, evecintks)
      if( allocated( evalinks)) deallocate( evalinks)
      if( allocated( evalinqp)) deallocate( evalinqp)
      
    end subroutine wfint_matchksgw_linreal
    !EOC

end module mod_wannier_interpolate
