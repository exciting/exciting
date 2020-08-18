module mod_wannier_bse
  use mod_wannier
  use mod_wannier_interpolate
  use m_linalg
  implicit none

! module variables
  type( k_set)  :: wfbse_kset
  type( G_set)  :: wfbse_Gset
  type( Gk_set) :: wfbse_Gkset
  integer       :: wfbse_ofst, wfbse_olst, wfbse_ufst, wfbse_ulst
  integer       :: wfbse_olstwf, wfbse_ufstwf
  logical       :: wfbse_initialized
  logical       :: wfbse_ordered
  real(8), allocatable :: wfbse_eval(:,:)

! methods
  contains

    function wfbse_usegwwannier() result( b)
      logical :: b

      b = .false.
      if( associated( input%properties)) then
        if( associated( input%properties%wannier)) then
          if( wf_initialized .and. input%properties%wannier%input .eq. 'gw') b = .true.
        end if
      end if
      return
    end function wfbse_usegwwannier

    subroutine wfbse_init
      integer :: ist, jst, ik, nocc
      logical :: check
      real(8) :: olst, ufst
      real(8), allocatable :: eval(:,:)

      if( wfbse_initialized) return

      call generate_k_vectors( wfbse_kset, wf_kset%bvec, input%xs%ngridk, input%xs%vkloff, input%xs%reducek, .false.)
      call generate_G_vectors( wfbse_Gset, wfbse_kset%bvec, intgv, input%groundstate%gmaxvr)
      call generate_Gk_vectors( wfbse_Gkset, wfbse_kset, wfbse_Gset, gkmax)

      ! check band ranges
      nocc = nint( chgval/occmax)
      wfbse_ofst = input%xs%bse%nstlbse(1)
      wfbse_olst = input%xs%bse%nstlbse(2)
      wfbse_ufst = input%xs%bse%nstlbse(3) + nocc
      wfbse_ulst = input%xs%bse%nstlbse(4) + nocc
      do ist = wfbse_ofst, wfbse_olst
        check = .false.
        do wf_group = 1, wf_ngroups
          if( wf_groups( wf_group)%nst .eq. wf_groups( wf_group)%nwf) then
            if( (ist .lt. wf_groups( wf_group)%fst) .or. (ist .gt. wf_groups( wf_group)%lst)) then
              check = .true.
            else
              check = .false.
              exit
            end if
          else
            do ik = 1, wf_kset%nkpt
              if( .not. any( wf_groups( wf_group)%win_ii( :, ik) .eq. ist)) then
                check = .true.
                exit
              else
                check = .false.
              end if
            end do
            if( .not. check) exit
          end if
        end do
        if( check) then
          write(*,*)
          write( *, '("Error (wfbse_init): The occupied bands considered in the BSE must be fully contained either in an isolated group or in the inner window of an entangled group in the Wannier function calculation.")')
          call terminate
        end if
      end do
      do ist = wfbse_ufst, wfbse_ulst
        check = .false.
        do wf_group = 1, wf_ngroups
          if( wf_groups( wf_group)%nst .eq. wf_groups( wf_group)%nwf) then
            if( (ist .lt. wf_groups( wf_group)%fst) .or. (ist .gt. wf_groups( wf_group)%lst)) then
              check = .true.
            else
              check = .false.
              exit
            end if
          else
            do ik = 1, wf_kset%nkpt
              if( .not. any( wf_groups( wf_group)%win_ii( :, ik) .eq. ist)) then
                check = .true.
                exit
              else
                check = .false.
              end if
            end do
            if( .not. check) exit
          end if
        end do
        if( check) then
          write(*,*)
          write( *, '("Error (wfbse_init): The unoccupied bands considered in the BSE must be fully contained either in an isolated group or in the inner window of an entangled group in the Wannier function calculation.")')
          call terminate
        end if
      end do

      ! map Wannier KS bands to original KS bands
      call wfhelp_geteval( eval, ist, jst, mode='gwks')
      call wfint_init( wf_kset, evalin=eval( wf_fst:wf_lst, :))
      olst = 0.d0
      ufst = 0.d0
      do ik = 1, wf_kset%nkpt
        do ist = wf_nwf, 1, -1
          if( abs( wfint_eval( ist, ik) - eval( wfbse_olst, ik)) .lt. 1.d-6) exit
        end do
        do jst = 1, wf_nwf
          if( abs( wfint_eval( jst, ik) - eval( wfbse_ufst, ik)) .lt. 1.d-6) exit
        end do
        olst = olst + dble( ist)
        ufst = ufst + dble( jst)
      end do
      olst = olst/wf_kset%nkpt
      ufst = ufst/wf_kset%nkpt
      wfbse_olstwf = nint( olst)
      wfbse_ufstwf = nint( ufst)

      ! generate bandmap
      call wfint_findbandmap

      wfbse_initialized = .true.
      
      return
    end subroutine wfbse_init

    subroutine wfbse_ordereval
      integer :: iq, ist, jst
      real(8), allocatable :: eval(:,:), evalint(:,:)
      complex(8), allocatable :: evecint(:,:,:), auxmat(:,:)

      if( wfbse_ordered) return

      if( allocated( wfbse_eval)) deallocate( wfbse_eval)
      allocate( wfbse_eval( nstsv, wfbse_kset%nkpt))
      wfbse_eval = 0.d0

      allocate( evecint( wf_nwf, wf_nwf, wfbse_kset%nkpt))
      allocate( evalint( wf_nwf, wfbse_kset%nkpt))
      allocate( auxmat( wf_nwf, wf_nwf))

      ! interpolate KS bands from GW grid onto BSE grid
      call wfhelp_geteval( eval, ist, jst, mode='gwks')
      call wfint_init( wfbse_kset, evalin=eval( wf_fst:wf_lst, :))
      evecint = wfint_transform
      evalint = wfint_eval
      ! interpolate GW bands from GW grid onto BSE grid
      call wfhelp_geteval( eval, ist, jst, mode='gw')
      call wfint_init( wfbse_kset, evalin=eval( wf_fst:wf_lst, :))
      ! compute overlap between both sets of bands
      ! and match GW energies to KS bands on BSE grid
      do iq = 1, wfbse_kset%nkpt
        call zgemm( 'c', 'n', wf_nwf, wf_nwf, wf_nwf, zone, &
               evecint( :, :, iq), wf_nwf, &
               wfint_transform( :, :, iq), wf_nwf, zzero, &
               auxmat, wf_nwf)
        do ist = 0, wfbse_olst - wfbse_ofst
          jst = maxloc( abs( auxmat( :, wfbse_olstwf - ist)), 1)
          wfbse_eval( wfbse_olst - ist, iq) = wfint_eval( jst, iq)
        end do
        do ist = 0, wfbse_ulst - wfbse_ufst
          jst = maxloc( abs( auxmat( :, wfbse_ufstwf + ist)), 1)
          wfbse_eval( wfbse_ufst + ist, iq) = wfint_eval( jst, iq)
        end do
        !write(*,*) iq
        !write(*,'(100f13.6)') wfbse_eval( :, iq)
      end do
      
      wfbse_ordered = .true.
      
      return
    end subroutine wfbse_ordereval
    
    subroutine wfbse_getevec( ik, evec, mode)
      integer, intent( in)                :: ik
      complex(8), intent( out)            :: evec( nmatmax_ptr, nstfv, nspinor)
      character(*), optional, intent( in) :: mode

      integer :: nmatmax_
      character(16) :: mode_

      mode_ = 'gs'
      if( present( mode)) mode_ = trim( mode)
      
      if( mode_ .eq. 'gs') then
        nmatmax_ = nmatmax
        nmatmax = wf_Gkset%ngkmax + nlotot
        call wfhelp_getevec( ik, evec)
        nmatmax = nmatmax_
      else if( mode_ .eq. 'bse') then
        call getevecfv( wfbse_kset%vkl(:, ik), wfbse_Gkset%vgkl( :, :, :, ik), evec)
      end if

      return
    end subroutine wfbse_getevec

    subroutine wfbse_destroy
      if( allocated( wfbse_eval)) deallocate( wfbse_eval)
      call delete_k_vectors( wfbse_kset)
      call delete_G_vectors( wfbse_Gset)
      call delete_Gk_vectors( wfbse_Gkset)
      wfbse_initialized = .false.
      wfbse_ordered = .false.
      return
    end subroutine wfbse_destroy

end module mod_wannier_bse

