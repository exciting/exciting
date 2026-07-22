module mod_wannier_spin
  use mod_wannier_variables

  use mod_spin,                  only: nspinor
  use mod_eigenvalue_occupancy,  only: nstfv, nstsv
  use mod_kpoint
  use mod_wannier_helper
  use m_linalg, only: zhediag

implicit none
  contains

    !> Set wannier global variable `wf_spin_dis` to true and `wf_spin_dis_method` to "collinear". Other methods are not implemented yet.
    subroutine wfspin_set_wf_spin_dis
      wf_spin_dis = .true.
      wf_spin_dis_method = "collinear"
    end subroutine wfspin_set_wf_spin_dis

    !> Double group size for spin disentanglement and allocate sz_eigenvector for transformation to "spin pure" states.
    !> Groups 1...wf_ngroups are spin down, wf_ngroups+1...2*wf_ngroups are spin up in wf_groups.
    subroutine wfspin_spindis_double_groups
      use constants, only : zzero
      type( wannier_group), allocatable :: wf_groups_copy(:)
      integer :: igroup, igroup_up, igroup_down

      wf_groups_copy = wf_groups
      deallocate( wf_groups )
      allocate( wf_groups(2*wf_ngroups) )
      
      do igroup = 1, wf_ngroups
        igroup_up = igroup
        igroup_down = igroup_up + wf_ngroups
        wf_groups(igroup_up) = wf_groups_copy(igroup)
        wf_groups(igroup_down) = wf_groups_copy(igroup)

        ! set KS-ist variables of wf_groups and allocate sz_eigenvectors
        wf_groups(igroup_up)%win_ni_ks = wf_groups_copy(igroup)%win_ni
        wf_groups(igroup_up)%win_no_ks = wf_groups_copy(igroup)%win_no
        wf_groups(igroup_up)%fst_ks = wf_groups_copy(igroup)%fst
        wf_groups(igroup_up)%lst_ks = wf_groups_copy(igroup)%lst
        wf_groups(igroup_up)%nst_ks = wf_groups_copy(igroup)%nst
        
        wf_groups(igroup_down)%win_ni_ks = wf_groups_copy(igroup)%win_ni
        wf_groups(igroup_down)%win_no_ks = wf_groups_copy(igroup)%win_no
        wf_groups(igroup_down)%fst_ks = wf_groups_copy(igroup)%fst
        wf_groups(igroup_down)%lst_ks = wf_groups_copy(igroup)%lst
        wf_groups(igroup_down)%nst_ks = wf_groups_copy(igroup)%nst

        !double wf_ngroups
        wf_ngroups = 2 * wf_ngroups

        ! allocate sz_eigenvectors for spin up and down groups
        allocate( wf_groups(igroup_up)%sz_eigenvector(wf_groups(igroup_up)%fst_ks:wf_groups(igroup_up)%lst_ks, 1:wf_groups(igroup_up)%nst_ks, wf_kset%nkpt), source=zzero )
        allocate( wf_groups(igroup_down)%sz_eigenvector(wf_groups(igroup_down)%fst_ks:wf_groups(igroup_down)%lst_ks, 1:wf_groups(igroup_down)%nst_ks, wf_kset%nkpt), source=zzero )
      end do
      
    end subroutine wfspin_spindis_double_groups

    !> Divide initial wf_groups into spin up and spin down groups employing method chosen by wf_spin_dis_method
    subroutine wfspin_spindis_divide_groups
      use constants, only : zzero
      use sorting, only : sort_index_1d

      type( wannier_group), allocatable :: wf_groups_copy(:)
      complex(8), allocatable :: sz_evec(:, :), spin_z(:, :)
      real(8), allocatable :: sz_eval(:)
      integer, allocatable :: idx(:)
      integer :: igroup, igroup_up, igroup_down, ik, i
      integer :: fst_ik, lst_ik, nst_ik
      integer :: no_down, no_up, ni_down, ni_up, n_down, n_up
      integer, allocatable :: window_mask(:)
     
      ! copy original wf_groups
      wf_groups_copy = wf_groups(1 : wf_ngroups/2)

      do igroup = 1, wf_ngroups/2
        igroup_up = igroup
        igroup_down = igroup_up + wf_ngroups/2

        ! choose spin disentanglement method
        if ( trim( wf_spin_dis_method ) == "collinear") then
          ! associate spin up/down groups with spin up/down states from collinear groundstate calculation
          ! in the gs output states 1...num_of_basis_functions_sv are spin up,
          ! while states num_of_basis_functions_sv+1...nstsv are spin down.
          ! In wannier_readinput, states are sorted according to their energies via wf_index_map,
          ! shuffling this ordering. Spin character will be extracted from wf_index_map entries.

          do ik = 1, wf_kset%nkpt

            ! define range of bands under consideration at ik, including both outer and inner energy window.
            fst_ik = min( minval(wf_groups_copy(igroup)%win_io(1:wf_groups_copy(igroup)%win_no(ik), ik)), &
                          minval(wf_groups_copy(igroup)%win_ii(1:wf_groups_copy(igroup)%win_ni(ik), ik)) )
            lst_ik = max( maxval(wf_groups_copy(igroup)%win_io(1:wf_groups_copy(igroup)%win_no(ik), ik)), &
                          maxval(wf_groups_copy(igroup)%win_ii(1:wf_groups_copy(igroup)%win_ni(ik), ik)) )
            nst_ik = lst_ik - fst_ik + 1

            ! set window_mask for later assigining state ist to outer/inner window. 
            ! entry of window_mask 0 for state in outer window, 1 for state in inner window.
            allocate( window_mask(nst_ik) )
            window_mask( wf_groups_copy(igroup)%win_io(1:wf_groups_copy(igroup)%win_no(ik), ik) - fst_ik + 1 ) = 0
            window_mask( wf_groups_copy(igroup)%win_ii(1:wf_groups_copy(igroup)%win_ni(ik), ik) - fst_ik + 1 ) = 1

            ! extract sz_eval by checking wf_index_map entries, and set sz_evec to identity matrix
            allocate( sz_eval(nst_ik), source = 0.d0 )
            allocate( sz_evec(nst_ik, nst_ik), source = zzero )
            do i = 1, nst_ik
              sz_evec(i, i) = 1.d0
              if ( wf_index_map(fst_ik + i - 1, ik) <= nstsv/2 ) then
                sz_eval(i) = 1.d0
              else
                sz_eval(i) = -1.d0
              end if
            end do

            ! sort sz_eval accoridng to spin-eigenvalue .
            ! sz_eval, window_mask entries sand sz_evec columns are rearranged, such that they
            ! correspond to states with spin-eigenvalues in ascending order.
            allocate( idx(nst_ik) )
            idx = sort_index_1d( nst_ik, sz_eval )
            sz_eval = sz_eval(idx)
            sz_evec = sz_evec(:, idx)
            window_mask = window_mask(idx)
            ! calculate number of down/up states
            n_down = count( sz_eval < 0.d0, dim=1 )
            n_up = count( sz_eval >= 0.d0, dim=1 )

            ! Sorting window_mask within spin up resp. down subarray, in ascending order (window_map entries: 0 for outer, 1 for inner).
            ! When sz_eval entries and sz_evec columns corresponding to states of spin/down
            ! are rearranged according to idx (which sorts "window mask"-subarray),
            ! they will be re-ordered such that they contain first outer, then inner window states.

            ! down
            idx(:n_down) = sort_index_1d( n_down, window_mask(:n_down) )
            sz_evec(:, :n_down) = sz_evec(:, idx(:n_down))
            sz_eval(:n_down) = sz_eval(idx(:n_down))

            ! up
            idx(:n_up) = sort_index_1d( n_up, window_mask(n_down+1:) )
            sz_evec(:, n_down+1:) = sz_evec(:, n_down+idx(:n_up))
            sz_eval(n_down+1:) = sz_eval(n_down+idx(:n_up))

            ! for new down/up wf_group, set no/ni & io/ii
            ! down group
            no_down = count( window_mask(:n_down) == 0, dim=1 )
            wf_groups(igroup_down)%win_no(ik) = no_down
            wf_groups(igroup_down)%win_io(:no_down, ik) = [(wf_groups(igroup_down)%fst+i-1, i=1, no_down)]
            ni_down = count( window_mask(:n_down) == 1, dim=1 )
            wf_groups(igroup_down)%win_ni(ik) = ni_down
            wf_groups(igroup_down)%win_ii(:ni_down, ik) = [(wf_groups(igroup_down)%fst+no_down+i-1, i=1, ni_down)]
            ! up group
            no_up = count( window_mask(n_down+1:) == 0, dim=1 ) 
            wf_groups(igroup_up)%win_no(ik) = no_up
            wf_groups(igroup_up)%win_io(:no_up, ik) = [(wf_groups(igroup_up)%fst+i-1, i=1, no_up)]
            ni_up = count( window_mask(n_down+1:) == 1, dim=1 ) 
            wf_groups(igroup_up)%win_ni(ik) = ni_up
            wf_groups(igroup_up)%win_ii(:ni_up, ik) = [(wf_groups(igroup_up)%fst+no_up+i-1, i=1, ni_up)]

            ! set group "sz eigenvectors"
            wf_groups(igroup_down)%sz_eigenvector(fst_ik:lst_ik, :n_down, ik) = sz_evec(:, :n_down)
            wf_groups(igroup_up)%sz_eigenvector(fst_ik:lst_ik, :n_up, ik) = sz_evec(:, n_down+1:)

            deallocate( idx, sz_eval, sz_evec, window_mask )
          end do
        end if!TODO error catch for unsupported methods
        
      
        !set new ranges
        wf_groups(igroup_down)%nst = maxval( wf_groups(igroup_down)%win_ni + wf_groups(igroup_down)%win_no)
        wf_groups(igroup_down)%lst = wf_groups(igroup_down)%fst + wf_groups(igroup_down)%nst - 1
        
        wf_groups(igroup_up)%nst = maxval( wf_groups(igroup_up)%win_ni + wf_groups(igroup_up)%win_no)
        wf_groups(igroup_up)%lst = wf_groups(igroup_up)%fst + wf_groups(igroup_up)%nst - 1
      end do

      ! set new nwf
      wf_nwf = 0
      write(*,*)
      write(*,'("Info (wfspin_spindis_divide_groups): Auto-setting new nwf for spin-disentangled groups.")')
      do igroup = 1, wf_ngroups
        ! re-set group fwf, lwf, nwf and wf_nwf
        ! auto select nwf for each group
        wf_groups( igroup)%nwf = nint( 0.5d0*(maxval( wf_groups( igroup)%win_ni) + minval( wf_groups( igroup)%win_ni + wf_groups( igroup)%win_no)))
        wf_groups( igroup)%fwf = wf_nwf + 1
        wf_nwf = wf_nwf + wf_groups( igroup)%nwf
        wf_groups( igroup)%lwf = wf_nwf
        ! print to stdout new group nwf
        write(*,'("Info (wfspin_spindis_divide_groups): Group ",I2," nwf = ",I3)') igroup, wf_groups(igroup)%nwf

        ! check compatibility of number of Wannier functions with windows
        do ik = 1, wf_kset%nkpt
          if( wf_groups( igroup)%win_no( ik) + wf_groups( igroup)%win_ni( ik) .lt. wf_groups( igroup)%nwf) then
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write( *, '("Error (wfspin_spindis_divide_groups): Outer window contains less than nwf (",I3,") bands for k-point ",3F13.6," in spin group ",I2,".")') &
                  wf_groups( igroup)%nwf, wf_kset%vkl( :, ik), igroup
            end if
            stop
          end if
          if( wf_groups( igroup)%win_ni( ik) .gt. wf_groups( igroup)%nwf) then
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write( *, '("Error (wfspin_spindis_divide_groups): Inner window contains more than nwf (",I3,") bands for k-point ",3F13.6," in spin group ",I2,".")') &
                  wf_groups( igroup)%nwf, wf_kset%vkl( :, ik), igroup
            end if
            stop
          end if
        end do
      end do
      return
    end subroutine wfspin_spindis_divide_groups

    !> calculate V^H M V at ik for neighbor idxn. V ist the spin-disentangled "S_z-eigevector", M is the "plane wave matrix" of initial ks-states
    subroutine wfspin_m( ik, idxn, m0_spin )
      use constants, only: zone, zzero

      integer, intent(in) :: ik
      integer, intent(in) :: idxn
      complex(8), intent(out) :: m0_spin(wf_groups(wf_group)%nst, wf_groups(wf_group)%nst)
      
      complex(8), allocatable :: auxmat(:,:)

      allocate( auxmat(  wf_groups( wf_group)%nst_ks, wf_groups( wf_group)%nst) )

      call zgemm( 'n', 'n', wf_groups( wf_group)%nst_ks, wf_groups( wf_group)%nst, wf_groups( wf_group)%nst_ks, zone, &
                  wf_m0( wf_groups( wf_group)%fst_ks, wf_groups( wf_group)%fst_ks, ik, idxn), wf_nst, &
                  wf_groups( wf_group)%sz_eigenvector( wf_groups( wf_group)%fst_ks, 1, wf_n_ik( idxn, ik)), wf_groups( wf_group)%nst_ks, zzero, &
                  auxmat, wf_groups( wf_group)%nst_ks)
      call zgemm( 'c', 'n', wf_groups( wf_group)%nst, wf_groups( wf_group)%nst, wf_groups( wf_group)%nst_ks, zone, &
                  wf_groups( wf_group)%sz_eigenvector( wf_groups( wf_group)%fst_ks, 1, ik), wf_groups( wf_group)%nst_ks, &
                  auxmat, wf_groups( wf_group)%nst_ks, zzero, &
                  m0_spin, wf_groups( wf_group)%nst)
      return
    end subroutine wfspin_m

    !> Back transformation from spin-disentangled states to original SV states.
    subroutine wfspin_undo_spindis
      use precision, only : dp
      use constants, only : zzero, zone

      integer :: igroup, ik

      complex(kind=dp), allocatable :: aux(:, :)

      if (.not. wf_spin_dis) return
      
      allocate( aux(wf_fst:wf_lst, wf_nwf) )
      do igroup = 1, wf_ngroups
        do ik = 1, wf_kset%nkpt
          aux = wf_transform(:, :, ik)
          call zgemm( 'n', 'n', wf_groups(igroup)%nst_ks, wf_groups(igroup)%nwf, wf_groups(igroup)%nst, zone, &
            wf_groups(igroup)%sz_eigenvector(:, :, ik), wf_groups(igroup)%nst_ks, &
            aux(wf_groups(igroup)%fst, wf_groups(igroup)%fwf), wf_nst, zzero, &
            wf_transform(wf_groups(igroup)%fst_ks, wf_groups(igroup)%fwf, ik), wf_nst )
        end do
      end do
    end subroutine wfspin_undo_spindis

    !> Calculates the $\hat{S}_z$ operator in the wannier gauge
    !> $$ (S_z)_{\mu\nu} = \sum_{{\bf k}} w_{\bf k} \sum_{\mu'\nu'}^{\mathcal{J_{\bf k}}}
    !> \sum_{m}^{N_{\text{BSV}}} (U_{\mu'\mu}^{{\bf k}})^\dagger  \left[
    !> (\widetilde{C_{m\mu'}^{{\bf k}\uparrow}})^\dagger \widetilde{C_{m\nu'}^{{\bf k}\uparrow}} 
    !> - (\widetilde{C_{m\mu'}^{{\bf k}\downarrow}})^\dagger \widetilde{C_{m\nu'}^{{\bf k}\downarrow}}
    !> \right] U_{\nu'\nu}^{{\bf k}}$$
    subroutine wfspin_sz_wannier( spin_z )
      use constants, only: zone, zzero
      complex(8), intent(out) :: spin_z(wf_nwf, wf_nwf)
      complex(8), allocatable :: evecsv(:, :, :), auxevecsv(:, :)
      complex(8), allocatable :: auxmat(:, :), spin_z_ik(:, :)
      integer :: ik, igroup, igroup_fwf, igroup_nwf

      ! SVLO failsafe
      if( issvlo() ) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *,'("Error (wfspin_sz_sv): no svlo capability yet!")')
          stop
        end if
      end if

      allocate( evecsv(nstfv, nstsv, nspinor), &
                auxevecsv(nstsv, nstsv) ) 
      allocate( spin_z_ik(wf_nwf, wf_nwf) )
      spin_z = zzero

      do ik = 1, wf_kset%nkpt
          spin_z_ik = zzero
          do igroup = 1, wf_ngroups
            call wfhelp_getevecsv(ik, auxevecsv)
            evecsv(:, :, 1) = auxevecsv(:nstsv/2, :)
            evecsv(:, :, 2) = auxevecsv(nstsv/2+1:, :)
            !"down" components
            ! B := C^down U
            igroup_fwf = wf_groups(igroup)%fwf
            igroup_nwf = wf_groups(igroup)%nwf
            allocate( auxmat(nstfv, igroup_nwf) )

            call zgemm('n', 'n', nstfv, igroup_nwf, wf_nst, zone, &
                      evecsv(1, wf_fst, 2), nstfv, &
                      wf_transform(wf_fst, igroup_fwf, ik), wf_nst, zzero, &
                      auxmat(:,:), nstfv)
            ! (B)* B
            call zgemm('c', 'n', igroup_nwf, igroup_nwf, nstfv, zone, &
                      auxmat(1, 1), nstfv, &
                      auxmat(1, 1), nstfv, zzero, &
                      spin_z_ik(igroup_fwf,igroup_fwf), wf_nwf)
            
            !"up" components
            ! A := C^up U
            call zgemm('n', 'n', nstfv, igroup_nwf, wf_nst, zone, &
                      evecsv(1, wf_fst, 1), nstfv, &
                      wf_transform(wf_fst, igroup_fwf, ik), wf_nst, zzero, &
                      auxmat(:,:), nstfv)
            ! (A)* A - (B)* B
            call zgemm('c', 'n', igroup_nwf, igroup_nwf, nstfv, zone, &
                      auxmat(1, 1), nstfv, &
                      auxmat(1, 1), nstfv, -zone, &
                      spin_z_ik(igroup_fwf,igroup_fwf), wf_nwf)

            deallocate( auxmat )
        end do
        spin_z = spin_z + wf_kset%wkpt( ik ) * spin_z_ik
      end do
      return
    end subroutine wfspin_sz_wannier
end module mod_wannier_spin
