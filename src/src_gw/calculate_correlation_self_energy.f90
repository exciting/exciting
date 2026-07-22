module calculate_correlation_self_energy
  use mod_APW_LO, only: apwordmax
  use mod_atoms, only: natmtot
  use mod_bands, only: eveck, eveckp, eveckalm, eveckpalm, n_states => nstse, evalfv
  use mod_dielectric_function, only: epsilon
  use mod_eigensystem, only: nmatmax
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_expand_products, only: expand_products_generic, split_interval
  use mod_get_eigenvectors_times_matchingcoefficients, only: get_eigenvectors_times_matchingcoefficients
  use mod_gw_degeneracies, only: ibgw_including_degeneracy, nbgw_including_degeneracy
  use mod_mpi_gw, only : indexes_parallelization
  use mod_misc_gw, only: Gamma
  use mod_muffin_tin, only: lmmaxapw
  use mod_product_basis, only: minmmat, mbsiz
  use mod_selfenergy, only: mwm, selfec, freq_selfc
  use modinput, only: input
  use modgw, only: kset, kqset, Gkqset, time_selfc, ibgw, nbgw, mblksiz, freq, fdebug
  use precision, only: i32, dp

  implicit none

  private

  public :: calcselfc

  type, public :: sigmac_indexes
    type(indexes_parallelization) :: k_points
    type(indexes_parallelization) :: bands
  end type

contains
!> Obtain the correlation part of the self energy for the given k-points, 
!> evaluating the one term (of a sum) corresponding to a given q-point
subroutine calcselfc( iq, indexes, offdiagonal )
    
#include "offload.fpp"

    !> index of the q-point term to evaluate
    integer(i32), intent(in) :: iq
    !> Set of indexes (k-points, bands) used to calculate sigmac
    type(sigmac_indexes), intent(in) :: indexes
    !> Compute the offdiagonal terms of the self-energy
    logical, optional, intent(in) :: offdiagonal
    
    ! local
    integer(i32) :: ik, ikp, jk, ie1, iom
    integer(i32) :: mdim, iblk, nblk, mstart, mend
    integer(i32) :: first_state, last_state, last_empty_state
    integer(i32) :: m_val_start, m_val_end, m_core_start, m_core_end
    real(dp) :: tstart, tend
    complex(dp), allocatable :: evec_aux(:, :)
    logical :: only_core_states_in_my_rank
    logical :: offdiagonal_local

    call timesec(tstart)

    ! Set the local value for offdiagonal    
    if (present(offdiagonal)) then
      offdiagonal_local = offdiagonal
    else
      offdiagonal_local = .false.
    end if

    ! Update data
    OMP_OFFLOAD target update to(epsilon) 

    !------------------------
    ! total number of states
    !------------------------
    first_state = indexes%bands%my_first
    ! last_state may include core states, depending on choice given in input.xml
    last_state = indexes%bands%my_last
    last_empty_state = min( last_state, n_states )
    mdim = last_state - first_state + 1
    ! Check if this MPI rank only treats core states along the m-dimension
    only_core_states_in_my_rank = (first_state>last_empty_state)

    !-------------------------------------------------------
    ! determine the number of blocks used in minm operation
    !-------------------------------------------------------
    if (mblksiz >= mdim) then
      nblk = 1
    else
      nblk = mdim / mblksiz
      if (mod(mdim,mblksiz) /= 0) nblk = nblk+1
    end if

    !-------------------------------------------
    ! products M*W^c*M
    !-------------------------------------------
    allocate( mwm(ibgw_including_degeneracy:nbgw_including_degeneracy, first_state:last_state, 1:freq%nomeg) )
    allocate( eveckalm(ibgw_including_degeneracy:nbgw_including_degeneracy, apwordmax, lmmaxapw, natmtot) )
    allocate( eveck(nmatmax, ibgw_including_degeneracy:nbgw_including_degeneracy) )
    if( .not. only_core_states_in_my_rank ) then
      allocate( eveckpalm(first_state:last_empty_state, apwordmax, lmmaxapw, natmtot) )
      allocate( eveckp(nmatmax, first_state:last_empty_state) )
      OMP_OFFLOAD target enter data map(alloc: eveckpalm, eveckp)
    end if
    allocate( evec_aux(nmatmax, nstfv) )

    OMP_OFFLOAD target enter data map(alloc: mwm, eveckalm, eveck)

    !================================
    ! loop over irreducible k-points
    !================================
    do ikp = indexes%k_points%my_first, indexes%k_points%my_last
      ! k vector
      ik = kset%ikp2ik(ikp)
      ! k-q vector
      jk = kqset%kqid(ik,iq)

      ! get KS eigenvectors
      if( .not. only_core_states_in_my_rank ) then
        call get_evec_gw( kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), evec_aux )
        eveckp = conjg( evec_aux(:, first_state:last_empty_state) )
      end if 
      call get_evec_gw( kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evec_aux )
      eveck = evec_aux(:, ibgw_including_degeneracy:nbgw_including_degeneracy)

      call get_eigenvectors_times_matchingcoefficients(ik, 't', eveck, eveckalm)
      if( .not. only_core_states_in_my_rank ) then 
        call get_eigenvectors_times_matchingcoefficients(jk, 'c', eveckp, eveckpalm)
        OMP_OFFLOAD target update to(eveckp, eveckpalm)
      end if
      OMP_OFFLOAD target update to(eveck, eveckalm)

      !=================================
      ! Loop over m-blocks in M^i_{nm}
      !=================================
      do iblk = 1, nblk

        mstart = first_state + (iblk-1)*mblksiz
        mend = min( last_state, mstart+mblksiz-1 )

        ! m-block M^i_{nm}
        allocate(minmmat(mbsiz,ibgw_including_degeneracy:nbgw_including_degeneracy,mstart:mend))
        OMP_OFFLOAD target data map(alloc: minmmat)
        call split_interval( mstart, mend, n_states, m_val_start, m_val_end, m_core_start, m_core_end )
        call expand_products_generic(ik, iq, ibgw_including_degeneracy, nbgw_including_degeneracy, 1, 0,  m_val_start, m_val_end, m_core_start, m_core_end, minmmat, .true.)
        ! For Gamma we retrieve the minmmat from the device
        ! because MWM corrections for head and wings are computed 
        ! in the host. This is because their memory layout is not
        ! friendly for the device.
        OMP_OFFLOAD target update from(minmmat) if(Gamma)

        !================================================================
        ! Calculate weight(q)*Sum_ij{M^i*W^c_{ij}(k,q;\omega)*conjg(M^j)}
        !================================================================
        call calcmwm(ikp, jk, ibgw_including_degeneracy, nbgw_including_degeneracy, mstart, mend, minmmat, offdiagonal_local)

        OMP_OFFLOAD end target data ! minmmat 
        deallocate(minmmat)

      end do ! iblk

      !=======================================
      ! Calculate the correlation self-energy
      !=======================================
      if (input%gw%taskname=='cohsex') then
        call calcselfc_cohsex(ikp, iq, mdim)
      else
        if (input%gw%selfenergy%method == 'cd') then
          ! Contour deformation technique
          call calcselfc_freqconv_cd(ikp, iq, mdim)
        else if (input%gw%selfenergy%method == 'ac') then
          ! Imaginary frequency formalism
          call calcselfc_freqconv_ac(ikp, iq)
        end if
      end if

      if (input%gw%debug) then
        write(fdebug,*) 'CORRELATION SELF-ENERGY: iq=', iq, ' ikp=', ikp
        write(fdebug,*) 'state iom  Sigma_c  enk'
        do ie1 = ibgw, nbgw
          do iom = 1, freq_selfc%nomeg
            write(fdebug,*) ie1, iom, abs(selfec(ie1,iom,ikp)), evalfv(ie1,ikp)
          end do
        end do
        write(fdebug,*)
      end if

    end do ! ikp

    OMP_OFFLOAD target exit data map(delete: eveck, eveckalm, mwm)

    deallocate(eveck)
    deallocate(eveckalm)
    deallocate(mwm)
    if( .not. only_core_states_in_my_rank ) then 
      OMP_OFFLOAD target exit data map(delete: eveckp, eveckpalm)
      deallocate(eveckp)
      deallocate(eveckpalm)
    end if

    ! timing
    call timesec(tend)
    time_selfc = time_selfc+tend-tstart

end subroutine
end module
