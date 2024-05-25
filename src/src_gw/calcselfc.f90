!> Obtain the correlation part of the self energy for the given k-points, 
!> evaluating the one term (of a sum) corresponding to a given q-point
subroutine calcselfc(iq, ikp_first, ikp_last)
    use modinput, only: input
    use modgw, only: time_selfc, kqset, kset, Gkqset, b2mb, ibgw, nbgw, freq, mblksiz, msize, fdebug
    use mod_APW_LO, only: apwordmax
    use mod_atoms, only: natmtot
    use mod_bands, only: eveckalm, eveckpalm, eveckp, eveck, nstse, evalfv
    use mod_core_states, only: ncg
    use mod_eigensystem, only: nmatmax
    use mod_eigenvalue_occupancy, only: nstfv
    use mod_gw_degeneracies, only: ibgw_including_degeneracy, nbgw_including_degeneracy
    use mod_muffin_tin, only: lmmaxapw
    use mod_product_basis, only: minmmat, mbsiz
    use mod_selfenergy, only: mwm, freq_selfc, selfec
    use precision, only: i32, dp
    
    implicit none

    !> index of the q-point term to evaluate
    integer(i32), intent(in) :: iq
    !> index of the first k-point (in the reduced BZ) to evaluate the self-energy
    integer(i32), intent(in) :: ikp_first
    !> index of the last k-point (in the reduced BZ) to evaluate the self-energy
    integer(i32), intent(in) :: ikp_last

    ! local
    integer(i32) :: ik, ikp, jk, ie1, iom
    integer(i32) :: mdim, iblk, nblk, mstart, mend
    real(dp) :: tstart, tend

    call timesec(tstart)

    !------------------------
    ! total number of states
    !------------------------
    if (input%gw%coreflag=='all') then
      mdim = nstse+ncg
    else
      mdim = nstse
    end if

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
    allocate(mwm(ibgw_including_degeneracy:nbgw_including_degeneracy,1:mdim,1:freq%nomeg))
    ! msize = sizeof(mwm)*b2mb
    ! write(*,'(" calcselfc: size(mwm) (Mb):",f12.2)') msize

    allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstfv))
    allocate(eveckp(nmatmax,nstfv))

    !================================
    ! loop over irreducible k-points
    !================================
    ! write(*,*)
    do ikp = ikp_first, ikp_last
      ! k vector
      ik = kset%ikp2ik(ikp)
      ! k-q vector
      jk = kqset%kqid(ik,iq)

      ! get KS eigenvectors
      call get_evec_gw(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), eveck)
      eveckp = conjg(eveck)
      call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), eveck)

      call expand_evec(ik, 't')
      call expand_evec(jk, 'c')

      !=================================
      ! Loop over m-blocks in M^i_{nm}
      !=================================
      do iblk = 1, nblk

        mstart = 1 + (iblk-1)*mblksiz
        mend = min(mdim, mstart+mblksiz-1)

        ! m-block M^i_{nm}
        allocate(minmmat(mbsiz,ibgw_including_degeneracy:nbgw_including_degeneracy,mstart:mend))
        msize = sizeof(minmmat)*b2mb

        call expand_products(ik, iq, ibgw_including_degeneracy, nbgw_including_degeneracy, -1, mstart, mend, nstse, minmmat)

        !================================================================
        ! Calculate weight(q)*Sum_ij{M^i*W^c_{ij}(k,q;\omega)*conjg(M^j)}
        !================================================================
        call calcmwm(ibgw_including_degeneracy, nbgw_including_degeneracy, mstart, mend, minmmat)

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
          call calcselfc_freqconv_ac(ikp, iq, mdim)
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

    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)

    ! delete MWM
    deallocate(mwm)

    ! timing
    call timesec(tend)
    time_selfc = time_selfc+tend-tstart

    return
end subroutine
