
subroutine calcselfx(iq)
!
! Calculate the q-dependent self-energy contribution
!
    use modinput, only: input
    use mod_atoms, only: idxas, natmtot
    use mod_eigenvalue_occupancy, only: nstfv
    use mod_APW_LO, only: apwordmax
    use mod_muffin_tin, only: lmmaxapw
    use mod_eigensystem, only: nmatmax
    use mod_bands, only: numin, nomax, eveckalm, eveckpalm, eveck, eveckp, evalfv
    use mod_product_basis, only: matsiz, mbsiz, minmmat
    use mod_core_states, only: ncg, corind
    use mod_coulomb_potential, only: barc, vccut, barcev, vmat
    use mod_misc_gw, only: vi, Gamma
    use mod_mpi_gw, only : myrank
    use modgw, only: kset, kqset, Gkqset, ciw, kiw, fdebug, time_selfx
    use mod_selfenergy, only: singc2, selfex
#include "mod_gw_degeneracies.inc"
    use mod_gw_degeneracies, only: get_degenerate_limits_qp_interval_ikp, &
                                   ibgw_including_degeneracy, &
                                   nbgw_including_degeneracy, &
                                   degenerate_subspaces
    use precision, only: i32, dp
    use constants, only: zone, zzero, real_zero, pi, fourpi

    implicit none

    ! input/output
    integer(i32), intent(in) :: iq

    ! local
    integer(i32) :: ik, ikp, jk, i
    integer(i32) :: mdim
    real(dp)    :: tstart, tend, t0, t1
    integer(i32) :: ie1, ie2, im
    integer(i32) :: ia, is, ias, ic, icg
    real(dp)    :: sxs2, fnk
    complex(dp) :: sx, vc
    complex(dp) :: mvm     ! Sum_{ij}{M^i*V^c_{ij}*conjg(M^j)}
    complex(dp), allocatable :: evecfv(:,:)
    ! For the averaging over degenerate states
    integer(i32) :: ispace_init, ispace_final, ispace, lowband, upband, size_deg
  
    ! external routines
    complex(dp), external :: zdotc

    call timesec(tstart)

    ! singular term prefactor (q->0)
    sxs2 = fourpi * vi

    !----------------------------------------
    ! Set v-diagonal mixed product basis set
    !----------------------------------------
    if (vccut) then
        sxs2 = real_zero
        mbsiz = matsiz
        if (allocated(barc)) deallocate(barc)
        allocate(barc(matsiz,mbsiz), source=zzero)
        do im = 1, matsiz
            if (barcev(im) > 0.0_dp) then
                vc = cmplx(barcev(im), 0.0_dp, kind=dp)
                barc(:,im) = vmat(:,im) * sqrt(vc)
            end if
        end do
    else
        call setbarcev(real_zero)
    end if

    !--------------------------------------------------
    ! total number of states (n->m + n->c transisions)
    !--------------------------------------------------
    if ((input%gw%coreflag=='all').or. &
    &   (input%gw%coreflag=='xal')) then
      mdim = nomax+ncg
    else
      mdim = nomax
    end if

    allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstfv))
    allocate(eveckp(nmatmax,nstfv))

    allocate(minmmat(mbsiz,ibgw_including_degeneracy:nbgw_including_degeneracy,1:mdim), source=zzero)
    ! msize = sizeof(minmmat)*b2mb
    ! write(*,'(" calcselfx: rank, size(minmmat) (Mb):",i4,f12.2)') myrank, msize

    !================================
    ! loop over irreducible k-points
    !================================
    ! write(*,*)
    do ikp = 1, kset%nkpt
      ! write(*,*) 'calcselfx: rank, (iq, ikp):', myrank, iq, ikp

      ! k vector
      ik = kset%ikp2ik(ikp)
      ! k-q vector
      jk = kqset%kqid(ik,iq)

      ! get KS eigenvectors
      allocate(evecfv(nmatmax,nstfv))
      call get_evec_gw(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), evecfv)
      eveckp = conjg(evecfv)
      call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
      call move_alloc(evecfv, eveck)

      call expand_evec(ik, 't')
      call expand_evec(jk, 'c')

      ! Obtain the limits for degenerate subspaces for the irreducible point
      call get_degenerate_limits_qp_interval_ikp(ikp, ispace_init, ispace_final)

      ! Calculate M^i_{nm}+M^i_{cm}
      call expand_products(ik, iq, ibgw_including_degeneracy, nbgw_including_degeneracy, -1, 1, mdim, nomax, minmmat)

      !========================================================
      ! Calculate the contribution to the exchange self-energy
      !========================================================

      ! We work with the degenerate subspaces and compute the average for 
      ! each of them. See init_dft_eigenvalues.f90 for further iformation on why
      ! this is important to preserve symmetry and, consequently, degeneracies.

      ! First we compute indices of the subspaces we are interested in
      call get_degenerate_limits_qp_interval_ikp(ikp, ispace_init, ispace_final)

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ie1,ie2,mvm,icg,is,ia,ias,ic,fnk,sx,lowband,upband,size_deg), & 
!$OMP SHARED(ispace_init,ispace_final,degenerate_subspaces,ikp,mdim,nomax,mbsiz,minmmat,kiw,corind,idxas), &
!$OMP SHARED(ciw,Gamma,kqset,singc2,selfex,ibgw,nbgw,jk,ik,sxs2)
! The different subspaces can have different sizes, and thus different computational cost. Therefore, the scheduler is set dynamic
!$OMP DO SCHEDULE(DYNAMIC)
#endif
      do ispace = ispace_init, ispace_final

        lowband  = degenerate_subspaces(1, ispace, ikp)
        upband   = degenerate_subspaces(2, ispace, ikp)
        size_deg = degenerate_subspaces(3, ispace, ikp)

        sx = zzero

        do ie1 = lowband, upband
          ! sum over occupied states
          do ie2 = 1, mdim
            !=======================
            ! Valence contribution
            !=======================
            if (ie2 <= nomax) then
              mvm = zdotc(mbsiz, minmmat(:,ie1,ie2), 1, minmmat(:,ie1,ie2), 1)
              sx = sx - kiw(ie2,jk)*mvm
            else
              !=============================
              ! Core electron contribution
              !=============================
              icg = ie2-nomax
              is = corind(icg,1)
              ia = corind(icg,2)
              ias = idxas(ia,is)
              ic = corind(icg,3)
              mvm = zdotc(mbsiz, minmmat(:,ie1,ie2), 1, minmmat(:,ie1,ie2), 1)
              sx = sx - ciw(ic,ias)*mvm
            end if ! occupied states
          end do ! ie2
          ! add singular term (q->0)
          if ((Gamma) .and. (ie1 <= nomax)) then
            ! occupation number
            fnk = kiw(ie1,ik) * kqset%nkpt
            sx  = sx - sxs2 * fnk * singc2
          end if

        end do 
        ! That ensures we are in the proper range (in that way states out of
        ! the print range are taken into account for degeneracy stuff, but 
        ! they are not printed). Macro defined in mod_gw_degeneracies.inc
        selfex(QP_ADJUST_RANGE(lowband,upband),ikp) = &
            selfex(QP_ADJUST_RANGE(lowband,upband),ikp) + sx / size_deg

      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

      ! debugging info
      if (input%gw%debug) then
        write(fdebug,*) 'EXCHANGE SELF-ENERGY: iq=', iq, ' ikp=', ikp
        write(fdebug,*) 'state   Sigma_x'
        do ie1 = ibgw, nbgw
          write(fdebug,'(i0,x,(3e16.8))') ie1, real(selfex(ie1,ikp)), aimag(selfex(ie1,ikp)), evalfv(ie1,ikp)
        end do
        write(fdebug,*)
      end if

    end do ! ikp

    deallocate(minmmat)
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)

    ! timing
    call timesec(tend)
    time_selfx = time_selfx+tend-tstart

    return
end subroutine
!EOC
