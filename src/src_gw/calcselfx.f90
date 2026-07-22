
module calculate_exchange_self_energy

  use constants, only: zzero, zone, pi, real_zero, fourpi
  use modinput, only: input
  use mod_atoms, only: idxas, natmtot
  use mod_bands, only: evalfv, eveck, eveckp, eveckalm, eveckpalm, numin, nomax
  use mod_core_states, only: ncg, corind
  use mod_coulomb_potential, only: barc, barcev, vccut, vmat 
  use mod_eigensystem, only: nmatmax
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_APW_LO, only: apwordmax
  use mod_get_eigenvectors_times_matchingcoefficients, only: get_eigenvectors_times_matchingcoefficients
  use mod_muffin_tin, only: lmmaxapw
  use mod_product_basis, only: matsiz, mbsiz, minmmat
  use mod_misc_gw, only: vi, Gamma
  use mod_selfenergy, only: singc2, selfex
  use modgw, only: ibgw, nbgw, kset, kqset, Gkqset, ciw, kiw, fdebug, time_selfx
  use modmpi, only: rank
  use mod_offdiagonal_selfenergy, only: add_q_contrib_to_offdiagonal_selfenergy_exchange_at_ik
#include "mod_gw_degeneracies.inc"
  use mod_gw_degeneracies, only: get_degenerate_limits_qp_interval_ikp, &
                                ibgw_including_degeneracy, &
                                nbgw_including_degeneracy, &
                                degenerate_subspaces
  use precision, only: i32, dp
  use mod_expand_products, only: expand_products_generic
#include "offload.fpp"

  implicit none 

  private 

  public :: calcselfx

contains
  !> Calculate the q-dependent self-energy contribution for the given k-points.
  !> Remark: To obtain the complete self-energy, it must be then summed over all q-points
  subroutine calcselfx(iq, ikp_first, ikp_last, offdiagonal)

      implicit none

      !> index of the q-point term to evaluate
      integer(i32), intent(in) :: iq
      !> index of the first k-point (in the reduced BZ) to evaluate the self-energy
      integer(i32), intent(in) :: ikp_first
      !> index of the last k-point (in the reduced BZ) to evaluate the self-energy
      integer(i32), intent(in) :: ikp_last
      !> Compute the offdiagonal terms of the self-energy
      logical, optional, intent(in) :: offdiagonal

      ! local
      integer(i32) :: ik, ikp, jk
      integer(i32) :: mdim
      integer(i32) :: ie1, ie2, im
      integer(i32) :: ia, is, ias, ic, icg
      integer(i32) :: ispace_init, ispace_final, ispace, lowband, upband, size_deg
      real(dp)     :: tstart, tend
      real(dp)     :: sxs2, fnk
      complex(dp)  :: sx, vc
      complex(dp)  :: mvm     ! Sum_ij{M^i*V^c_{ij}*conjg(M^j)}
      logical      :: offdiagonal_local

      ! external routine
      complex(dp), external :: zdotc

      call timesec(tstart)

      if (present(offdiagonal)) then
        offdiagonal_local = offdiagonal
      else
        offdiagonal_local = .false.
      end if

      ! singular term prefactor (q->0)
      sxs2 = fourpi * vi
      if( vccut ) sxs2 = real_zero

      !----------------------------------------
      ! Set v-diagonal mixed product basis set
      !----------------------------------------
      if (vccut) then
        call setbarcev( real_zero, .false. )
      else
        call setbarcev( real_zero, Gamma )
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
      
      OMP_OFFLOAD target enter data map(alloc: eveckalm, eveckpalm, eveck, eveckp)

      allocate(minmmat(mbsiz,ibgw_including_degeneracy:nbgw_including_degeneracy,1:mdim), source=zzero)

      !================================
      ! loop over irreducible k-points
      !================================
      do ikp = ikp_first, ikp_last
        ! k vector
        ik = kset%ikp2ik(ikp)
        ! k-q vector
        jk = kqset%kqid(ik,iq)

        ! get KS eigenvectors
        call get_evec_gw(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), eveckp)
        eveckp = conjg(eveckp)
        call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), eveck)

        call get_eigenvectors_times_matchingcoefficients(ik, 't', eveck, eveckalm)
        call get_eigenvectors_times_matchingcoefficients(jk, 'c', eveckp, eveckpalm)

        OMP_OFFLOAD target update to(eveck, eveckp, eveckalm, eveckpalm)
        
        !========================================================
        ! Calculate the contribution to the exchange self-energy
        !========================================================

        ! We work with the degenerate subspaces and compute the average for 
        ! each of them. See init_dft_eigenvalues.f90 for further iformation on why
        ! this is important to preserve symmetry and, consequently, degeneracies.

        ! First we compute indices of the subspaces we are interested in
        call get_degenerate_limits_qp_interval_ikp(ikp, ispace_init, ispace_final)

        ! Calculate M^i_{nm}+M^i_{cm}
        OMP_OFFLOAD target data map(alloc: minmmat)
        call expand_products_generic(ik, iq, ibgw_including_degeneracy, nbgw_including_degeneracy, 1, 0, 1, nomax, 1, mdim-nomax, minmmat, .true.)
        OMP_OFFLOAD target update from(minmmat)
        OMP_OFFLOAD end target data


        !$omp parallel default(none) private(ie1,ie2,mvm,icg,is,ia,ias,ic,fnk,sx,lowband,upband,size_deg), & 
        !$omp shared(ispace_init,ispace_final,degenerate_subspaces,ikp,mdim,nomax,mbsiz,minmmat,kiw,corind,idxas), &
        !$omp shared(ciw,Gamma,kqset,singc2,selfex,ibgw,nbgw,jk,ik,sxs2)
        ! The different subspaces can have different sizes, and thus different computational cost. Therefore, the scheduler is set dynamic
        !$omp do schedule(dynamic)
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
                mvm = dot_product(minmmat(:,ie1,ie2), minmmat(:,ie1,ie2))
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
                mvm = dot_product(minmmat(:,ie1,ie2), minmmat(:,ie1,ie2))
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
          ! they are not printed). Macro defined in mod_gw_degeneracies.inc (uses ibgw,nbgw)
          selfex(QP_ADJUST_RANGE(lowband,upband),ikp) = &
              selfex(QP_ADJUST_RANGE(lowband,upband),ikp) + sx / size_deg

        end do
        !$omp end do
        !$omp end parallel

        ! If required compute the offdiagonal terms
        ! of the exchange self-energy
        if (offdiagonal_local) then
            call add_q_contrib_to_offdiagonal_selfenergy_exchange_at_ik(ikp, ispace_init, ispace_final, mdim, nomax, jk, &
                                                                        degenerate_subspaces, minmmat, kiw, corind, idxas, ciw)
        end if

        ! debugging info for the diagonal terms
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

      OMP_OFFLOAD target exit data map(delete: eveck, eveckp, eveckalm, eveckpalm)

      deallocate(eveck)
      deallocate(eveckp)
      deallocate(eveckalm)
      deallocate(eveckpalm)

      ! timing
      call timesec(tend)
      time_selfx = time_selfx+tend-tstart

  end subroutine calcselfx

end module calculate_exchange_self_energy
