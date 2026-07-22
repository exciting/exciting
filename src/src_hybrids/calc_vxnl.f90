!
!BOP
! !ROUTINE: calc_vxnl
! !INTERFACE:
!
subroutine calc_vxnl()
! !USES:
    use mod_atoms, only: idxas, natmtot
    use mod_misc_gw, only: vi, gammapoint, Gamma
    use mod_core_states, only: ncg, corind
    use modinput, only: input
    use constants, only: pi, zzero
    use mod_product_basis, only: mpwipw, locmatsiz, mbsiz, matsiz
    use modgw, only: kset, kqset, Gqset, Gkqset, b2mb, mblksiz, kiw, ciw
    use mod_selfenergy, only: singc2
    use mod_hybrids, only: vxnl, exnl, vnlmat
    use mod_coulomb_potential, only: barc
    use mod_APW_LO, only: apwordmax
    use mod_get_eigenvectors_times_matchingcoefficients, only: get_eigenvectors_times_matchingcoefficients
    use mod_muffin_tin, only: lmmaxapw
    use mod_eigensystem, only: nmatmax
    use mod_eigenvalue_occupancy, only: efermi, nstfv, occsv
    use mod_coulomb_potential, only: delete_coulomb_potential
    use modmpi
    use mod_bands, only: evalfv, nomax, numin, ikvbm, ikcbm, ikvcm, eveck, eveckalm, eveckp, eveckpalm
    use cdft, only: cdft_input_keys
    use general_find_vbm_cbm, only: find_vbm_cbm
    use precision, only: i32, dp
    use mod_expand_products, only: expand_products_generic, split_interval
#include "offload.fpp"

!
! !DESCRIPTION:
!   Calculates the non-local exchange potential
!   and the non-local exchange energy for Hartree-Fock based hybrid functionals.
!
!EOP
!BOC
    implicit none

    integer(i32) :: ikp, ik, jk, iq, ikq, jkp
    integer(i32) :: ie12, ie12tot, ie1, ie2, ie3, icg
    integer(i32) :: ist, l, im
    integer(i32) :: i, j, k, jst, ispn
    integer(i32) :: is, ia, ias, ic
    integer(i32) :: n, nmdim, m, mdim
    integer(i32) :: iblk, nblk, mstart, mend
    real(dp)     :: tstart, tend, sxs2, msize
    complex(dp)  :: mvm
    integer(i32) :: ikfirst, iklast

    integer(i32), allocatable :: idxpair(:,:)
    complex(dp), allocatable :: minm(:,:,:)
    complex(dp), allocatable :: evecsv(:,:)
    
    integer(i32) :: m_val_start, m_val_end, m_core_start, m_core_end
    type(cdft_input_keys) :: cdft_calculation
    real(dp), parameter :: tolerance = 1.0e-8_dp

    call cpu_time(tstart)

    !----------------------------------------
    ! Read KS eigenvalues from file EVALSV.OUT
    !----------------------------------------
    if (allocated(evalfv)) deallocate(evalfv)
    allocate(evalfv(nstfv,kset%nkpt))
    evalfv(:,:) = 0.0_dp
    do ik = 1, kset%nkpt
      call getevalfv(kset%vkl(:,ik), evalfv(:,ik))
    end do

    call cdft_calculation%read_input_keys( input%groundstate )
    if (cdft_calculation%is_on()) then
      call find_vbm_cbm(1, nstfv, kset%nkpt, occsv, evalfv, nomax, numin, ikvbm, ikcbm, ikvcm)
    else
      ! VB / CB state index in ground state
      call find_vbm_cbm(1, nstfv, kset%nkpt, evalfv, efermi, nomax, numin, ikvbm, ikcbm, ikvcm)
    end if

    ! BZ integration weights
    call kintw()
    deallocate(evalfv)

    ! singular term prefactor
    sxs2 = 4.0_dp*pi*vi*singc2*kqset%nkpt

    if ((input%gw%coreflag=='all').or. &
        (input%gw%coreflag=='xal')) then
      mdim = nomax+ncg
    else
      mdim = nomax
    end if

    allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstfv))
    allocate(eveckp(nmatmax,nstfv))
    OMP_OFFLOAD target enter data map(alloc: eveckalm, eveckpalm, eveck, eveckp)

    !------------------------------------------!
    ! Matrix elements of non-local potential   !
    !------------------------------------------!
    if (allocated(vxnl)) deallocate(vxnl)
    allocate(vxnl(nstfv,nstfv,kset%nkpt))
    vxnl(:,:,:) = zzero

    ie12tot = nstfv*(nstfv+1)/2
    allocate(idxpair(2,ie12tot))
    ie12 = 0
    do ie1 = 1, nstfv
      do ie2 = ie1, nstfv
        ie12 = ie12+1
        idxpair(1,ie12) = ie1
        idxpair(2,ie12) = ie2
      end do
    end do

    !-------------------------------------------------------
    ! determine the number of blocks used in minm operation
    !-------------------------------------------------------
    if (mblksiz >= mdim) then
      nblk = 1
    else
      nblk = mdim / mblksiz
      if (mod(mdim, mblksiz) /= 0) nblk = nblk+1
    end if
    if ((input%groundstate%outputlevelnumber>1) .and. (rank==0)) then
      write(60,*) 'Info(calc_vxnl):'
      write(60,'(a,3i8)') '    mdim, nblk, mblksiz: ', mdim, nblk, mblksiz
    end if

    !---------------------------------------
    ! Loop over k-points
    !---------------------------------------
    ikq = 0
    do ikp = 1, kset%nkpt
      !---------------------------------------
      ! Integration over BZ
      !---------------------------------------
      do iq = 1, kqset%nkpt
        Gamma = gammapoint(kqset%vqc(:,iq))
        ik  = kset%ikp2ik(ikp)
        jk  = kqset%kqid(ik,iq)
        jkp = kset%ik2ikp(jk)

        !=======================================
        ! distribute (k,q)-pair over processors
        !=======================================
        ikq = ikq + 1
        if (mod(ikq, procs) == rank) then

          matsiz = locmatsiz+Gqset%ngk(1,iq)
          call diagsgi(iq)
          call calcmpwipw(iq)

          !------------------------------------
          ! Calculate the bare Coulomb matrix
          !------------------------------------
          call calcbarcmb(iq)
          call setbarcev(0.0_dp, Gamma)

          if ((input%groundstate%outputlevelnumber>1) .and. (rank==0)) then
            write(60,'(a,3i8)')  '    ---> rank, ikp, iq = ', rank, ikp, iq
            write(60,'(a,4i8)') '    locmatsiz, ngk, matsiz, mbsiz:', &
                                locmatsiz, Gqset%ngk(1,iq), matsiz, mbsiz
          end if

          !------------------------------------------------------------
          ! k-q vector
          call getevecfv(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), eveck)
          eveckp = conjg(eveck)
          ! k vector
          call getevecfv(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), eveck)
          call get_eigenvectors_times_matchingcoefficients(ik, 't', eveck, eveckalm)
          call get_eigenvectors_times_matchingcoefficients(jk, 'c', eveckp, eveckpalm)
          
          OMP_OFFLOAD target update to(eveckalm, eveckpalm, eveck, eveckp)

          !=================================
          ! Loop over m-blocks in M^i_{nm}
          !=================================
          do iblk = 1, nblk

            mstart = 1 + (iblk-1)*mblksiz
            mend = min(mdim, mstart+mblksiz-1)

            ! m-block M^i_{nm}
            allocate(minm(mbsiz,1:nstfv,mstart:mend))
            OMP_OFFLOAD target enter data map(alloc: minm)
            if ((input%groundstate%outputlevelnumber>1) .and. (rank==0)) then
              msize = sizeof(minm)*b2mb
              write(60,'(a,3i8,f14.2)') '    iblk, mstart, mend, size(minm) (Mb):', &
                                        iblk, mstart, mend, msize
            end if

            !---------------------
            ! Calculate M^i_{nm}
            !---------------------
            call split_interval( mstart, mend, nomax, m_val_start, m_val_end, m_core_start, m_core_end)
            call expand_products_generic(ik, iq, 1, nstfv, 1, 0,  m_val_start, m_val_end, m_core_start, m_core_end, minm, .true.)
            OMP_OFFLOAD target update from(minm)

            ! sum over occupied states
            do ie3 = mstart, mend
              if (ie3 <= nomax) then
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ie12, ie1, ie2, mvm)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
                ! TODO (mrm): 
                !
                ! Port this to the device. Please use the following partition scheme:
                ! - ie3 and ie12 : teams as they operate on non-contiguous memory
                ! - dot_product  : use internal execution unit paralellism here (i.e. do parallel). 
                ! See src_gw/calcminm2.f90 MT loop for further inspiration on how to fully port
                ! calc_vxnl to the device, beyond the "expand_products" and "calcminm2".
                ! 
                do ie12 = 1, ie12tot
                  ie1 = idxpair(1,ie12)
                  ie2 = idxpair(2,ie12)
                  mvm = dot_product(minm(1:mbsiz,ie1,ie3), minm(1:mbsiz,ie2,ie3))
                  if (cdft_calculation%is_on()) then
                    ! 1/kqset%nkpt needs to be replaced by weight of kq point in the future.
                    vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp) - 0.5_dp*mvm*occsv(ie3,jkp)/kqset%nkpt
                  else
                    vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp) - kiw(ie3,jk)*mvm
                  end if
                end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
              else
                ! Core electron contribution
                icg = ie3 - nomax
                is  = corind(icg,1)
                ia  = corind(icg,2)
                ias = idxas(ia,is)
                ic  = corind(icg,3)
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ie12, ie1, ie2, mvm)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
                do ie12 = 1, ie12tot
                  ie1 = idxpair(1,ie12)
                  ie2 = idxpair(2,ie12)
                  mvm = dot_product(minm(1:mbsiz,ie1,ie3), minm(1:mbsiz,ie2,ie3))
                  vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp) - ciw(ic,ias)*mvm
                end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
              end if ! core
            end do ! ie3

            OMP_OFFLOAD target exit data map(delete: minm)
            deallocate(minm)

          end do ! iblk

          !--------------------------
          ! add singular term (q->0)
          !--------------------------
          if (Gamma) then
            do ie1 = 1, nomax
              if (cdft_calculation%is_on()) then
                vxnl(ie1,ie1,ikp) = vxnl(ie1,ie1,ikp) - 0.5_dp*sxs2*occsv(ie1,ikp)/kqset%nkpt
              else
                vxnl(ie1,ie1,ikp) = vxnl(ie1,ie1,ikp) - sxs2*kiw(ie1,ik)
              end if
            end do
          end if

        end if ! rank

        call delete_coulomb_potential()

      end do ! iq

    end do ! ikp

    ! clear memory
    OMP_OFFLOAD target exit data map(delete: eveckalm, eveckpalm, eveck, eveckp)
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)
    deallocate(idxpair)

    !! We need to free those for the device accelerated
    !! version.
    if (allocated(mpwipw)) then
      OMP_OFFLOAD target exit data map(delete: mpwipw)
      deallocate(mpwipw)
    end if
    if (allocated(barc)) then
      OMP_OFFLOAD target exit data map(delete: barc)
      deallocate(barc)
    end if

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE, vxnl, nstfv*nstfv*kset%nkpt,  &
                       MPI_DOUBLE_COMPLEX,  MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
#endif

    exnl = 0.d0
    do ikp = 1, kset%nkpt
      do ie1 = 1, nstfv
        do ie2 = ie1+1, nstfv
          vxnl(ie2,ie1,ikp) = conjg(vxnl(ie1,ie2,ikp))
        end do
      end do
      do ie1 = 1, nstfv
        exnl = exnl + 0.5_dp*kset%wkpt(ikp)*occsv(ie1,ikp)*vxnl(ie1,ie1,ikp)
      end do
    end do ! ikp

    call cpu_time(tend)

    return
end subroutine
