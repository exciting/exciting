!
!BOP
! !ROUTINE: calc_vxnl
! !INTERFACE:
!
subroutine calc_vxnl()
! !USES:
    use modmain
    use mod_hybrids,           only: kset, vxnl, eveckalm, eveckpalm, eveck, eveckp, &
                                     vxnl, kqset, minmmat, exnl, fgw, Gamma, ikcbm, ikvbm, ikvcm, nomax, numin, &
                                     matsiz, mbsiz, ncg, vi, kiw, ciw, corind, Gqset, Gkqset, singc2, gammapoint, &
                                     locmatsiz, evalfv
    use mod_coulomb_potential, only: barc, delete_coulomb_potential, vccut, vmat, barcev
    use modfvsystem
    use modmpi
!
! !DESCRIPTION:
!   Calculates the non-local exchange potential
!   and the non-local exchange energy for Hartree-Fock based hybrid functionals.
!
!EOP
!BOC
    implicit none

    integer(4) :: ikp, ik, jk, iq, ikp_
    integer(4) :: ie1, ie2, ie3, icg, icg1
    integer(4) :: ist, l, im
    integer(4) :: i, j, k, jst, ispn
    integer(4) :: is, ia, ias, ic
    integer(4) :: n, nmdim, m, mdim
    real(8)    :: tstart, tend, sxs2
    complex(8) :: mvm, zt1, vc
    integer    :: ikfirst, iklast

    complex(8), allocatable :: minm(:,:,:)
    complex(8), allocatable :: evecsv(:,:)

    complex(8), external :: zdotc

    call cpu_time(tstart)
    if (rank == 0) call boxmsg(fgw,'-','Calculate Vx_NL')

    if (allocated(evalfv)) deallocate(evalfv)
    allocate(evalfv(nstfv,kset%nkpt))
    evalfv(:,:) = 0.d0

    !----------------------------------------
    ! Read KS eigenvalues from file EVALSV.OUT
    !----------------------------------------
    do ik = 1, kset%nkpt
      call getevalfv(kset%vkl(:,ik), evalfv(:,ik))
    end do

#ifdef MPI
    ikfirst = firstofset(rank,kset%nkpt)
    iklast  = lastofset(rank,kset%nkpt)
#else
    ikfirst = 1
    iklast  = kset%nkpt
#endif

    !------------------------------------------!
    ! Matrix elements of non-local potential   !
    !------------------------------------------!
    if (allocated(vxnl)) deallocate(vxnl)
    allocate(vxnl(nstfv,nstfv,ikfirst:iklast))
    vxnl(:,:,:) = zzero

    ! VB / CB state index
    call find_vbm_cbm(1,nstfv,nkpt,evalfv,efermi,nomax,numin,ikvbm,ikcbm,ikvcm)
    if (rank==0) then
      write(fgw,'(a,i4)') " Band index of VBM:", nomax
      write(fgw,'(a,i4)') " Band index of CBM:", numin
      write(fgw,*)
    end if
    ! BZ integration weights
    call kintw()
    deallocate(evalfv)

    ! singular term prefactor
    sxs2 = 4.d0*pi*vi*singc2*kqset%nkpt

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
    !---------------------------------------
    ! Loop over k-points
    !---------------------------------------
    do ikp = ikfirst, iklast

      write(fgw,*) 'vxnl: ', rank, ikfirst, iklast, ikp

      !---------------------------------------
      ! Integration over BZ
      !---------------------------------------
      do iq = 1, kqset%nkpt

        Gamma = gammapoint(kqset%vqc(:,iq))

        ! Set the size of the basis for the corresponding q-point
        matsiz = locmatsiz+Gqset%ngk(1,iq)
        call diagsgi(iq)
        call calcmpwipw(iq)

        !------------------------------------
        ! Calculate the bare Coulomb matrix
        !------------------------------------
        call calcbarcmb(iq)
        call setbarcev(0.d0)

        !------------------------------------------------------------
        ! Calculate the M^i_{nm}(k,q) matrix elements for given k and q
        !------------------------------------------------------------
        ik  = kset%ikp2ik(ikp)
        jk  = kqset%kqid(ik,iq)

        ! k-q vector
        call getevecfv(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), eveck)
        eveckp = conjg(eveck)
        ! k vector
        call getevecfv(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), eveck)

        call expand_evec(ik,'t')
        call expand_evec(jk,'c')

        ! M^i_{nm}+M^i_{cm}
        allocate(minmmat(mbsiz,nstfv,1:mdim))
        minmmat(:,:,:) = zzero
        call expand_products(ik, iq, 1, nstfv, -1, 1, mdim, nomax, minmmat)
        call delete_coulomb_potential

        do ie1 = 1, nstfv
          do ie2 = ie1, nstfv
            zt1 = zzero
            do ie3 = 1, mdim
              if (ie3 <= nomax) then
                mvm = zdotc(mbsiz, minmmat(:,ie1,ie3), 1, minmmat(:,ie2,ie3), 1)
                zt1 = zt1 - kiw(ie3,jk)*mvm
              else
                !=============================
                ! Core electron contribution
                !=============================
                icg = ie3 - nomax
                is  = corind(icg,1)
                ia  = corind(icg,2)
                ic  = corind(icg,3)
                mvm = zdotc(mbsiz, minmmat(:,ie1,ie3), 1, minmmat(:,ie2,ie3), 1)
                ias = idxas(ia,is)
                zt1 = zt1 - ciw(ic,ias)*mvm
              end if
            end do ! ie3
            vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp) + zt1
          end do ! ie2
          !__________________________
          ! add singular term (q->0)
          if (Gamma.and.(ie1<=nomax) ) then
            vxnl(ie1,ie1,ikp) = vxnl(ie1,ie1,ikp) - sxs2*kiw(ie1,ik)
          end if
        end do ! ie1

        deallocate(minmmat)

      end do ! iq

      do ie1 = 1, nstfv
        do ie2 = ie1+1, nstfv
          vxnl(ie2,ie1,ikp) = conjg(vxnl(ie1,ie2,ikp))
        end do
      end do

      !------------------------------------------------------------
      ! Debugging Info
      !------------------------------------------------------------
      if (input%gw%debug) then
        call linmsg(fgw,'-','')
        call linmsg(fgw,'-',' Diagonal elements of Vx_NL_nn ')
        write(fgw,*) 'for k-point ', ikp
        do ie1 = 1, nstfv
          write(fgw,'(i4,2f12.4)') ie1, vxnl(ie1,ie1,ikp)
        end do
      end if

    end do ! ikp

    ! clear memory
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)

    !---------------------------------------------------------------
    ! Calculate the non-local exchange energy
    !---------------------------------------------------------------
    exnl = 0.d0
    do ikp = ikfirst, iklast
      do ie1 = 1, nomax
        exnl = exnl+kset%wkpt(ikp)*vxnl(ie1,ie1,ikp)
      end do
    end do
#ifdef MPI
    call MPI_AllReduce(MPI_IN_PLACE, exnl, 1, &
    &                  MPI_DOUBLE_PRECISION, MPI_SUM, &
    &                  MPI_COMM_WORLD, ierr)
    call barrier
#endif

    call cpu_time(tend)
    if (rank==0) call write_cputime(fgw,tend-tstart, 'CALC_VXNL')

    return
end subroutine
