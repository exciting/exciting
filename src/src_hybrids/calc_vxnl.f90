!
!BOP
! !ROUTINE: calc_vxnl
! !INTERFACE:
!
subroutine calc_vxnl()
! !USES:
    use modmain
    use mod_hybrids
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

    integer(4) :: ikp, ik, jk, iq, ikq
    integer(4) :: ie12, ie12tot, ie1, ie2, ie3, icg
    integer(4) :: ist, l, im
    integer(4) :: i, j, k, jst, ispn
    integer(4) :: is, ia, ias, ic
    integer(4) :: n, nmdim, m, mdim
    integer(4) :: iblk, nblk, mstart, mend
    real(8)    :: tstart, tend, sxs2
    complex(8) :: mvm, zt1, vc
    integer    :: ikfirst, iklast

    integer(4), allocatable :: idxpair(:,:)
    complex(8), allocatable :: minm(:,:,:)
    complex(8), allocatable :: evecsv(:,:)

    complex(8), external :: zdotc

    call cpu_time(tstart)

    !----------------------------------------
    ! Read KS eigenvalues from file EVALSV.OUT
    !----------------------------------------
    if (allocated(evalfv)) deallocate(evalfv)
    allocate(evalfv(nstfv,kset%nkpt))
    evalfv(:,:) = 0.d0
    do ik = 1, nkpt
      call getevalfv(kset%vkl(:,ik), evalfv(:,ik))
    end do

    ! VB / CB state index
    call find_vbm_cbm(1, nstfv, kset%nkpt, evalfv, efermi, nomax, numin, ikvbm, ikcbm, ikvcm)
    ! write(*,*) 'calc_vxnl: ', nomax, numin, efermi

    ! BZ integration weights
    call kintw()
    deallocate(evalfv)

    ! singular term prefactor
    sxs2 = 4.d0*pi*vi*singc2*kqset%nkpt

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

    !------------------------------------------!
    ! Matrix elements of non-local potential   !
    !------------------------------------------!
    if (allocated(vxnl)) deallocate(vxnl)
    allocate(vxnl(nstfv,nstfv,kset%nkpt))
    vxnl(:,:,:) = zzero

    ie12tot = 0.5d0*nstfv*(nstfv+1)
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
    if ((input%groundstate%outputlevelnumber>1) .and.rank==0) then
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
          call setbarcev(0.d0)

          if ((input%groundstate%outputlevelnumber>1) .and.rank==0) then
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
          call expand_evec(ik,'t')
          call expand_evec(jk,'c')

          !=================================
          ! Loop over m-blocks in M^i_{nm}
          !=================================
          do iblk = 1, nblk

            mstart = 1 + (iblk-1)*mblksiz
            mend = min(mdim, mstart+mblksiz-1)

            ! m-block M^i_{nm}
            allocate(minmmat(mbsiz,1:nstfv,mstart:mend))
            if ((input%groundstate%outputlevelnumber>1) .and.rank==0) then
              msize = sizeof(minmmat)*b2mb
              write(60,'(a,3i8,f14.2)') '    iblk, mstart, mend, size(minm) (Mb):', &
                                        iblk, mstart, mend, msize
            end if

            !---------------------
            ! Calculate M^i_{nm}
            !---------------------
            call expand_products(ik, iq, 1, nstfv, -1, mstart, mend, nomax, minmmat)

            ! sum over occupied states
            do ie3 = mstart, mend
              if (ie3 <= nomax) then
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ie12, ie1, ie2, mvm)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
                do ie12 = 1, ie12tot
                  ie1 = idxpair(1,ie12)
                  ie2 = idxpair(2,ie12)
                  mvm = zdotc(mbsiz, minmmat(:,ie1,ie3), 1, minmmat(:,ie2,ie3), 1)
                  vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp) - kiw(ie3,jk)*mvm
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
                  mvm = zdotc(mbsiz, minmmat(:,ie1,ie3), 1, minmmat(:,ie2,ie3), 1)
                  vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp) - ciw(ic,ias)*mvm
                end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
              end if ! core
            end do ! ie3

            deallocate(minmmat)

          end do ! iblk

          !--------------------------
          ! add singular term (q->0)
          !--------------------------
          if (Gamma) then
            do ie1 = 1, nomax
              vxnl(ie1,ie1,ikp) = vxnl(ie1,ie1,ikp) - sxs2*kiw(ie1,ik)
            end do
          end if

        end if ! rank

        call delete_coulomb_potential()

      end do ! iq

    end do ! ikp

    ! clear memory
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)
    deallocate(idxpair)

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
      do ie1 = 1, nomax
        exnl = exnl + kset%wkpt(ikp)*vxnl(ie1,ie1,ikp)
      end do
    end do ! ikp

    call cpu_time(tend)

    return
end subroutine
