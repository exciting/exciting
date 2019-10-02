!
!BOP
! !ROUTINE: calc_vxnl
! !INTERFACE:
!
subroutine calc_vxnl_test()
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
      call setbarcev(input%gw%barecoul%barcevtol)

      ! M^i_{nm}+M^i_{cm}
      allocate(minmmat(mbsiz,nstfv,mdim))

      !---------------------------------------
      ! Loop over k-points
      !---------------------------------------
      ikq = 0

      do ikp = 1, kset%nkpt

          ikq = ikq + 1

          if (mod(ikq, procs) == rank) then

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

            call expand_products(ik, iq, 1, nstfv, -1, 1, mdim, nomax, minmmat)

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
              if ( Gamma .and. (ie1<=nomax) ) then
                vxnl(ie1,ie1,ikp) = vxnl(ie1,ie1,ikp) - sxs2*kiw(ie1,ik)
              end if
            end do ! ie1

          end if ! rank

        end do ! ikp

        deallocate(minmmat)
        call delete_coulomb_potential()

    end do ! iq

    ! clear memory
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)

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
