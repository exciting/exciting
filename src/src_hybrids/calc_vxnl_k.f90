
subroutine calc_vxnl_k(ikp, COMM_LEVEL_2)

    use modinput
    use modmain,               only: pi, nmatmax, nstfv, apwordmax, lmmaxapw, natmtot, zzero, idxas
    use mod_hybrids,           only: vxnl, kset, kqset, nomax, ncg, vi, singc2, matsiz, Gamma, &
                                     mbsiz, locmatsiz, Gqset, Gkqset, eveck, eveckp, kiw, ciw, &
                                     eveckalm, eveckpalm, minmmat, corind
    use mod_coulomb_potential, only: delete_coulomb_potential
    use mod_misc_gw,           only: gammapoint
    use modmpi
    implicit none

    ! input/output
    integer(4), intent(in) :: ikp
    integer(4), intent(in) :: COMM_LEVEL_2

    ! local
    integer(4) :: ik, jk, iq
    integer(4) :: ie1, ie2, ie3, icg
    integer(4) :: is, ia, ias, ic
    integer(4) :: n, nmdim, m, mdim
    real(8)    :: sxs2
    complex(8) :: mvm, zt1
    integer(4) :: level2rank, level2procs
    complex(8), external :: zdotc

    vxnl(1:nstfv,1:nstfv,ikp) = zzero

    ! singular term prefactor
    sxs2 = 4.d0*pi*vi*singc2*kqset%nkpt

    if ((input%gw%coreflag=='all').or. &
    &   (input%gw%coreflag=='xal')) then
      mdim = nomax+ncg
    else
      mdim = nomax
    end if

    !---------------------------------------
    ! Integration over BZ
    !---------------------------------------

#ifdef MPI
    call MPI_COMM_SIZE(COMM_LEVEL_2, level2procs, ierr)
    call MPI_COMM_RANK(COMM_LEVEL_2, level2rank, ierr)
    call MPI_BARRIER(COMM_LEVEL_2, ierr)
#else
    level2rank = 0
    level2procs = 1
#endif

    do iq = 1, kqset%nkpt

      if ( (mod(iq-1, level2procs) == level2rank) .and. (level2rank < kqset%nkpt) ) then

        write(*,*) "for ikp ", ikp, " do iq ", iq, " on proc", rank

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

        ! k-q/k eigenvectors
        allocate(eveck(nmatmax,nstfv))
        allocate(eveckp(nmatmax,nstfv))
        call getevecfv(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), eveck)
        eveckp = conjg(eveck)
        call getevecfv(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), eveck)
        allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot))
        allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot))
        call expand_evec(ik,'t')
        call expand_evec(jk,'c')

        ! M^i_{nm}+M^i_{cm}
        allocate(minmmat(mbsiz,nstfv,1:mdim))
        minmmat(:,:,:) = zzero
        call expand_products(ik, iq, 1, nstfv, -1, 1, mdim, nomax, minmmat)
        call delete_coulomb_potential()

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

        ! clear memory
        deallocate(minmmat)
        deallocate(eveck)
        deallocate(eveckp)
        deallocate(eveckalm)
        deallocate(eveckpalm)

      end if

    end do ! iq

#ifdef MPI
    if (level2rank >= kqset%nkpt) vxnl(1:nstfv,1:nstfv,ikp) = zzero
    call MPI_ALLREDUCE(MPI_IN_PLACE, vxnl(1:nstfv,1:nstfv,ikp), nstfv*nstfv, MPI_DOUBLE_COMPLEX,  MPI_SUM, &
                       COMM_LEVEL_2, ierr)
#endif

    ! complete the matrix employing the matrix hermiticity
    do ie1 = 1, nstfv
      do ie2 = ie1+1, nstfv
        vxnl(ie2,ie1,ikp) = conjg(vxnl(ie1,ie2,ikp))
      end do
    end do

end subroutine