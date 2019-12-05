
subroutine init_dft_eigenvalues()

    use modinput
    use modmain, only : nstfv, nmatmax, efermi, evalcr, &
                        occmax, chgval, filext
    use modgw
    use mod_mpi_gw, only : myrank
    use mod_hdf5
    implicit none

    integer :: ikp, ik, ib
    real(8) :: e0, egap

    if (allocated(evalfv)) deallocate(evalfv)
    allocate(evalfv(nstfv,kset%nkpt))
    evalfv(:,:) = 0.d0

    if (allocated(occfv)) deallocate(occfv)
    allocate(occfv(nstfv,kset%nkpt))
    occfv(:,:) = 0.d0

    !---------------------------------------------
    ! Read KS eigenvalues from file EVALFV_GW.OUT
    !---------------------------------------------
    do ikp = 1, kset%nkpt
      ik = kset%ikp2ik(ikp)
      if (ldapu == 0) then
        call getevalfv(kqset%vkl(:,ik), evalfv(:,ikp))
      else
        call getevalsv(kqset%vkl(:,ik), evalfv(:,ikp))
      end if
    end do

    !----------------------------------------
    ! find Fermi energy (LIBBZINT routine)
    !----------------------------------------
    call fermi_exciting(.false., &
                        chgval, &
                        nstfv, kset%nkpt, evalfv, &
                        kset%ntet, kset%tnodes, kset%wtet, kset%tvol, &
                        efermi, egap, fermidos)

    ! Calculate state occupation numbers
    call tetiw(kset%nkpt, kset%ntet, nstfv, evalfv, kset%tnodes, &
               kset%wtet, kset%tvol, efermi, occfv)
    do ik = 1, kset%nkpt
      do ib = 1, nstfv
        occfv(ib,ik) = 2.d0/kset%wkpt(ik)*occfv(ib,ik) ! prefactor 2 due to spin degeneracy in FV
      end do
    end do

    ! Setup the energy scale: Ef_KS = 0
    evalfv(:,:) = evalfv(:,:)-efermi
    evalcr(:,:) = evalcr(:,:)-efermi
    efermi = 0.d0

    !------------------------------------------------------------------

    nvelgw = chgval-2.d0*dble(ibgw-1)

    ! initialize the number of states to calculate the dielectric function
    nstdf = int(chgval/2.d0)+input%gw%nempty+1
    if (nstdf > nstfv) then
      nstdf = nstfv
      write(fgw,*)
      write(fgw,*)'WARNING(init_dft_eigenvalues) nstdf > nstfv !'
      write(fgw,*)
    end if

    ! initialize the number of states to calculate the correlation self energy
    if (input%gw%selfenergy%nempty>0) then
        nstse = int(chgval/2.d0)+input%gw%selfenergy%nempty+1
        if (nstse > nstfv) then
          nstse = nstfv
          write(fgw,*)
          write(fgw,*)'WARNING(init_dft_eigenvalues) nstse > nstfv !'
          write(fgw,*)
        end if
    else
        nstse = nstdf
    end if

    !----------------------------------------
    ! Output band structure summary
    !----------------------------------------
    if (myrank==0) then
      call boxmsg(fgw,'-',"Kohn-Sham eigenstates summary")
      write(fgw,*)'Maximum number of LAPW states:             ', nmatmax
      write(fgw,*)'Minimal number of LAPW states:             ', minval(nmat(1,:))
      write(fgw,*)'Number of states used in GW:'
      write(fgw,*)'    - total KS                             ', nstfv
      write(fgw,*)'    - occupied                             ', int(chgval/2.d0)
      write(fgw,*)'    - unoccupied                           ', input%gw%nempty
      write(fgw,*)'    - dielectric function                  ', nstdf
      write(fgw,*)'    - self energy                          ', nstse
      e0 = maxval(evalfv(nstfv,:))
      write(fgw,'(a,f12.6)')' Energy of the highest unoccupied state:    ', e0
      write(fgw,*)'Number of valence electrons:               ', int(chgval)
      write(fgw,*)'Number of valence electrons treated in GW: ', int(nvelgw)
      if (nstfv<=input%gw%nempty) then
        write(fgw,*)
        write(fgw,*)'WARNING(init_dft_eigenvalues) One uses the maximum number of available states!'
        write(fgw,*)
      end if
      call flushifc(fgw)
    end if

    !---------------------------------------------------------
    ! Search for the indices of VBM and CBM (nomax and numin)
    !---------------------------------------------------------
    call bandstructure_analysis('Kohn-Sham band structure', 1, nstfv, kset%nkpt, evalfv, efermi)

    !-----------------------------------------------------------------
    ! Check for consistency with specified QP bands range [ibgw,nbgw]
    !-----------------------------------------------------------------
    ! lower QP band index
    if ( (ibgw<1) .or. (ibgw>nstfv) ) ibgw = 1
    if (ibgw >= numin) then
        if (myrank==0) then
          write(*,*) "ERROR(init_dft_eigenvalues): Wrong QP bands interval!"
          write(*,*) "  ibgw = ", ibgw, " >= CBM = ", numin
        end if
        stop
    end if
    ! upper QP band index
    if ((nbgw < 1) .or. (nbgw > nstfv)) then
        ! use just a limited range of states where QP corrections are applied
        nbgw = nstfv
    end if
    if (nbgw <= nomax) then
        if (myrank==0) then
          write(*,*) "ERROR(init_dft_eigenvalues): Wrong QP bands interval!"
          write(*,*) "  nbgw = ", nbgw, " <= VBM = ", nomax
        end if
        stop
    endif
    nbandsgw = nbgw-ibgw+1

    !------------------------
    ! If symmetry is used
    !------------------------
    if (input%gw%reduceq) call sym_state_degeneracy

    return
end subroutine
