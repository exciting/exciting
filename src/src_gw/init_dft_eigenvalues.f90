
subroutine init_dft_eigenvalues()

    use modinput
    use modmain, only : nstsv, nmatmax, evalsv, efermi, evalcr, &
    &                   occmax, chgval, nspinor
    use modgw
    use mod_mpi_gw, only : myrank
    use mod_hdf5
    implicit none
    
    integer :: ik, ik0, n
    integer :: ib
    real(8) :: e0, e1, egap
    
    integer, allocatable :: idx(:)
    real(8), allocatable :: ev(:)
    
    if (allocated(evalsv)) deallocate(evalsv)
    allocate(evalsv(nstsv,kset%nkpt))
    evalsv(:,:) = 0.d0

    !----------------------------------------
    ! Read KS eigenvalues from file EVALSV.OUT
    !----------------------------------------
    do ik = 1, kset%nkpt
      ik0 = kset%ikp2ik(ik)
      !call getevalsvgw(ik0,evalsv(:,ik))
      call getevalsvgw_new('GW_EVALSV.OUT',ik0,kqset%vkl(:,ik0), &
      &                     nstsv,evalsv(1,ik))
    end do
    
    !----------------------------------------
    ! find Fermi energy (LIBBZINT routine)
    !----------------------------------------
    !n = min( int(chgval/2.d0)+10, nstsv)
    !call fermi( kset%nkpt, n, &
    !&           evalsv(1:n,:), &
    !&           kset%ntet, &
    !&           kset%tnodes, &
    !&           kset%wtet, &
    !&           kset%tvol, &
    !&           chgval,nspinor,efermi,egap)
    call fermi_exciting(input%groundstate%tevecsv, &
    &                   chgval, &
    &                   nstsv,nkpt,evalsv, &
    &                   ntet,tnodes,wtet,tvol, &
    &                   efermi,egap,fermidos)
    
    !evalsv(:,:) = evalsv(:,:)-efermi
    !evalcr(:,:) = evalcr(:,:)-efermi
    !efermi = 0.d0
    
    ! apply scissor shift
    !e0 = 0.0037 ! 0.1 eV
    !do ik = 1, kset%nkpt
    !  do ib = 1, nstsv
    !    if (evalsv(ib,ik)>efermi) evalsv(ib,ik) = evalsv(ib,ik)+e0
    !  end do
    !end do    

    !---------------------------------------------------------
    ! Search for the indices of VBM and CBM (nomax and numin)
    !---------------------------------------------------------
    call bandstructure_analysis('Kohn-Sham bandstructure analysis', &
    &  1,nstsv,kset%nkpt,evalsv(1:nstsv,:),efermi)
       
    !-----------------------------------------------------------------
    ! Check for consistency with specified QP bands range [ibgw,nbgw]
    !-----------------------------------------------------------------
    ! lower QP band index
    if ((ibgw<1) .or. (ibgw>nstsv)) ibgw = 1
    if (ibgw >= numin) then
        if (myrank==0) then
          write(*,*) "ERROR(init_dft_eigenvalues): Wrong QP bands interval is chosen!"
          write(*,*) "  ibgw = ", ibgw, " >= CBM = ", numin
        end if
        stop
    end if
  
    ! upper QP band index
    if ((nbgw < 1) .or. (nbgw > nstsv)) then
        ! use just a limited range of states where QP corrections are applied
        nbgw = nstsv
    end if
    if (nbgw <= nomax) then
        if (myrank==0) then
          write(*,*) "ERROR(init_dft_eigenvalues): Wrong QP bands interval is chosen!"
          write(*,*) "  nbgw = ", nbgw, " <= VBM = ", nomax
        end if
        stop
    endif

    !------------------------------------------------------------------

    nbandsgw = nbgw-ibgw+1
    nvelgw = chgval-occmax*dble(ibgw-1)
    
    ! initialize the number of states to calculate the dielectric function
    nstdf = int(chgval/2.d0)+input%gw%nempty+1
    if (nstdf > nstsv) then
      nstdf = nstsv
      write(fgw,*)
      write(fgw,*)'WARNING(init_dft_eigenvalues) nstdf > nstsv !'
      write(fgw,*)
    end if
    
    ! initialize the number of states to calculate the correlation self energy
    if (input%gw%selfenergy%nempty>0) then
        nstse = int(chgval/2.d0)+input%gw%selfenergy%nempty+1
        if (nstse > nstsv) then
          nstse = nstsv
          write(fgw,*)
          write(fgw,*)'WARNING(init_dft_eigenvalues) nstse > nstsv !'
          write(fgw,*)
        end if
    else
        nstse = nstdf
    end if

    !----------------------------------------
    ! Output band structure summary 
    !----------------------------------------
    if (myrank==0) then   
      write(fgw,*)'Maximum number of LAPW states:             ', nmatmax
      write(fgw,*)'Minimal number of LAPW states:             ', minval(nmat(1,:))
      write(fgw,*)'Number of states used in GW:'
      write(fgw,*)'    - total KS                             ', nstsv
      write(fgw,*)'    - occupied                             ', int(chgval/2.d0)
      write(fgw,*)'    - unoccupied                           ', input%gw%nempty
      write(fgw,*)'    - dielectric function                  ', nstdf
      write(fgw,*)'    - self energy                          ', nstse
      e0 = maxval(evalsv(nstsv,:))
      write(fgw,'(a,f12.6)')' Energy of the highest unoccupied state:    ', e0
      write(fgw,*)'Number of valence electrons:               ', int(chgval) 
      write(fgw,*)'Number of valence electrons treated in GW: ', int(nvelgw)
      if (nstsv<=input%gw%nempty) then
        write(fgw,*)
        write(fgw,*)'WARNING(init_dft_eigenvalues) One uses the maximum number of available states!'
        write(fgw,*)
      end if
      call flushifc(fgw)
    end if
   
    !------------------------
    ! If symmetry is used
    !------------------------
    if (input%gw%reduceq) call sym_state_degeneracy
    
    return
end subroutine

