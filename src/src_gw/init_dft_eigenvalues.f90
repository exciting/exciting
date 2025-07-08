
subroutine init_dft_eigenvalues()
    use constants,  only: real_zero
    use mod_charge_and_moment, only: chgval
    use mod_corestate, only: evalcr
    use mod_eigenvalue_occupancy, only: nstfv, efermi, occmax, fermidos 
    use mod_eigensystem, only: nmat, nmatmax
    use mod_bands, only: numin, nomax, nstdf, nstse, evalfv, occfv, bandstructure_analysis
    use mod_gw_degeneracies, only: initialize_degeneracy_module, absolute_tolerance_gw_degeneracy, &
                                   relative_tolerance_gw_degeneracy
    use mod_LDA_LU, only: ldapu
    use mod_hdf5
    use mod_misc, only: filext
    use modinput, only: input
    use modgw, only: kset, kqset, ibgw, nvelgw, fgw, nbandsgw, nbgw
    use modmpi, only: rank
    use precision,  only: i32, dp

    implicit none

    integer(i32) :: ikp, ik, ib
    real(dp) :: e0, egap
    logical :: enforce_degeneracy

    enforce_degeneracy = input%gw%enforceDegeneracy

    if (allocated(evalfv)) deallocate(evalfv)
    allocate(evalfv(nstfv,kset%nkpt), source = real_zero)

    if (allocated(occfv)) deallocate(occfv)
    allocate(occfv(nstfv,kset%nkpt), source = real_zero)

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
        occfv(ib,ik) = 2.0_dp/kset%wkpt(ik)*occfv(ib,ik) ! prefactor 2 due to spin degeneracy in FV
      end do
    end do

    ! Setup the energy scale: Ef_KS = 0
    evalfv(:,:) = evalfv(:,:) - efermi
    evalcr(:,:) = evalcr(:,:) - efermi
    efermi = real_zero

    !------------------------------------------------------------------------------------------------
    ! Computing degenerate subspaces, this is needed to enforce symmetry in the self-energy operator
    !------------------------------------------------------------------------------------------------
    
    call initialize_degeneracy_module(kset, nstfv, evalfv, occfv, enforce_degeneracy)

    !------------------------------------------------------------------

    nvelgw = chgval - 2.0_dp * real(ibgw-1, kind=dp)

    ! initialize the number of states to calculate the dielectric function
    nstdf = int(chgval/2.0_dp, kind=i32) + input%gw%nempty + 1
    if (nstdf > nstfv) then
      nstdf = nstfv
      if (rank==0) then
        write(fgw,*)
        write(fgw,*)'WARNING(init_dft_eigenvalues) nstdf > nstfv !'
        write(fgw,*)
      end if
    end if

    ! Checking for truncation of degenerate subspaces
    if( enforce_degeneracy ) call check_degenerate_subspaces(nstdf)

    ! initialize the number of states to calculate the correlation self energy
    if (input%gw%selfenergy%nempty>0) then
        nstse = int(chgval/2.0_dp, kind=i32) + input%gw%selfenergy%nempty + 1
        if (nstse > nstfv) then
          nstse = nstfv
          if (rank==0) then
            write(fgw,*)
            write(fgw,*)'WARNING(init_dft_eigenvalues) nstse > nstfv !'
            write(fgw,*)
          end if
        end if

        ! Again check for the truncation of degenerate subspaces
        if( enforce_degeneracy ) call check_degenerate_subspaces(nstse)

    else
        nstse = nstdf
    end if

    !----------------------------------------
    ! Output band structure summary
    !----------------------------------------
    if (rank==0) then
      call boxmsg(fgw,'-',"Kohn-Sham eigenstates summary")
      write(fgw,*)'Maximum number of LAPW states:             ', nmatmax
      write(fgw,*)'Minimal number of LAPW states:             ', minval(nmat(1,:))
      write(fgw,*)'Number of states used in GW:'
      write(fgw,*)'    - total KS                             ', nstfv
      write(fgw,*)'    - occupied                             ', int(chgval/2.0_dp, kind=i32)
      write(fgw,*)'    - unoccupied                           ', input%gw%nempty
      write(fgw,*)'    - dielectric function                  ', nstdf
      write(fgw,*)'    - self energy                          ', nstse
      e0 = maxval(evalfv(nstfv,:))
      write(fgw,'(a,f12.6)')' Energy of the highest unoccupied state:    ', e0
      write(fgw,*)'Number of valence electrons:               ', int(chgval, kind=i32)
      write(fgw,*)'Number of valence electrons treated in GW: ', int(nvelgw, kind=i32)
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
    call bandstructure_analysis('Kohn-Sham band structure', 1, evalfv, efermi, .true.)

    !-----------------------------------------------------------------
    ! Check for consistency with specified QP bands range [ibgw,nbgw]
    !-----------------------------------------------------------------
    ! lower QP band index
    if ( (ibgw<1) .or. (ibgw>nstfv) ) ibgw = 1
    if (ibgw >= numin) then
        if (rank==0) then
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
        if (rank==0) then
          write(*,*) "ERROR(init_dft_eigenvalues): Wrong QP bands interval!"
          write(*,*) "  nbgw = ", nbgw, " <= VBM = ", nomax
        end if
        stop
    end if
    nbandsgw = nbgw - ibgw + 1
    
    return

contains

    !> When we truncate the number of bands, there's a risk of cutting a degenerate
    !> subspace. In such a case the dielectric matrix and the self-energy are not ensured 
    !> to retain the full symmetry of the crystal, which can lead to the lift of the degeneracies.
    !> To prevent this symmetry breaking, we adjust the number of empty bands
    !> downward, a choice made for hybrid compatibility. For a more
    !> detailed discussion, see Section 6 of https://doi.org/10.1016/j.cpc.2011.12.006.
    subroutine check_degenerate_subspaces(ntruncation)

      use math_utils, only: get_degeneracies

      !> The truncation limit
      integer(i32), intent(inout) :: ntruncation
      ! Elements to check if we are truncating a degenerated subspace
      integer(i32) :: n_spaces, ispace
      !idx_degeneracies is allocated/deallocated inside the get_degeneracies procedure 
      integer(i32), allocatable :: idx_degeneracies(:,:)
      logical, allocatable      :: degeneracies(:,:)
      logical :: my_rank_writes_to_output


      my_rank_writes_to_output = ( rank == 0 )
      ! Here we check that we are not working with the full space (i.e. the maximum number
      ! of states available)
      if (ntruncation <= nstfv .and. minval(nmat(1,:)) /= nstfv) then
        if( my_rank_writes_to_output ) then
          write(fgw,*)
          call boxmsg(fgw,'-',"Checking for possible degenerate subspace truncation")
          write(fgw,*) "   -Initial band truncation was ", ntruncation
        end if

        ! Get a complete list of the degenerate states
        ! So degeneracies(ib,ik) is true if the ib-th state of ik-th kpoint
        ! is inside a degenerate subspace bigger in size that itself. 
        allocate(degeneracies(nstfv,kset%nkpt), source = .true.)
        do ikp = 1, kset%nkpt
          ! Get degeneracies and the number of degenerate subspaces
          ! Call get_degeneracies, which returns a list of the degenerate subspaces [each subspace info is
          ! giving in individual rows providing the starting [row 1] and ending [row 2]
          ! indices for each them, as well their size [row 3].
          idx_degeneracies = get_degeneracies(evalfv(:nstfv,ikp), absolute_tolerance_gw_degeneracy, &
                                              relative_tolerance_gw_degeneracy)
          n_spaces = size(idx_degeneracies, 2)

          ! Notice that we ignore the highest eigenvalue subspace as we cannot ensure that
          ! we haven't cut a degenerate subspace without computing a higher energy eigenvalue.
          ! Something that in the case of Hybrids is not possible.
          do ispace = n_spaces - 1, 1, -1
            if (idx_degeneracies(3, ispace) == 1) degeneracies(idx_degeneracies(1,ispace),ikp) = .false.
          end do
        end do

        ! Now check the first state in the list that is not degenerate
        ! and set the truncation to that value
        do ib = ntruncation, 1, -1
          if (.not. any(degeneracies(ib,:))) then
            ntruncation = ib
            exit
          end if
        end do
        ! Free space
        deallocate(degeneracies, idx_degeneracies)

        ! Print new truncation
        if( my_rank_writes_to_output ) then
          write(fgw,*) "   -Final band truncation has been adjusted to prevent the cutting of the degenerate subspaces to ", ntruncation
          write(fgw,*)
        end if
      else
        ! We are working with the full space, so no correction is needed.
        if( my_rank_writes_to_output ) then
          write(fgw,*)
          call boxmsg(fgw,'-',"Checking for possible degenerate subspace truncation")
          write(fgw,*) "   -Working with the maximum of our space, nothing to check."
          write(fgw,*)
        end if
      end if

    end subroutine check_degenerate_subspaces

end subroutine init_dft_eigenvalues
