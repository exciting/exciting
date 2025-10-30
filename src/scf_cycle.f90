subroutine scf_cycle(verbosity)
    use cdft, only: cdft_input_keys, deallocate_cdft_global_arrays, determine_cdft_occupations, &
      file_extension_GS, initialize_cdft_global_arrays, update_occupations_with_the_maximum_overlap_method
    use exciting_mpi, only: xmpi_bcast, xmpi_allreduce, xmpi_allgatherv
    use lo_recommendation, only: recommend_local_orbital_trial_energies
    use mod_APW_LO, only: apwn, apwe0, lorbe0, lorbl, lorbord, lorbn, maxapword, maxlapw, nlorb
    use mod_atoms, only: atposc, idxas, natoms, natmtot, nspecies, spr, spsymb
    use mod_charge_and_moment, only: chgdst, chgpart, momtot
    use mod_convergence, only: currentconvergence, vchgdst, vcurrentconvergence, vdeltae
    use mod_eigenvalue_occupancy, only: evalsv, fermidos, occsv, nstfv, nstsv
    use mod_eigensystem, only: mt_hscf, MTInitAll, MTNullify, nmatmax
    use mod_energy, only: engytot, engyknst
    use mod_force, only: forcemax, forcetot
    use mod_Gvector, only: ngrid, ngrtot
    use mod_Gkvector, only: vgkl
    use mod_kpoint, only: nkpt, vkl, wkpt
    use mod_LDA_LU, only: ldapu, lmmaxlu
    use mod_misc, only: filext, task, tlast, tstop
    use mod_muffin_tin, only: lmmaxvr, nrmt, nrmtmax
    use mod_OEP_HF, only: resoep
    use mod_potential_and_density, only: generate_density_and_magnetization, m2effig, magir, magmt, meffig, rhomt, rhoir, veffmt, veffir, vxcmt, vxcir, exmt, exir, ecmt, ecir, xctype
    use mod_spin, only: ndmag, nspinor, nspnfv
    use mod_timing, only: stopwatch, time_density_init, time_pot_init, timefor, timefv, &
      timeinit, timeio, timemat, timemixer, timemt, timepot, timerho, timesv
    use modinput, only: input, getfixspinnumber
    use modmpi, only: barrier, firstofset, lastofset, mpiglobal, &
      procs, rank, splittfile
    use precision, only: dp, i32
    use scl_xml_out_Module, only: deltae, dforcemax, iscl, scl_iter_xmlout, scl_xml_out_write, scl_xml_write_moments
    use secular_equation, only: seceqn
    use sirius_api,    only: set_radial_functions_sirius, solve_seceqn_sirius, get_eval_sirius, get_evec_sirius, &
                             put_occ_sirius, generate_density_sirius, get_periodic_function_sirius
    use sirius_init,   only: sirius_options
    use total_energy, only: energy
    use to_char_conversion, only: to_char
    use trial_energy_selection, only: select_apw_trial_energies, select_local_orbital_trial_energies
    use TS_vdW_module, only: C6ab, R0_eff_ab
    use mod_gen_lo, only: genlofr
    use kinetic_energy_density, only: gen_ked, ked_mt, ked_cr, ked_ir, ked_magmt, ked_magir
    use kinetic_energy_density_vars, only: ked_var_init, ked_var_free, timeked
    use mgga_potxc
    use mgga_poteff
    use mgga_init
    use mGGA_eigensystem
    
    Implicit None

    integer(i32), intent(IN) :: verbosity
    Real(dp) :: et, fm
    Real(dp), Allocatable :: evalfv(:, :)
    Complex(dp), Allocatable :: evecfv(:, :, :)
    Complex(dp), Allocatable :: evecsv(:, :)
    Logical :: exist
    Integer(i32) :: ik, is, ia, idm, id, lmax, nodesmax
    Integer(i32) :: n, nwork
    Real(dp), Allocatable :: v(:)
    Real(dp) :: timetot, ts0, ts1, tin1, tin0, ta,tb
    
    character(len=77) :: string, acoord

    Real(dp), Allocatable :: rhomtref(:,:,:) ! muffin-tin charge density (reference)
    Real(dp), Allocatable :: rhoirref(:)     ! interstitial real-space charge density (reference)

    type(cdft_input_keys) :: cdft_calculation
    logical :: spin_polarization, use_mGGA
    integer(i32) :: first_k, last_k
    complex (dp), allocatable :: evecfv_store(:, :, :)
    real (dp), allocatable :: occsv_ref(:, :)

    first_k = firstofset(rank, nkpt)
    last_k = lastofset(rank, nkpt)
    acoord = "lattice"
    spin_polarization = associated( input%groundstate%spin )
    if (input%structure%cartesian) acoord = "cartesian"

    If ((verbosity>-1).and.(rank==0)) Then
        write(string,'("Self-consistent loop started")')
        call printbox(60,"+",string)
    End If

! reset density-dependent dispersion coefficients of Tkatchenko-Scheffler method
    If (Allocated(C6ab)) Deallocate(C6ab)
    If (Allocated(R0_eff_ab)) Deallocate(R0_eff_ab)
!
    call MTNullify(mt_hscf)

!_______________________________________________________________
! initialise or read the charge density and potentials from file
    call stopwatch("exciting:init_scf", 1)
!! TIME - Begin of initialisation segment
    Call timesec (ts0)
    If ((task == 1) .or. (task == 3)) Then
        Call readstate
        If ((verbosity>-1).and.(rank==0)) write(60,'(" Potential read in from STATE.OUT")')
    Else If (task == 7) Then
        ! restart from previous HYBRID iteration
        call genmeffig()

       continue
    Else If (task == 200) Then
        Call phveff
        If ((verbosity>-1).and.(rank==0)) write(60,'(" Supercell potential constructed from STATE.OUT")')
    Else
        Call timesec(tin0)
        Call rhoinit
        Call timesec(tin1)
        time_density_init=tin1-tin0
        call timesec(tin0)
        if ( associated(input%groundstate%mgga) ) then 
            call init_mgga()
            call calc_poteff_gga(veffmt, veffir, 2, xctype, rhomt, rhoir, vxcmt, vxcir, & 
                                 exmt, ecmt, ecir, exir)                
        else 
            call poteff( .true. )
        end if 
        Call genveffig
        Call timesec(tin1)
        time_pot_init=tin1-tin0
        If ((verbosity>-1).and.(rank==0)) write(60,'(" Density and potential initialised from atomic data")')
    End If

    Call genmeffig
    If ((verbosity>-1).and.(rank==0)) then
        write (60, *)
        Call flushifc (60)
    end if
    call stopwatch("exciting:init_scf", 0)

    call MTNullify(mt_hscf)

!_____________________________
! convergence vectors

    If (allocated(vcurrentconvergence)) deallocate(vcurrentconvergence)
    Allocate(vcurrentconvergence(input%groundstate%niterconvcheck))
    vcurrentconvergence=0.d0
    If (allocated(vdeltae)) deallocate(vdeltae)
    Allocate(vdeltae(input%groundstate%niterconvcheck))
    vdeltae=0.d0
    If (allocated(vchgdst)) deallocate(vchgdst)
    Allocate(vchgdst(input%groundstate%niterconvcheck))
    vchgdst=0.d0

!_____________________________
! reference density

    If (allocated(rhomtref)) deallocate(rhomtref)
    Allocate(rhomtref(lmmaxvr,nrmtmax,natmtot))
    rhomtref(:,:,:) = rhomt(:,:,:)
    If (allocated(rhoirref)) deallocate(rhoirref)
    Allocate (rhoirref(ngrtot))
    rhoirref(:) = rhoir(:)

    Call timesec (ts1)
    timeinit = timeinit+ts1-ts0
!! TIME - End of initialisation segment

!-----------------------------------------------------
! CDFT
    call cdft_calculation%read_input_keys( input%groundstate )
    if( cdft_calculation%is_on() ) then
      string = filext
      filext = file_extension_GS
      do ik = 1, nkpt
        call getoccsv( vkl(:, ik), occsv(:, ik) )
      end do
      if ( cdft_calculation%read_density_potential_from_file() ) call readstate
      call determine_cdft_occupations( cdft_calculation, wkpt, occsv )
      if ( cdft_calculation%is_maximum_overlap_method_required() ) then
        occsv_ref = occsv
        allocate( evecfv_store(nmatmax, nstfv, first_k:last_k) )
        splittfile = .false.
        do ik = first_k, last_k
          call getevecfv( vkl(:, ik), vgkl(:, :, :, ik), evecfv_store(:, :, ik) )
        end do
        call initialize_cdft_global_arrays( evecfv_store, first_k )
      end if
      filext = string
    end if

!----------------------------------------------------
!! TIME - Mixer segment
    if ( associated(input%groundstate%mgga) ) then
        iscl = 0 
    else 
        Call timesec (ts0)
        ! size of mixing vector
        n = lmmaxvr*nrmtmax*natmtot+ngrtot
        If (spin_polarization) n = n*(1+ndmag)
        If (ldapu .Ne. 0) n = n + 2*lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot
        ! allocate mixing arrays
        Allocate (v(n))
        ! call mixing array allocation functions by setting
        nwork = -1
        ! and call interface
        iscl = 0
        Call packeff (.True., n, v)
        If (rank .Eq. 0) Call mixerifc(input%groundstate%mixernumber, n, v, currentconvergence, nwork)
        Call packeff (.False., n, v)
        Call timesec (ts1)
        timemixer = ts1-ts0+timemixer
    end if 
!! TIME - End of mixer segment
!----------------------------------------------------

! set last iteration flag
    tlast = .False.
! set stop flag
    tstop = .False.
    engytot = 0.d0
    fm = 0.d0
! delete any existing eigenvector files
    If ((rank .Eq. 0) .And. ((task .Eq. 0) .Or. (task .Eq. 2))) Call delevec()

!! TIME - First IO segment
!----------------------------------------!
! begin the self-consistent loop
!----------------------------------------!
    call stopwatch("exciting:scf", 1)
    Do iscl = 1, input%groundstate%maxscl
        call timesec (ts0)
! exit self-consistent loop if last iteration is complete
        if (tlast) then
            If ((verbosity>-1).and.(rank==0)) Then
                call printline(60," ")
                call printline(60,"+")
                If (input%groundstate%niterconvcheck.ge.2) Then
                   write(string,'("Convergency criteria checked for the last", I2," iterations ")') input%groundstate%niterconvcheck
                Else
                   write(string,'("Convergency criteria checked for the last iteration ")')
                End If
                call printtext(60,"+",string)
                write(string,'("Convergence targets achieved. Performing final SCF iteration")')
                call printtext(60,"+",string)
                call printline(60,"+")
                Call flushifc(60)
            End If
        else
            If (iscl .Ge. input%groundstate%maxscl) Then
                If ((verbosity>-1).and.(rank==0)) Then
                    write(string,'("Reached self-consistent loops maximum : ", I4)') &
                   &  input%groundstate%maxscl
                    call printbox(60,"+",string)
                    call warning('Warning(gndstate): Reached self-consistent loops maximum')
                    Call flushifc(60)
                End If
                tlast = .True.
                goto 10
            End If
            If ((verbosity>-1).and.(rank==0)) Then
                write(string,'("SCF iteration number : ", I4)') iscl
                call printbox(60,"+",string)
                Call flushifc(60)
            End If
        end if

10      continue

        call timesec (ts1)
        timeio = ts1 - ts0 + timeio

!! TIME - End of first IO segment

!! TIME - Muffin-tin segment
        Call timesec (ts0)

        if (task /= 7) then
          ! No updates of core and valence radial functions during hybrids run
          call gencore          ! generate the core wavefunctions and densities
          ! find the first linearization energies 
          if (iscl==1) then 
             call select_local_orbital_trial_energies(nlorb, lorbord, lorbl, lorbn, nspecies, idxas, nrmt, spr, veffmt(1,:,:), lorbe0) 
             call select_apw_trial_energies(maxapword, maxlapw, apwn(:, 0:, :), nspecies, idxas, nrmt, spr, veffmt(1,:,:), apwe0(:, 0:, :))
          endif
          call linengy          ! find the new linearization energies
          if (rank==0) call writelinen
          call genapwfr         ! generate the APW radial functions
          call genlofr          ! generate the local-orbital radial functions
          ! compute recommendations for local orbital trial energies
          if ((associated(input%groundstate%lorecommendation)) .and. (tlast)) then
            nodesmax = input%groundstate%lorecommendation%nodesmaxlo
            lmax = input%groundstate%lorecommendation%lmaxlo
            call recommend_local_orbital_trial_energies(nodesmax, lmax, nspecies, spsymb, idxas, nrmt, spr, veffmt(1,:,:))
          endif
          call olprad           ! compute the overlap radial integrals
        end if
        if ( associated(input%groundstate%sirius) ) then
          call set_radial_functions_sirius( input )
        end if
!------------------------------------------------------------
! Effective Hamiltonian Setup: Radial and Angular integrals
!------------------------------------------------------------
        call stopwatch("exciting:rad_int", 1)
        if ( associated(input%groundstate%mgga) ) then
            use_mGGA = ( ((task == 1) .or. (task == 3)) .and. mgga_read_in ) .or. (iscl >= 2)
            if (use_mGGA) then
                call mGGA_eig_init(veffmt, veffir, vxcmt_mgga_nonmult, vxcir_mgga_nonmult)
            else
                call mGGA_eig_init(veffmt, veffir)
            end if
        else 
            call MTInitAll(mt_hscf) 
            call hmlint(mt_hscf)
        end if 
        call stopwatch("exciting:rad_int", 0)
!________________
! partial charges

        if (input%groundstate%tpartcharges) then
            allocate(chgpart(lmmaxvr,natmtot,nstsv))
            chgpart(:,:,:)=0.d0
        end if

        call timesec (ts1)
        timemt = ts1 - ts0 + timemt

!! TIME - End of muffin-tin segment

!! TIME - Second IO segment

        call timesec (ts0)

!-----------------------------------------------
! Solve Secular Equation for each k-point
!-----------------------------------------------
        call timesec(ta)

        if ( associated(input%groundstate%sirius) .and. sirius_options%use_eigen_states ) then
          call solve_seceqn_sirius()
          call get_eval_sirius()
          call get_evec_sirius()
        else
! start k-point loop
#ifdef MPI
            call barrier()
            If (rank == 0) then 
              if( input%groundstate%solver%type /= 'Davidson' .or. procs > 1 ) Call delevec()
            end if
            splittfile = .True.
#else
            splittfile = .False.
#endif
            Do ik = first_k, last_k

!____________________________________________
! every thread should allocate its own arrays

                Allocate (evalfv(nstfv, nspnfv))
                Allocate (evecfv(nmatmax, nstfv, nspnfv))
                Allocate (evecsv(nstsv, nstsv))

                if (iscl.le.1) then
                  evecfv=0d0
                elseif (input%groundstate%solver%type.eq.'Davidson') then
                  Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
                endif
!! TIME - seceqn does not belong to IO

                call timesec(ts1)
                timeio = ts1 - ts0 + timeio

!__________________________________________________________
! solve the first- and second-variational secular equations
                call seceqn (ik, evalfv, evecfv, evecsv, cdft_calculation%is_maximum_overlap_method_required() )

                call timesec(ts0)

!______________________________________
! write the eigenvalues/vectors to file

                Call putevalfv (ik, evalfv)
                Call putevalsv (ik, evalsv(:, ik))
                Call putevecfv (ik, evecfv)
                Call putevecsv (ik, evecsv)

!__________________________
! store evecfv for the case of maximum overlap in cdft
                if ( cdft_calculation%is_maximum_overlap_method_required() &
                    .and. (ik >= first_k) .and. (ik<=last_k) ) &
                    evecfv_store(1:nmatmax, 1:nstfv, ik) = evecfv(1:nmatmax, 1:nstfv, 1)

!__________________________
! calculate partial charges
                if (input%groundstate%tpartcharges) call genpchgs(ik,evecfv,evecsv)
                deallocate (evalfv, evecfv, evecsv)

            End Do ! ik

! end k-point loop -------------------------------------------------------------
            call xmpi_allgatherv( mpiglobal, evalsv, nstsv * (last_k - first_k + 1) )
            if ( task == 7 ) call xmpi_allgatherv( mpiglobal, engyknst, nstfv * (last_k - first_k + 1) )
        end if

        call timesec(tb)
! Release memory used by the MT Hamiltonian
!        call mt_hscf%release()

        call stopwatch("exciting:occ", 1)
!-----------------------------------------------
! find the occupation numbers and Fermi energy
!-----------------------------------------------
        if( cdft_calculation%is_on() ) then
          if ( cdft_calculation%is_maximum_overlap_method_required() ) then
            occsv = occsv_ref
            call update_occupations_with_the_maximum_overlap_method( evecfv_store, occsv(:, first_k:) )
          end if
        else 
          call occupy
        end if 

        If (rank==0) Then
! write out the eigenvalues and occupation numbers
            Call writeeval
! write the Fermi energy to file
            Call writefermi
        End If
!write the occupancies to file
        Do ik = first_k, last_k
            Call putoccsv (ik, occsv(:, ik))
        End Do
        if ( associated(input%groundstate%sirius) ) then
          call put_occ_sirius()
        end if
        call timesec(ts1)
        timeio = ts1 - ts0 + timeio
        call stopwatch("exciting:occ", 0)

        call stopwatch("exciting:rhomag", 1)
!-----------------------------------------------
! Calculate density and magnetization
!-----------------------------------------------
        if (associated(input%groundstate%sirius) .and. sirius_options%use_density ) then
          call generate_density_sirius()
          call get_periodic_function_sirius(rhoir, ngrid)
          call timesec(ts0)
        else
          call generate_density_and_magnetization
          call timesec(ts0)
#ifdef MPI
        ! EXX case
          If (input%groundstate%xctypenumber.Lt.0) Call mpiresumeevecfiles()
#endif

          if ((input%groundstate%tpartcharges).and.(rank==0)) then
              ! write out partial charges
              call writepchgs(69,input%groundstate%lmaxvr)
              call flushifc(69)
          end if
          call timesec(ts1)
          timeio = ts1 - ts0 + timeio
!! TIME - End of second IO segment

          call timesec(ts0)
! symmetrise the density
          call symrf(input%groundstate%lradstep, rhomt, rhoir)
! symmetrise the magnetisation
          If (spin_polarization) Call symrvf(input%groundstate%lradstep, magmt, magir)
! convert the density from a coarse to a fine radial mesh
          call rfmtctof (rhomt)
! convert the magnetisation from a coarse to a fine radial mesh
          Do idm = 1, ndmag
              Call rfmtctof (magmt(:, :, :, idm))
          End Do
        end if

! add the core density to the total density
        Call addrhocr
! calculate the charges
        Call charge( 's.c.f. loop iteration ' // to_char( iscl ) )
! calculate the moments
        If (spin_polarization) Call moment
! normalise the density
        Call rhonorm
        call stopwatch("exciting:rhomag", 0)
! LDA+U
        If (ldapu .Ne. 0) Then
! generate the LDA+U density matrix
            Call gendmatlu
! generate the LDA+U potential matrix
            Call genvmatlu
! write the LDA+U matrices to file
            if (rank .eq. 0) Call writeldapu
        End If
! generate charge distance
        call chgdist(rhomtref,rhoirref)
        do id=1, input%groundstate%niterconvcheck-1
           vchgdst(id) = vchgdst(id+1)
        end do
        vchgdst(input%groundstate%niterconvcheck) = chgdst

! store density to reference
        rhoirref(:)=rhoir(:)
        rhomtref(:,:,:)=rhomt(:,:,:)

! compute kinetic energy density 
        if ( associated(input%groundstate%mgga)  ) then 
            call timesec (ts0)
            call gen_ked()
 
            ! symmetrise the kinetic energy density
            call symrf(input%groundstate%lradstep, ked_mt, ked_ir)
            ! convert the density from a coarse to a fine radial mesh
            call rfmtctof (ked_mt)
            if (associated(input%groundstate%spin)) Call symrvf(input%groundstate%lradstep, ked_magmt, ked_magir)
            call timesec (ts1)
            timeked = timeked + ts1-ts0
        end if 

!-----------------------------------
! Compute the effective potential
!-----------------------------------
        call timesec (ts0)
        if ( associated(input%groundstate%mgga)) then 
            call calc_poteff_mgga(veffmt, veffir, 3, xctype_mgga, rhomt, rhoir, exmt, ecmt, ecir, exir, &
                                vxcmt, vxcmt_mgga_nonmult, vxcir, vxcir_mgga_nonmult, ked_ir, ked_mt)
            call calc_poteff_gga(veffmt_gga, veffir_gga, 2, xctype, rhomt, rhoir, vxcmt_gga, vxcir_gga, &
                                exmt_gga, ecmt_gga, ecir_gga, exir_gga)
        else 
            call poteff( .true. )
        end if 
        call timesec (ts1)
        timepot = ts1-ts0+timepot

!---------------
! Mixing
!---------------
        if ( associated(input%groundstate%mgga) ) then  
            call mgga_mixer(iscl, v, nwork, currentconvergence, vcurrentconvergence)
        else 
            Call timesec (ts1)
            ! pack interstitial and muffin-tin effective potential and field into one array
            Call packeff (.True., n, v)
            ! mix in the old potential and field with the new
            If (rank .Eq. 0) Then
                Call mixerifc (input%groundstate%mixernumber, n, v, currentconvergence, nwork)
                do id=1, input%groundstate%niterconvcheck-1
                    vcurrentconvergence(id) = vcurrentconvergence(id+1)
                end do
                vcurrentconvergence(input%groundstate%niterconvcheck) = currentconvergence
            End If
        call xmpi_bcast(mpiglobal, v)
            ! unpack potential and field
            Call packeff (.False., n, v)
        end if 
!---------------
! Fourier transform effective potential to G-space
        Call genveffig
        if (allocated(meffig)) deallocate(meffig)
        if (allocated(m2effig)) deallocate(m2effig)
! add the fixed spin moment effect field
        If (getfixspinnumber() .Ne. 0) Call fsmfield
        Call genmeffig
! reduce the external magnetic fields if required
        If (spin_polarization) Then
            If (input%groundstate%spin%reducebf .Lt. 1.d0) Then
                input%groundstate%spin%bfieldc(:) = &
               &  input%groundstate%spin%bfieldc(:) * input%groundstate%spin%reducebf
                Do is = 1, nspecies
                    Do ia = 1, natoms (is)
                        input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:) = &
                       &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:) * input%groundstate%spin%reducebf
                    End Do
                End Do
            End If
        End If

!--------------------------------
! compute the energy components
!--------------------------------
        et = engytot
        Call energy
        Call timesec(ts1)
        timepot=ts1-ts0+timepot

!----------------------------------------------
! compute the forces (without IBS corrections)
!----------------------------------------------
        if (input%groundstate%tforce) then
            Call force(.false.)
        end if

!-----------------------------
! Print results
!-----------------------------

!! TIME - Third IO segment
        call timesec(ts0)
        deltae=abs(et-engytot)
        do id=1, input%groundstate%niterconvcheck-1
           vdeltae(id)=vdeltae(id+1)
        end do
        vdeltae(input%groundstate%niterconvcheck)=abs(et-engytot)

        If ((verbosity>-1).and.(rank==0)) Then
! output energy components
            call writeengy(60)
            if (verbosity>0) Write (60,*)
            Write (60, '(" DOS at Fermi energy (states/Ha/cell)",T45, ": ", F18.8)') fermidos
! write DOS at Fermi energy to FERMIDOS.OUT and flush
!            Write (62, '(G18.10)') fermidos
!            Call flushifc (62)
! output charges and moments
            Call writechg (60,input%groundstate%outputlevelnumber)
! write total moment to MOMENT.OUT and flush
            If (spin_polarization) Then
                Write (63, '(3G18.10)') momtot (1:ndmag)
                Call flushifc (63)
            End If
! output effective fields for fixed spin moment calculations
            If (getfixspinnumber() .Ne. 0) Call writefsm (60)
! output forces to INFO.OUT
!            if (input%groundstate%tforce) call writeforce(60,input%relax%outputlevelnumber)
! write band-gap if the dos at the Fermi energy is smaller than the given threshold
            if ( .not. cdft_calculation%is_on() ) then
              if ( fermidos < 1.0d-4 ) call printbandgap(60)
            end if
! check for WRITE file
            Inquire (File='WRITE', Exist=exist)
            If (exist) Then
                Write (60,*)
                Write (60, '(" WRITE file exists - writing STATE.OUT")')
                Call writestate
                Open (50, File='WRITE')
                Close (50, Status='DELETE')
            End If
            Call scl_iter_xmlout ()
            If (spin_polarization) Call scl_xml_write_moments()
            Call scl_xml_out_write()
        End If
! write STATE.OUT file if required
        If (input%groundstate%nwrite .Ge. 1) Then
            If (Mod(iscl, input%groundstate%nwrite) .Eq. 0) Then
                Call writestate
                if ((verbosity>-1).and.(rank==0)) Then
                    write(60,*)
                    write(60, '(" Wrote STATE.OUT")')
                end if
            End If
        End If
        call timesec(ts1)
        timeio = ts1 - ts0 + timeio
!! TIME - End of third IO segment

! exit self-consistent loop if last iteration is complete
        If (tlast) goto 20

! update convergence criteria
        et = engytot
        if (input%groundstate%tforce) then
            dforcemax=abs(fm-forcemax)
            if (dforcemax .lt. 1.d-10) dforcemax=0.
        end if

!! TIME - Fourth IO segment
        call timesec(ts0)

! output the current total time
        timetot = timeinit+timemat+timefv+timesv+timerho+timepot+timefor+timeio+timemt+timemixer
        if ( associated(input%groundstate%mgga) ) timetot = timetot + timeked
        if ((verbosity>-1).and.(rank==0)) then
            write(60,*)
            write(60, '(" Wall time (seconds)",T45, ": ", F12.2)') timetot
        end if

! write TOTENERGY.OUT
        if ((verbosity>-1).and.(rank==0)) then
            Write (61, '(G22.12)') engytot
            Call flushifc (61)
        end if

!----------------------
! Convergence tests
!----------------------
        If (iscl .Ge. 2) Then

!...write convergence if only energy
            if ((verbosity>-1).and.(rank==0).and.(input%groundstate%scfconv.eq.'energy')) then
                write(60,*)
                write(60,'(" Absolute change in total energy   (target) : ",G13.6,"  (",G13.6,")")') &
                &     deltae, input%groundstate%epsengy
            end if

!...write convergence if only potential
            if ((verbosity>-1).and.(rank==0).and.(input%groundstate%scfconv.eq.'potential')) then
                if (associated(input%groundstate%OEP)) then
                    write(60,*)
                    write(60, '(" Magnitude of OEP residual",T45 ": ", F18.8)') resoep
                end if
                write(60,*)
                Write(60,'(" RMS change in effective potential (target) : ",G13.6,"  (",G13.6,")")') &
                &     currentconvergence, input%groundstate%epspot
            end if

!...write convergence if only energy
            if ((verbosity>-1).and.(rank==0).and.(input%groundstate%scfconv.eq.'charge')) then
                write(60,*)
                write(60,'(" Charge distance                   (target) : ",G13.6,"  (",G13.6,")")') &
                &     chgdst, input%groundstate%epschg
            end if

!...write convergence if multiple convergence
            if ((verbosity>-1).and.(rank==0).and.(input%groundstate%scfconv.eq.'multiple')) then
                write(60,*)
                Write(60,'(" RMS change in effective potential (target) : ",G13.6,"  (",G13.6,")")') &
               &    currentconvergence, input%groundstate%epspot
                write(60,'(" Absolute change in total energy   (target) : ",G13.6,"  (",G13.6,")")') &
               &    deltae, input%groundstate%epsengy
                write(60,'(" Charge distance                   (target) : ",G13.6,"  (",G13.6,")")') &
               &    chgdst, input%groundstate%epschg
            end if

!...if tforce=true write also convergence of non-IBS forces
            if ((verbosity>-1).and.(rank==0).and.(input%groundstate%tforce)) then
                write(60,'(" Abs. change in max-nonIBS-force   (target) : ",G13.6,"  (",G13.6,")")') &
               &    dforcemax, input%groundstate%epsforcescf
            end if

!...write in RMSDVEFF.OUT and DFSCFMAX.OUT
            if ((verbosity>-2).and.(rank==0)) then
                Write (65, '(G18.10)') currentconvergence
                Call flushifc(65)
                if (input%groundstate%tforce) then
                    Write (67, '(G22.12)') dforcemax
                    Call flushifc(67)
                end if
            end if

!-----------------------
! check for convergence
!-----------------------

            if (input%groundstate%scfconv .eq. 'energy') tlast = all(abs(vdeltae) .lt. input%groundstate%epsengy)

            if (input%groundstate%scfconv .eq. 'potential') tlast = all(abs(vcurrentconvergence) .lt. input%groundstate%epspot)

            if (input%groundstate%scfconv .eq. 'charge') tlast = all(abs(vchgdst) .lt. input%groundstate%epschg)

            if (input%groundstate%scfconv .eq. 'multiple') then
                tlast = all(abs(vcurrentconvergence) .lt. input%groundstate%epspot) .and. &
                &       all(abs(vdeltae) .lt. input%groundstate%epsengy) .and. &
                &       all(abs(vchgdst) .lt. input%groundstate%epschg)
            end if

            if (input%groundstate%tforce) then
                tlast = tlast .and. (dforcemax .lt. input%groundstate%epsforcescf)
                fm = forcemax
            end if

! check for STOP file
            if (rank==0) then
                Inquire (File='STOP', Exist=Exist)
                If (exist) Then
                    write(string,'("STOP file exists - stopping self-consistent loop")')
                    call printbox(60,"+",string)
                    tstop = .True.
                    tlast = .True.
                    Open (50, File='STOP')
                    Close (50, Status='DELETE')
                End If
            end if

        End If ! iscl>2

        call xmpi_bcast(mpiglobal, tstop)
        call xmpi_bcast(mpiglobal, tlast)

        call timesec(ts1)
        timeio = ts1 - ts0 + timeio
!! TIME - End of fourth IO segment
    End Do ! iscl
! end the self-consistent loop
20  Continue
    call stopwatch("exciting:scf", 0)

    Call timesec(ts0)
    If ((verbosity>-1).and.(rank==0)) Then
        write(string,'("Self-consistent loop stopped")')
        call printbox(60,"+",string)
    end if
    ! write density and potentials to file only if maxscl > 1
    If ((input%groundstate%maxscl > 1)) Then
        Call writestate
        If ((verbosity>-1).and.(rank==0)) Then
            Write (60, '(" STATE.OUT is written")')
        end if
    End If
! delete BROYDEN.OUT
    If (rank==0) then
        Inquire (File='BROYDEN.OUT', Exist=Exist)
        If (exist) Then
            Open (23, File='BROYDEN.OUT')
            Close (23, Status='DELETE')
        End If
    End If
    Call timesec(ts1)
    timeio = ts1 - ts0 + timeio

!------------------
! Compute forces
!------------------
    If (( .Not. tstop) .And. (input%groundstate%tforce)) Then
        Call force(input%groundstate%tfibs)
! For whatever reason each MPI process may produce very slightly different forces.
! At this spot, we equalise them, so that we do not end up with a different geometry
! for every process.
        call xmpi_allreduce( forcetot, mpiglobal )
        forcetot(1:3,1:natmtot)=forcetot(1:3,1:natmtot)/dble(procs)
! output forces to INFO.OUT
        if ((verbosity>-1).and.(rank==0)) then
           call printbox(60,"-","Writing atomic positions and forces")
           idm = 0
           write(60,*)
           write(60,'(" Atomic positions (",A,") :")') trim(acoord)
           do is = 1, nspecies
               do ia = 1, natoms (is)
                   idm = idm+1
                   if (input%structure%cartesian) then
                       write(60,'(" atom ",I5,2x,A2,T18,": ",3F14.8)') &
                      &  idm, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
                      &  atposc(:,ia,is)
                   else
                       write(60,'(" atom ",I5,2x,A2,T18,": ",3F14.8)') &
                      &  idm, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
                      &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
                   end if
               end do
           end do
           call writeforce(60,2)
        end if
    End If ! compute forces

    call stopwatch("exciting:post_scf", 1)
    ! set nwork to -2 to tell interface to call the deallocation functions
    If (rank .Eq. 0) Call mixerifc(input%groundstate%mixernumber, n, v, currentconvergence, -2)
    Deallocate(v)
    Call mpiresumeevecfiles()

    if (allocated(rhomtref)) deallocate(rhomtref)
    if (allocated(rhoirref)) deallocate(rhoirref)
    if( cdft_calculation%is_on() )  call deallocate_cdft_global_arrays()

    If ((verbosity>-1).and.(rank==0)) Then
! add blank line to TOTENERGY.OUT, FERMIDOS.OUT, MOMENT.OUT and RMSDVEFF.OUT
!      Write (62,*)
      If (spin_polarization) write (63,*)
! add blank line to DTOTENERGY.OUT, DFORCEMAX.OUT, CHGDIST.OUT and PCHARGE.OUT
!      Write (66,*)
!      If (input%groundstate%tforce) Write (67,*)
!      Write (68,*)
      if (input%groundstate%tpartcharges) write(69,*)
    End If

    If ((verbosity>-2).and.(rank==0)) Then
! write last total energy and add blank line to TOTENERGY.OUT
      Write (61, '(G22.12)') engytot
      Write (61,*)
      Call flushifc(61)
! add blank line to RMSDVEFF.OUT and DFSCFMAX.OUT
      Write (65,*)
      Call flushifc(65)
      if (input%groundstate%tforce) then
          Write (67,*)
          Call flushifc(67)
      end if
    End If

    if (associated(input%groundstate%xsLO).and.(rank==0)) then
      call genxsLOs()
    end if

    call mt_hscf%release()
    call stopwatch("exciting:post_scf", 0)
end subroutine
