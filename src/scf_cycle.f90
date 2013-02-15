!BOP
! !ROUTINE: scf_cycle
! !INTERFACE:
!
!
subroutine scf_cycle
! !USES:
    Use modinput
    Use modmain
    Use modmpi
    Use scl_xml_out_Module
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!   Created February 2013 (DIN)
!EOP
!BOC
    Implicit None

    Real (8) :: et, fm
    Real (8), Allocatable :: evalfv(:, :)
    Complex (8), Allocatable :: evecfv(:, :, :)
    Complex (8), Allocatable :: evecsv(:, :)
    Logical :: force_converged, tibs, exist
    Integer :: ik, is, ia, idm
    Integer :: n, nwork
    Real(8), Allocatable :: v(:)
    Real(8) :: timetot, ts0, ts1, tin1, tin0

    If (rank .Eq. 0) Then
        Write (60,*)
        Write (60, '("+------------------------------+")')
        Write (60, '("| Self-consistent loop started |")')
        Write (60, '("+------------------------------+")')
        Write (60,*)
    End If

! initialize reference density
    If (allocated(rhomtref)) deallocate(rhomtref)
    Allocate(rhomtref(lmmaxvr,nrmtmax,natmtot))
    If (allocated(rhoirref)) deallocate(rhoirref)
    Allocate (rhoirref(ngrtot))

!! TIME - Begin of initialisation segment 
    Call timesec (ts0)
! initialise or read the charge density and potentials from file
    If ((task .Eq. 1) .Or. (task .Eq. 3)) Then
        Call readstate
        If (rank .Eq. 0) write (60, '("Potential read in from STATE.OUT")')
    Else If (task .Eq. 200) Then
        Call phveff
        If (rank .Eq. 0) write (60, '("Supercell potential constructed from STATE.OUT")')
    Else
        Call timesec(tin0)
        Call rhoinit
        Call timesec(tin1)
        time_density_init=tin1-tin0
! store density to reference
        rhoirref(:)=rhoir(:)
        rhomtref(:,:,:)=rhomt(:,:,:)
        Call timesec(tin0)
        Call poteff
        Call genveffig
        Call timesec(tin1)
        time_pot_init=tin1-tin0
        If (rank .Eq. 0) write (60, '("Density and potential initialised from atomic data")')
    End If
    If (rank .Eq. 0) then
        write (60, *)
        Call flushifc (60)
    end if
    Call timesec (ts1)
    timeinit = timeinit+ts1-ts0
!! TIME - End of initialisation segment    

!! TIME - Mixer segment
    Call timesec (ts0)
! size of mixing vector
    n = lmmaxvr*nrmtmax*natmtot+ngrtot
    If (associated(input%groundstate%spin)) n = n*(1+ndmag)
    If (ldapu .Ne. 0) n = n + 2*lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot
! allocate mixing arrays
    Allocate (v(n))
! call mixing array allocation functions by setting
    nwork = -1
! and call interface
    If (rank .Eq. 0) Call mixerifc(input%groundstate%mixernumber, n, v, currentconvergence, nwork)
    Call timesec (ts1)
    timemixer = ts1-ts0+timemixer
!! TIME - End of mixer segment

! set last iteration flag
    tlast = .False.
! set stop flag
    tstop = .False.
    et = 0.d0
    fm = 0.d0

! delete any existing eigenvector files
    If ((rank .Eq. 0) .And. ((task .Eq. 0) .Or. (task .Eq. 2))) Call delevec
! begin the self-consistent loop
    If (rank .Eq. 0) Then
        Write (60,*)
        Write (60, '("+------------------------------+")')
        Write (60, '("| Self-consistent loop started |")')
        Write (60, '("+------------------------------+")')
    End If

!! TIME - First IO segment
    Call timesec (ts0)
!----------------------------------------!
! begin the self-consistent loop
!----------------------------------------!
    iscl = 0
    Do iscl = 1, input%groundstate%maxscl
!
        If (rank .Eq. 0) Then
            Write (60,*)
            Write (60, '("+-------------------------+")')
            Write (60, '("| Iteration number : ", I4, " |")') iscl
            Write (60, '("+-------------------------+")')
        End If
        If (iscl .Ge. input%groundstate%maxscl) Then
            If (rank .Eq. 0) Then
                Write (60,*)
                Write (60, '("Reached self-consistent loops maximum")')
                call warning('Warning(gndstate): Reached self-consistent loops maximum')
            End If
            tlast = .True.
        End If
        If (rank .Eq. 0) Call flushifc (60)
        Call timesec (ts1)
        timeio=ts1-ts0+timeio
!! TIME - End of first IO segment

!! TIME - Muffin-tin segment
        Call timesec (ts0)
! generate the core wavefunctions and densities
        Call gencore
! find the new linearization energies
        Call linengy
        if (rank .eq. 0) Call writelinen
! generate the APW radial functions
        Call genapwfr
! generate the local-orbital radial functions
        Call genlofr(tlast)
! compute the overlap radial integrals
        Call olprad
! compute the Hamiltonian radial integrals
        Call hmlrad
! zero partial charges
        chgpart(:,:,:)=0.d0
        Call timesec (ts1)
        timemt=ts1-ts0+timemt
!! TIME - End of muffin-tin segment

!! TIME - Second IO segment       
        Call timesec (ts0)
#ifdef MPI
        Call MPI_barrier (MPI_COMM_WORLD, ierr)
        If (rank .Eq. 0) Call delevec ()
#endif
#ifdef MPISEC
        splittfile = .True.
        Do ik = firstk (rank), lastk (rank)
#endif
#ifdef NEVERDEFINED
        End Do
#endif
#ifndef MPISEC
        splittfile = .False.
! begin parallel loop over k-points
#ifdef KSMP
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv)
!$OMP DO
#endif
        Do ik = 1, nkpt
#endif
! every thread should allocate its own arrays
            Allocate (evalfv(nstfv, nspnfv))
            Allocate (evecfv(nmatmax, nstfv, nspnfv))
            Allocate (evecsv(nstsv, nstsv))
!! TIME - seceqn does not belong to IO
            Call timesec(ts1)
            timeio=ts1-ts0+timeio            
! solve the first- and second-variational secular equations
            Call seceqn (ik, evalfv, evecfv, evecsv)
            Call timesec(ts0)
! write the eigenvalues/vectors to file
            Call putevalfv (ik, evalfv)
            Call putevalsv (ik, evalsv(:, ik))
            Call putevecfv (ik, evecfv)
            Call putevecsv (ik, evecsv)
! calculate partial charges
            !call genpchgs(ik,evecfv,evecsv)
            Deallocate (evalfv, evecfv, evecsv)
        End Do
#ifdef KSMP
!$OMP END DO
!$OMP END PARALLEL
#endif
#ifdef MPISEC
        Call mpi_allgatherv_ifc(nkpt,nstsv,rbuf=evalsv)
        Call MPI_barrier(MPI_COMM_WORLD, ierr)
#endif
! find the occupation numbers and Fermi energy
        Call occupy
        If (rank .Eq. 0) Then
! write out the eigenvalues and occupation numbers
            Call writeeval
! write the Fermi energy to file
            Call writefermi
        End If
! set the charge density and magnetisation to zero
        rhomt (:, :, :) = 0.d0
        rhoir (:) = 0.d0
        If (associated(input%groundstate%spin)) Then
            magmt (:, :, :, :) = 0.d0
            magir (:, :) = 0.d0
        End If
#ifdef MPIRHO
        Do ik = firstk (rank), lastk (rank)
!write the occupancies to file
            Call putoccsv (ik, occsv(:, ik))
        End Do
        Do ik = firstk (rank), lastk (rank)
#endif
#ifndef MPIRHO
        If (rank .Eq. 0) Then
            Do ik = 1, nkpt
!write the occupancies to file
                Call putoccsv (ik, occsv(:, ik))
            End Do
        End If
#ifdef KSMP
! begin parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecfv,evecsv)
!$OMP DO
#endif
        Do ik = 1, nkpt
#endif
            Allocate (evecfv(nmatmax, nstfv, nspnfv))
            Allocate (evecsv(nstsv, nstsv))
! get the eigenvectors from file
            Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
            Call getevecsv (vkl(:, ik), evecsv)
!! TIME - rhovalk does not belong to IO
            Call timesec(ts1)
            timeio=ts1-ts0+timeio
! add to the density and magnetisation
            Call rhovalk (ik, evecfv, evecsv)
            Deallocate (evecfv, evecsv)
            Call timesec(ts0)
        End Do
#ifndef MPIRHO
#ifdef KSMP
!$OMP END DO
!$OMP END PARALLEL
#endif
#endif
#ifdef MPIRHO
        Call mpisumrhoandmag
#endif
#ifdef MPI
        If ((input%groundstate%xctypenumber.Lt.0).Or.(xctype(2).Ge.400).Or.(xctype(1).Ge.400)) &
       &    Call mpiresumeevecfiles()
#endif
        If (rank .Eq. 0) Then
! write out partial charges
            !call writepchgs(69,input%groundstate%lmaxvr)
            call flushifc(69)
        end if

        Call timesec(ts1)
        timeio=ts1-ts0+timeio
!! TIME - End of second IO segment

!! TIME - potential
        Call timesec(ts0)
! symmetrise the density
        Call symrf(input%groundstate%lradstep, rhomt, rhoir)
! symmetrise the magnetisation
        If (associated(input%groundstate%spin)) Call symrvf(input%groundstate%lradstep, magmt, magir)
! convert the density from a coarse to a fine radial mesh
        Call rfmtctof (rhomt)
! convert the magnetisation from a coarse to a fine radial mesh
        Do idm = 1, ndmag
            Call rfmtctof (magmt(:, :, :, idm))
        End Do
! add the core density to the total density
        Call addrhocr
! calculate the charges
        Call charge
! calculate the moments
        If (associated(input%groundstate%spin)) Call moment
! normalise the density
        Call rhonorm
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
        call chgdist
! store density to reference
        rhoirref(:)=rhoir(:)
        rhomtref(:,:,:)=rhomt(:,:,:)
! compute the effective potential
        Call poteff
! pack interstitial and muffin-tin effective potential and field into one array
        Call packeff (.True., n, v)
! mix in the old potential and field with the new
        If (rank .Eq. 0) Call mixerifc (input%groundstate%mixernumber, n, v, currentconvergence, nwork)
#ifdef MPI
        Call MPI_bcast (v(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
! unpack potential and field
        Call packeff (.False., n, v)
! add the fixed spin moment effect field
        If (getfixspinnumber() .Ne. 0) Call fsmfield
! Fourier transform effective potential to G-space
        Call genveffig
! reduce the external magnetic fields if required
        If (associated(input%groundstate%spin)) Then
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
! compute the energy components
        Call energy
        Call timesec(ts1)
        timepot=ts1-ts0+timepot
!! TIME - End of potential

!! TIME - Third IO segment
! compute the forces (without IBS corrections)
        Call timesec(ts0)
        If (input%groundstate%tforce) Then
            tibs=input%groundstate%tfibs
            input%groundstate%tfibs=.false.
!! forces do no belong to IO
            Call timesec(ts1)
            timeio=ts1-ts0+timeio
            Call force
            Call timesec(ts0)
            input%groundstate%tfibs=tibs
            If (rank .Eq. 0) Then
! output forces to INFO.OUT
                Call writeforce (60)
! write maximum force magnitude to FORCEMAX.OUT
                Write (64, '(G18.10,"   (without IBS correction)")') forcemax
                Call flushifc (64)
            end if
        end if
        If (rank .Eq. 0) Then
! output energy components
            Call writeengy (60)
            Write (60,*)
            Write (60, '("Density of states at Fermi energy : ", G18.10)') fermidos
            Write (60, '(" (states/Hartree/unit cell)")')
! write total energy to TOTENERGY.OUT and flush
            Write (61, '(G22.12)') engytot
            Call flushifc (61)
! write DOS at Fermi energy to FERMIDOS.OUT and flush
            Write (62, '(G18.10)') fermidos
            Call flushifc (62)
! output charges and moments
            Call writechg (60)
! write total moment to MOMENT.OUT and flush
            If (associated(input%groundstate%spin)) Then
                Write (63, '(3G18.10)') momtot (1:ndmag)
                Call flushifc (63)
            End If
! output effective fields for fixed spin moment calculations
            If (getfixspinnumber() .Ne. 0) Call writefsm (60)
! check for WRITE file
            Inquire (File='WRITE', Exist=exist)
            If (exist) Then
                Write (60,*)
                Write (60, '("WRITE file exists - writing STATE.OUT")')
                Call writestate
                Open (50, File='WRITE')
                Close (50, Status='DELETE')
            End If
! write STATE.OUT file if required
            If (input%groundstate%nwrite .Ge. 1) Then
                If (Mod(iscl, input%groundstate%nwrite) .Eq. 0) Then
                    Call writestate
                    Write (60,*)
                    Write (60, '("Wrote STATE.OUT")')
                End If
            End If
! update convergence criteria
            deltae=abs(et-engytot)
            dforcemax=abs(fm-forcemax)
            Call scl_iter_xmlout ()
            If (associated(input%groundstate%spin)) Call scl_xml_write_moments()
            Call scl_xml_out_write()
        End If
        Call timesec(ts1)
        timeio=ts1-ts0+timeio
!! TIME - End of third IO segment

! exit self-consistent loop if last iteration is complete
        If (tlast) Go To 20

!! TIME - Fourth IO segment
        Call timesec(ts0)
        If (rank .Eq. 0) Then
! check for convergence
            If (iscl .Ge. 2) Then
                Write (60,*)
                Write (60, '("RMS change in effective potential (target) : ", G18.10, " (", G18.10, ")")') &
               &    currentconvergence, input%groundstate%epspot
                write(60,'("Absolute change in total energy (target)   : ",G18.10," (",G18.10,")")') &
               &    deltae, input%groundstate%epsengy
                if (input%groundstate%tforce) then
                    write(60,'("Absolute change in |max. force| (target)   : ",G18.10," (",G18.10,")")') &
                   &    dforcemax, input%groundstate%epsforce
                end if
                write(60,'("Charge distance (target)                   : ",G18.10," (",G18.10,")")') &
                &   chgdst, input%groundstate%epschg
                write(66,'(G18.10)') deltae
                call flushifc(66)
                if (input%groundstate%tforce) then
                    write(67,'(G18.10)') dforcemax
                    call flushifc(67)
                end if
                write(68,'(G18.10)') chgdst
                call flushifc(68)
                Write (65, '(G18.10)') currentconvergence
                Call flushifc (65)
                If ((currentconvergence .Lt. input%groundstate%epspot).and. &
               &    (deltae .lt. input%groundstate%epsengy).and. &
               &    (chgdst .lt. input%groundstate%epschg).and. &
               &    (dforcemax .lt. input%groundstate%epsforce)) Then
                    Write (60,*)
                    Write (60, '("Convergence targets achieved")')
                    tlast = .True.
                End If
            End If
            et = engytot
            fm = forcemax
            If ((input%groundstate%xctypenumber .Lt. 0).Or.(xctype(2) .Ge. 400).Or.(xctype(1) .Ge. 400)) Then
                Write (60,*)
                Write (60, '("Magnitude of OEP residual : ", G18.10)') resoep
            End If
! check for STOP file
            Inquire (File='STOP', Exist=Exist)
            If (exist) Then
                Write (60,*)
                Write (60, '("STOP file exists - stopping self-consistent loop")')
                tstop = .True.
                tlast = .True.
                Open (50, File='STOP')
                Close (50, Status='DELETE')
            End If
! output the current total time
            timetot = timeinit + timemat + timefv + timesv + timerho &
           &        + timepot + timefor+timeio+timemt+timemixer
            Write (60,*)
            Write (60, '("Wall time (seconds) : ", F12.2)') timetot
        End If
#ifdef MPI
        Call MPI_bcast (tstop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        Call MPI_bcast (tlast, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
        Call timesec(ts1)
        timeio=ts1-ts0+timeio
!! TIME - End of fourth IO segment
    End Do
! end the self-consistent loop

20  Continue

    Call timesec(ts0)
    If (rank .Eq. 0) Then
        Write (60,*)
        Write (60, '("+------------------------------+")')
        Write (60, '("| Self-consistent loop stopped |")')
        Write (60, '("+------------------------------+")')
! write density and potentials to file only if maxscl > 1
        If (input%groundstate%maxscl .Gt. 1) Then
            Call writestate
            Write (60,*)
            Write (60, '("Wrote STATE.OUT")')
        End If
    end if
    Call timesec(ts1)
    timeio=ts1-ts0+timeio   
    
! compute forces
    If (( .Not. tstop) .And. (input%groundstate%tforce)) Then
        Call force
        If (rank .Eq. 0) Then
! output forces to INFO.OUT
            Call writeforce (60)
! write maximum force magnitude to FORCEMAX.OUT
            Write (64, '(G18.10)') forcemax
            Call flushifc (64)
        End If
    End If

! set nwork to -2 to tell interface to call the deallocation functions
    If (rank .Eq. 0) Call mixerifc (input%groundstate%mixernumber, n, v, currentconvergence, -2)
    Deallocate (v)
    Call mpiresumeevecfiles ()

    if (allocated(rhomtref)) deallocate(rhomtref)
    if (allocated(rhoirref)) deallocate(rhoirref)

end subroutine
