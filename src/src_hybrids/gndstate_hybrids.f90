!
Subroutine gndstate_hybrids

    Use modinput
    Use modmain
    Use modmpi
    Use scl_xml_out_Module
    use modgw
    use mod_hartreefock

    Implicit None
! time measurements
    Real(8) :: timetot, ts0, ts1, tsg0, tsg1, tin1, tin0
    real(8) :: et
    integer :: ik

!-------------------------------------------------------------------------------

    fgw = 600
    if (rank .Eq. 0) then
! open INFO.OUT file
        open(60, File='INFO'//trim(filext), Action='WRITE', Form='FORMATTED')
! Hartree-Fock related debugging info
        open(fgw, File='HARTFOCK.DBG', Action='WRITE', Form='FORMATTED')
    end if
! only tetrahedron method is implemented (LIBBZINT)
    input%groundstate%stypenumber = -1
    
!! TIME - Initialisation segment
    Call timesec (tsg0)
    Call timesec (ts0)

! initialise global variables
    Call timesec (tin0)
    Call init0
    Call timesec (tin1)
    time_init0=tin1-tin0
    Call timesec (tin0)
    Call init1
    Call timesec (tin1)
    time_init1=tin1-tin0
    
!-------------------
! print info
!-------------------
    if (rank==0) then
! write out general information to INFO.OUT
        call writeinfo(60)
! write the real and reciprocal lattice vectors to file
        Call writelat
! write interatomic distances to file
        Call writeiad (.False.)
! write symmetry matrices to file
        Call writesym
#ifdef XS
! write relation to inverse symmetries
        Call writesymi
! write advanced information on symmetry group
        Call writesym2
! write out symmetrization matrix for rank 2 tensors
        Call writesymt2
#endif
! output the k-point set to file
        Call writekpts
! write lattice vectors and atomic positions to file
        Call writegeometryxml(.False.)
! open TOTENERGY.OUT
        Open (61, File='TOTENERGY'//trim(filext), Action='WRITE', Form='FORMATTED')
! open FERMIDOS.OUT
        Open (62, File='FERMIDOS'//trim(filext), Action='WRITE', Form='FORMATTED')
! open MOMENT.OUT if required
        If (associated(input%groundstate%spin)) open (63, file='MOMENT'//trim(filext), action='WRITE', form='FORMATTED')
! open RMSDVEFF.OUT
        Open (65, File='RMSDVEFF'//trim(filext), Action='WRITE', Form='FORMATTED')
! open DTOTENERGY.OUT
        open(66,file='DTOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! open CHGDIST.OUT
        open(68,file='CHGDIST'//trim(filext),action='WRITE',form='FORMATTED')
! open PCHARGE.OUT
        if (input%groundstate%tpartcharges) &
       &  open(69,file='PCHARGE'//trim(filext),action='WRITE',form='FORMATTED')
    end if ! rank

    Call timesec (ts1)
    timeinit = timeinit+ts1-ts0
!! TIME - End of initialisation segment

    if (rank==0) then
        write(60,*)
        write(60, '("+-----------------------------------------------------------+")')
        write(60, '("| Hybrids module started ")')
        write(60, '("+-----------------------------------------------------------+")')
    end if

!---------------------------------------
! initialise HF related input parameters
!---------------------------------------
    call init_hf_input_params

!---------------------------------------
!   Initialize k/q grids
!---------------------------------------
    Call init_kqpts

! non-local exchange energy
    allocate(exnlk(nkpt))
    exnlk(:) = 0.d0

! BZ integrated non-local energy
    exnl = 0.d0

! non-local exchange potential
    allocate(vxnl(nstfv,nstfv,nkpt))
    vxnl(:,:,:) = zzero

! eigenvector from previous interation
    allocate(evecfv0(nmatmax,nstfv,nspnfv,nkpt))
    evecfv0(:,:,:,:) = zzero

!----------------------------------------!
! begin the self-consistent loop
!----------------------------------------!
    et = 0.d0

    ihyb = 0
    do ihyb = 1, input%groundstate%maxscl

        If (rank==0) Then
            Write (60,*)
            Write (60, '("+-------------------------+")')
            Write (60, '("| (Hybrids) Iteration number : ", I4, " |")') ihyb
            Write (60, '("+-------------------------+")')
            Call flushifc(60)
        End If
        If (ihyb >= input%groundstate%maxscl) Then
            If (rank==0) Then
                Write (60,*)
                Write (60, '("Reached hybrids scf loops maximum")')
                call warning('Warning(gndstate_hybrids): Reached self-consistent loops maximum')
                Call flushifc(60)
            End If
            ! exit i-loop
            exit
        End If

!------------------------------------!
!   SCF cycle
!------------------------------------!
        if (ihyb>1) task = 1
        call scf_cycle
        
!--------------------------------------------!
!   Matrix elements of non-local potential   !
!--------------------------------------------!
        if (rank .Eq. 0) call boxmsg(60,'-','Non-local exchange potential')
        call cpu_time(ts0)
        call calc_vxnl
        call cpu_time(ts1)
        if (rank == 0) then
            call write_cputime(60,ts1-ts0,'calc_vxnl')
            call boxmsg(60,'-','')
        end if

! update the non-local energy
        engytot = engytot-ex_coef*exnl
        exnl = 0.d0
        do ik = 1, nkpt
            exnl = exnl+wkpt(ik)*exnlk(ik)
        end do
        engytot = engytot+ex_coef*exnl

! update convergence criteria
        deltae=abs(et-engytot)

!------------------------------------------------
! Store the eigenvectors from previous iteration
!------------------------------------------------
        do ik = 1, nkpt
            call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv0(:,:,:,ik))
        end do

! Check for convergence
        if (ihyb > 1) then

            write(60,'("(Hybrids) Absolute change in total energy (target)   : ",G18.10," (",G18.10,")")') &
           &  deltae, input%groundstate%epsengy
            write(60,'("(Hybrids) Charge distance (target)                   : ",G18.10," (",G18.10,")")') &
           &  chgdst, input%groundstate%epschg

            if (deltae < input%groundstate%epsengy) Then
                write(60,*)
                write(60,'("(Hybrids) Convergence targets achieved")')
                ! exit ihyb-loop
                exit
            End If
            et = engytot
            
! output the current total time
            timetot = timeinit + timemat + timefv + timesv + timerho &
           &        + timepot + timefor + timeio + timemt + timemixer
            
            if (rank==0) then
                write(60,*)
                write(60, '("Wall time (seconds) : ", F12.2)') timetot
            end if

       end if ! ihyb>1

    end do ! i

    if (rank==0) then
        write(60,*)
        write(60, '("+-----------------------------------------------------------+")')
        write(60, '("| Hybrids module stopped ")')
        write(60, '("+-----------------------------------------------------------+")')
        write(60,*)
    end if




! generate the new species files with the optimized linearization energies
    If ((rank .Eq. 0).and.(input%groundstate%tspecies)) Call updatespecies

!-------------------------------------!
!   output timing information
!-------------------------------------!
!! TIME - Fifth IO segment
    call timesec(ts0)
    if (rank==0) then
! close the TOTENERGY.OUT file
        close(61)
! close the FERMIDOS.OUT file
        close(62)
! close the MOMENT.OUT file
        if (associated(input%groundstate%spin)) close(63)
! close the RMSDVEFF.OUT file
        close(65)
! close the DTOTENERGY.OUT file
        close(66)
! close the CHGDIST.OUT file
        close(68)
! close the PCHARGE.OUT file
        if (input%groundstate%tpartcharges) close(69)
    end if
    
! xml output
    if (rank==0) then
        call structure_xmlout
        call scl_xml_setGndstateStatus("finished")
        call scl_xml_out_write
    end if
    call timesec(ts1)
    timeio=ts1-ts0+timeio
!! TIME - End of fifth IO segment

    If (rank .Eq. 0) then
        Write (60,*)
        Write (60, '("Timings (seconds) :")')
        Write (60, '(" initialisation", T40, ": ", F12.2)') timeinit
        Write (60, '("            - init0", T40,": ", F12.2)') time_init0
        Write (60, '("            - init1", T40,": ", F12.2)') time_init1
        Write (60, '("            - rhoinit", T40,": ", F12.2)') time_density_init
        Write (60, '("            - potential initialisation", T40,": ", F12.2)') time_pot_init
        Write (60, '("            - others", T40,": ", F12.2)') timeinit-time_init0-time_init1-time_density_init-time_pot_init
        Write (60, '(" Hamiltonian and overlap matrix set up", T40, ": ", F12.2)') timemat
        Write (60, '("            - hmlaan", T40,": ", F12.2)') time_hmlaan
        Write (60, '("            - hmlalon", T40,": ", F12.2)') time_hmlalon
        Write (60, '("            - hmllolon", T40,": ", F12.2)') time_hmllolon
        Write (60, '("            - olpaan", T40,": ", F12.2)') time_olpaan
        Write (60, '("            - olpalon", T40,": ", F12.2)') time_olpalon
        Write (60, '("            - olplolon", T40,": ", F12.2)') time_olplolon
        Write (60, '("            - hmlistln", T40,": ", F12.2)') time_hmlistln
        Write (60, '("            - olpistln", T40,": ", F12.2)') time_olpistln
        Write (60, '(" first-variational secular equation", T40, ": ", F12.2)') timefv
        If (associated(input%groundstate%spin)) Then
            Write (60, '(" second-variational calculation", T40, ": ", F12.2)') timesv
        End If
        Write (60, '(" charge density calculation", T40, ": ", F12.2)') timerho
        Write (60, '(" potential calculation", T40, ": ", F12.2)') timepot
        Write (60, '(" muffin-tin manipulations", T40, ": ", F12.2)') timemt
        Write (60, '(" APW matching", T40, ": ", F12.2)') timematch
        Write (60, '(" disk reads/writes", T40, ": ", F12.2)') timeio
        Write (60, '(" mixing efforts", T40, ": ", F12.2)') timemixer
        If (input%groundstate%tforce) Write (60, '(" force calculation", T40, ": ", F12.2)') timefor
        timetot = timeinit + timemat + timefv + timesv + timerho + &
        &         timepot + timefor+timeio+timemt+timemixer+timematch
        Write (60, '(" sum", T40, ": ", F12.2)') timetot
        Call timesec(tsg1)
        Write (60, '(" total", T40, ": ", F12.2)') tsg1-tsg0
        Write (60,*)
        Write (60, '("More timings (seconds)")')
        Write (60, '(" Dirac eqn solver", T40, ": ", F12.2)') time_rdirac
        Write (60, '(" Rel. Schroedinger eqn solver", T40, ": ", F12.2)') time_rschrod
        Write (60, '(" Total time spent in radial solvers", T40, ": ", F12.2)') time_rdirac+time_rschrod
        If (input%groundstate%xctypenumber .Lt. 0) Then 
            Write (60, '(" Time spent for oepvnl ", T40,": ", F12.2)') time_oepvnl
            Write (60, '(" Time spent for oep iteration ", T40,": ", F12.2)') time_oep_iter
        End If
        Write (60,*)
        Write (60, '("+-----------------------------------------------------------+")')
        Write (60, '("| EXCITING ", a, " stopped")') trim(versionname)
        Write (60, '("+-----------------------------------------------------------+")')
        close (60)
    endif
   
    Return
   End Subroutine
