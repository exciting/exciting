!
!BOP
! !ROUTINE: hybrids
! !INTERFACE:
!
Subroutine hybrids
! !USES:
    Use modinput
    Use modmain
    Use modmpi
    Use scl_xml_out_Module
    Use mod_hybrids
!
! !DESCRIPTION:
!   Main routine for Hartree-Fock based hybrid functionals.
!
! !REVISION HISTORY:
!   Created September 2013 (DIN)
!   Modified October 2013
!   Modified Februar 2014
!EOP
!BOC
    Implicit None

    integer :: ik, Recl
    real(8) :: et
! time measurements
    Real(8) :: timetot, ts0, ts1, tsg0, tsg1, tin1, tin0, time_hyb
    character*(77) :: string
    
    ! Charge distance
    Real (8), Allocatable :: rhomtref(:,:,:) ! muffin-tin charge density (reference)
    Real (8), Allocatable :: rhoirref(:)     ! interstitial real-space charge density (reference)
    Logical :: restart
    integer::xctype_(3)    
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
         
! require forces for structural optimisation
    If ((task .Eq. 2) .Or. (task .Eq. 3)) input%groundstate%tforce = .True.    
! chech if restart should be performed
    restart = input%groundstate%Hybrid%restart
    
!-------------------
! print info
!-------------------
    fgw = 600
    if (rank==0) then
! open INFO.OUT file
        open(60, File='INFO'//trim(filext), Action='WRITE', Form='FORMATTED')
! Hartree-Fock related debugging info
        open(fgw, File='HYBRIDS.OUT', Action='WRITE', Form='FORMATTED')
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

!___________________
! Open support files

! open TOTENERGY.OUT
        Open (61, File='TOTENERGY'//trim(filext), Action='WRITE', Form='FORMATTED')

! open RMSDVEFF.OUT
        Open (65, File='RMSDVEFF'//trim(filext), Action='WRITE', Form='FORMATTED')

! open MOMENT.OUT if required
        If (associated(input%groundstate%spin)) then
            open (63, file='MOMENT'//trim(filext), action='WRITE', form='FORMATTED')
        end if

    end if ! rank

    Call timesec (ts1)
    timeinit = timeinit+ts1-ts0
!! TIME - End of initialisation segment

!------------------------------------!
!   Hybrids cycle
!------------------------------------!
    time_hyb = 0.d0
    call timesec(ts0)

    if (rank==0) then
        write(string,'("Hybrids module started")') 
        call printbox(60,"*",string)
        call flushifc(60)
    end if

!---------------------------------------
! initialise HF related input parameters
!---------------------------------------
    call init_hybrids

!---------------------------------------
!   Initialize k/q grids
!---------------------------------------
    Call init_kqpts_hybrids

!--------------------------------------------------------------
! Calculate the integrals to treat the singularities at G+q->0
!--------------------------------------------------------------
    call setsingc
    
    call timesec(ts1)
    time_hyb = time_hyb+ts1-ts0
    
!------------------------------------------!
! begin the (external) self-consistent loop
!------------------------------------------!

!_____________________________
! reference density

    If (allocated(rhomtref)) deallocate(rhomtref)
    Allocate(rhomtref(lmmaxvr,nrmtmax,natmtot))
    If (allocated(rhoirref)) deallocate(rhoirref)
    Allocate (rhoirref(ngrtot))

    do ihyb = 0, input%groundstate%Hybrid%maxscl
    
        ! exit self-consistent loop if last iteration is complete
        If (ihyb >= input%groundstate%Hybrid%maxscl) Then
            If (rank==0) Then
                write(string,'("Reached hybrids self-consistent loops maximum : ", I4)') &
               &  input%groundstate%Hybrid%maxscl
                call printbox(60,"+",string)
                call warning('Warning(hybrids): Reached self-consistent loops maximum')
                Call flushifc(60)
            End If
            exit ! exit ihyb-loop
        End If
        If (rank==0) Then
            write(string,'("Hybrids iteration number : ", I4)') ihyb
            call printbox(60,"+",string)
            Call flushifc(60)
        End If

        ! hybrids always start after normal DFT self-consistent run
        if (ihyb==0) then
          ex_coef = 0.d0
          ec_coef = 1.d0
          ! for Libxc use PBE from Libxc
          if (xctype(1)==100) then
              xctype_=xctype
              xctype=(/100,101,130/)
          end if
        else
          task = 7
          ! restore settings for hybrid functional
          ex_coef = input%groundstate%Hybrid%excoeff
          ec_coef = input%groundstate%Hybrid%eccoeff
          if (xctype(1)==100) then
              xctype=xctype_
          end if
          ! settings for convergence and mixing
          input%groundstate%mixerswitch = 1
          input%groundstate%scfconv = 'charge'
          rhomtref(:,:,:) = rhomt(:,:,:)
          rhoirref(:) = rhoir(:)
        end if
        
!---------------------------
! KS self-consistent run
!---------------------------
        call scf_cycle(-1)
       
        ! some output        
        if (rank==0) then
            call writeengy(60)
            write(60,*)
            write(60,'(" DOS at Fermi energy (states/Ha/cell)",T45 ": ", F18.8)') fermidos
            call writechg(60,input%groundstate%outputlevelnumber)
            if (fermidos<1.0d-4) call printbandgap(60)
            call flushifc(60)
        end if

! check for convergence
        if (ihyb>0) then
            deltae = dabs(et-engytot)
            call chgdist(rhomtref,rhoirref)
            if (rank==0) Then
              write(60,*)
              write(60,'(" Absolute change in total energy   (target) : ",G18.10," (",G18.10,")")') &
              &     deltae, input%groundstate%epsengy
              write(60,'(" Charge distance                   (target) : ",G18.10," (",G18.10,")")') &
              &     chgdst, input%groundstate%epschg
            end if
            if (chgdst .lt. input%groundstate%epschg) then
                if (rank==0) Then
                    write(string,'("Convergence target is reached")')
                    call printbox(60,"+",string)
                    Call flushifc(60)
                end If
                exit ! exit ihyb-loop
            end if
        end if
        et = engytot

!-----------------------------------
! calculate the non-local potential
!-----------------------------------
        if ((ihyb==0).and.restart) then
          call timesec(ts0)
          Call getvnlmat 
          call timesec(ts1)
          if (rank==0) then
              call write_cputime(60,ts1-ts0, 'READ_VNLMAT')
              write(60,*)
          end if    
          time_hyb = time_hyb+ts1-ts0
          
        else if  (ihyb < input%groundstate%Hybrid%maxscl-1) Then
          !------------------------------------------
          call timesec(ts0)
          call calc_vxnl
          if (rank==0) write(*,*) 'vxnl=', sum(vxnl)
          if (input%groundstate%tevecsv.and.(rank==0)) write(*,*) 'bxnl=', sum(bxnl)
          call timesec(ts1)
          if (rank==0) then
              write(60,*)
              call write_cputime(60,ts1-ts0, 'CALC_VXNL')
          end if
          !------------------------------------------
          call timesec(ts0)
          call calc_vnlmat
          if (rank==0) write(*,*) 'calc_vnlmat=', sum(vnlmat)
          call timesec(ts1)
          if (rank==0) then
              call write_cputime(60,ts1-ts0, 'CALC_VNLMAT')
              write(60,*)
          end if
!-----------------------------------
! update radial functions
!-----------------------------------
          if ((input%groundstate%Hybrid%updateRadial).and.(ihyb>0)) call updateradial

          time_hyb = time_hyb+ts1-ts0
          
        end if
        
! output the current total time
        timetot = timeinit + timemat + timefv + timesv + timerho  &
        &       + timepot + timefor + timeio + timemt + timemixer &
        &       + time_hyb
            
        if (rank==0) then
            write(60, '(" Wall time (seconds) : ", F12.2)') timetot
        end if

    end do ! ihyb

    if (rank==0) then
        write(string,'("Hybrids module stopped")') 
        call printbox(60,"+",string)
        call flushifc(60)
    end if

!-------------------------------------!
!   output timing information
!-------------------------------------!
!! TIME - Fifth IO segment
    call timesec(ts0)
    if (rank==0) then

!____________________
! Close support files

! close the TOTENERGY.OUT file
        close(61)
! close the RMSDVEFF.OUT file
        close(65)
! close the MOMENT.OUT file
        if (associated(input%groundstate%spin)) close(63)

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
    Call timesec(tsg1)
    If ((rank .Eq. 0).and.(input%groundstate%outputlevelnumber>1)) then
        write(string,'("Timings (seconds)")') 
        call printbox(60,"-",string)
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
        Write (60, '(" Dirac eqn solver", T40, ": ", F12.2)') time_rdirac
        Write (60, '(" Rel. Schroedinger eqn solver", T40, ": ", F12.2)') time_rschrod
        Write (60, '(" Total time spent in radial solvers", T40, ": ", F12.2)') time_rdirac+time_rschrod
        write (60,*)
    end if

    If (rank .Eq. 0) then
        Write (60, '(" Total time spent (seconds)", T40, ": ", F12.2)') tsg1-tsg0
        if (lwarning) call printbox(60,"-","CAUTION! Warnings have been written in file WARNING.OUT !")
        write(string,'("EXCITING ", a, " stopped")') trim(versionname)
        call printbox(60,"=",string)
        close (60)
    endif
   
!----------------------------------------
! Save HF energies into binary file
!----------------------------------------

    Inquire (IoLength=Recl) nkptnr, nstsv, vklnr(:,1), evalsv(:,1)
    Open (70, File='EVALHF.OUT', Action='WRITE', Form='UNFORMATTED', &
   &   Access='DIRECT', status='REPLACE', Recl=Recl)
    do ik = 1, nkptnr
        write(70, Rec=ik) nkptnr, nstsv, vklnr(:,ik), evalsv(:,indkp(ik))-efermi
    end do ! ik
    Close(70)

!----------------------------------------
! Clean data
!----------------------------------------
    call exit_hybrids
      
    Return
    
End Subroutine
!EOC

