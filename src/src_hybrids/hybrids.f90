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
    integer :: ihyb, maxscl
    real(8) :: et
! time measurements
    real(8) :: timetot, ts0, ts1, tsg0, tsg1, tin1, tin0, time_hyb
    character(80) :: string

    ! Charge distance
    Real (8), Allocatable :: rhomtref(:,:,:) ! muffin-tin charge density (reference)
    Real (8), Allocatable :: rhoirref(:)     ! interstitial real-space charge density (reference)
    real(8) :: conv_old, conv_emp ! convergence for more empty orbitals 
    logical :: exist
!! TIME - Initialisation segment
    call timesec(tsg0)
    call timesec(ts0)

! Tetrahedron method
    input%groundstate%stypenumber = -1

! Initialise global variables
    Call timesec (tin0)
    Call init0
    Call timesec (tin1)
    time_init0 = tin1-tin0
    Call timesec (tin0)
    Call init1
    Call timesec (tin1)
    time_init1 = tin1-tin0
   ! Call init2

!-------------------
! print info
!-------------------

    if (rank==0) then
! open INFO.OUT file
        open(60, File='INFO'//trim(filext), Action='WRITE', Form='FORMATTED')
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
!
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

    if (rank==0) then
        write(string,'("* Hybrids module started")')
        call printbox(60,"*",string)
        call flushifc(60)
    end if

    !---------------------------------------
    ! initialise HF related input parameters
    !---------------------------------------
    time_hyb = 0.d0
    call timesec(ts0)
    call init_hybrids()
    call timesec(ts1)
    time_hyb = time_hyb+ts1-ts0

    ! reference density
    If (allocated(rhomtref)) deallocate(rhomtref)
    Allocate(rhomtref(lmmaxvr,nrmtmax,natmtot))
    If (allocated(rhoirref)) deallocate(rhoirref)
    Allocate (rhoirref(ngrtot))

    !----------------------
    ! restart is requested
    !----------------------
    select case (task)

      case(0)
        !_____________________________________________________________
        ! Initialization: PBE-DFT self-consistent run
        !
        if (rank==0) then
          write(60,*)
          write(60,*) 'Performing PBE self-consistent run'
        end if
        ex_coef = 0.d0
        ec_coef = 1.d0
        call scf_cycle(-1)
        ex_coef = input%groundstate%Hybrid%excoeff
        ec_coef = input%groundstate%Hybrid%eccoeff

        ! Kinetic energy of core states (is not going to be updated)
        call energykncr()

        ! save PBE potential and density
        string = filext
        filext='_PBE.OUT'
        call writestate()
        filext = string
        if (rank == 0) Then
          write(60,*) "Storing STATE_PBE.OUT"
          call flushifc(60)
        end if
        !________________________________
        ! Inizialize mixed product basis
        call init_product_basis()
      
      case(1)
write(*,*)"fromfile"
        !----------------------------------------------
        ! Initialization from previous hybrid SCF run
        !----------------------------------------------

        if (rank == 0) then
          write(60,*)
          write(60,*) 'Restart from previous hybrid calculations'
          call flushifc(60)
        end if
        !_______________________________________________________________________________________________
        ! step 1: read PBE potential to initialize radial parts of the LAPW basis and the product-basis
        inquire(File='STATE_PBE.OUT', Exist=exist)
        if (exist) then
          string = filext
          filext = '_PBE.OUT'
          call readstate()
          filext = string
        else
          stop 'ERROR(hybrids): Restart is not possible, STATE_PBE.OUT is missing!'
        end if
        !copy to vrelmt because STATE.OUT only has veffmt saved
        vrelmt=veffmt


        ! Core/Valence radial functions and integrals required for scf_cycle() + task=7
        ! (initialzed from PBE)
        call gencore()          ! generate the core wavefunctions and densities
        call linengy()          ! find the new linearization energies
        call genapwfr()         ! generate the APW radial functions
        call genlofr(.false.)   ! generate the local-orbital radial functions
        call olprad()           ! compute the overlap radial integrals
        call energykncr()       ! core kinetic energy
        !
        ! Inizialize mixed-product basis
        call init_product_basis()
        ! Initialize mt_hcsf parameter
        call MTInitAll(mt_hscf)
        !_____________________________________________________________________________________
        ! step 2: Read hybrid density and potential and get prepared for scf_cycle() + task=7
        call readstate()
        call readfermi()
        call hmlint(mt_hscf)
        call genmeffig()
        call MTNullify(mt_hscf)
      if (input%groundstate%do.eq."extraempty") then

        conv_emp = 0.d0
        conv_old = 0.d0
        splittfile = .False.
        task=7
        ex_coef = input%groundstate%Hybrid%excoeff
        ec_coef = input%groundstate%Hybrid%eccoeff
        input%groundstate%mixerswitch = 1
        input%groundstate%scfconv = 'charge'

        
        do ihyb=1, input%groundstate%Hybrid%maxscl
                write(*,*) "----------restart with orbitals----------", ihyb

                splittfile = .False.
                call MTInitAll(mt_hscf)
                call hmlint(mt_hscf)
                call calc_vxnl()
!                write(*,*)"afer calc_vxnl------"
                do ik = 1, 14
                        write(*,'(14F13.9)') dble(vxnl(ik,1:14, :))
                end do

                call calc_vnlmat()
                call mt_hscf%release
                call timesec(ts0)
                call scf_cycle(-1)
                call timesec(ts1)
                do ik = 1, nstfv
                     !  write(*,*)dble(vxnl(ik,ik, 1))
                       conv_emp = conv_emp+ dble(vxnl(ik,ik, 1))**2
                end do
                write(*,*)conv_emp-conv_old, "convergence criteria--------------"
                if ((conv_emp-conv_old).lt.1e-5) then
                        exit
                end if
                call flushifc(60)
                conv_old = conv_emp
                conv_emp = 0.d0

        end do

      end if 
      case default
        stop 'ERROR(hybrids): Not supported task!'

    end select
    !--------------------------------------------------
    ! External SCF loop
    !--------------------------------------------------

    ! hybrid-functional switch
    task = 7

    ex_coef = input%groundstate%Hybrid%excoeff
    ec_coef = input%groundstate%Hybrid%eccoeff

    ! Reference density
    rhomtref(:,:,:) = rhomt(:,:,:)
    rhoirref(:) = rhoir(:)

    ! settings for convergence and mixer
    input%groundstate%mixerswitch = 1
    input%groundstate%scfconv = 'charge'





    do ihyb = 1, input%groundstate%Hybrid%maxscl

      If (rank==0) Then
        write(string,'("+ Hybrids iteration number : ", I4)') ihyb
        call printbox(60,"+",string)
        if (rank==0) write(60,*)
        Call flushifc(60)
      End If

      ! to insure correct reading from GS files
      splittfile = .False.


      call MTInitAll(mt_hscf)
      call hmlint(mt_hscf)


      !-----------------------------------
      ! Calculate the non-local potential
      !-----------------------------------
      call timesec(ts0)
      call calc_vxnl()
write(*,*)"afer calc_vxnl------"
!                do ik = 1, 14
 !                       write(*,'(14F13.9)') dble(vxnl(ik,1:14, :))
  !              end do

      call timesec(ts1)

      if ((input%groundstate%outputlevelnumber>1) .and.rank==0) then
        write(60,*)
        write(60,'(" CPU time for vxnl (seconds)",T45 ": ", F12.2)') ts1-ts0
      end if
      !------------------------------------------
      call timesec(ts0)
      call calc_vnlmat()
      call timesec(ts1)
      If ((input%groundstate%outputlevelnumber>1) .and.rank==0) Then
        write(60,'(" CPU time for vnlmat (seconds)",T45 ": ", F12.2)') ts1-ts0
      end if
      time_hyb = time_hyb+ts1-ts0
      call mt_hscf%release

      !--------------------------------------------------
      ! Internal SCF loop
      !--------------------------------------------------
      call timesec(ts0)
      call scf_cycle(-1)
      call timesec(ts1)

      If ((input%groundstate%outputlevelnumber>1) .and.rank==0) Then
        write(60, '(" CPU time for scf_cycle (seconds)",T45 ": ", F12.2)') ts1-ts0
        write(60,*)
      end if
      if (rank == 0) then
          call writeengy(60)
          write(60,*)
          write(60,'(" DOS at Fermi energy (states/Ha/cell)",T45 ": ", F18.8)') fermidos
          call writechg(60,input%groundstate%outputlevelnumber)
          if (fermidos<1.0d-4) call printbandgap(60)
          call flushifc(60)
      end if

      !---------------------------
      ! Convergence check
      !---------------------------
      call chgdist(rhomtref, rhoirref)
      if (rank==0) Then
        write(60,*)
        write(60,'(" Charge distance                   (target) : ",G13.6," (",G13.6,")")') &
        &     chgdst, input%groundstate%epschg
      end if
      if (chgdst < input%groundstate%epschg) then
        if (rank==0) Then
            write(string,'("Convergence target is reached")')
            call printbox(60,"+",string)
            call flushifc(60)
            call write_vxnl()
        end if
        exit ! exit ihyb-loop
      else
        ! update the reference
        rhomtref(:,:,:) = rhomt(:,:,:)
        rhoirref(:) = rhoir(:)
      end if
      et = engytot

      ! maximum iterations is reached
      If (ihyb == input%groundstate%Hybrid%maxscl) Then
          If (rank==0) Then
              write(string,'("Reached hybrids self-consistent loops maximum : ", I4)') &
              &  input%groundstate%Hybrid%maxscl
              call printbox(60,"+",string)
              call warning('Warning(hybrids): Reached self-consistent loops maximum')
              Call flushifc(60)
          End If
      End If
    end do ! ihyb

    if (rank==0) then
        write(string,'("+ Hybrids module stopped")')
        call printbox(60,"+",string)
        call flushifc(60)
    end if

!-------------------------------------!
!   output timing information
!-------------------------------------!
! TIME - Fifth IO segment
    call timesec(ts0)
    if (rank==0) then
!____________________
! Close support files
!
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

! TIME - End of fifth IO segment
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
! Clean data
!----------------------------------------
    close(600) ! HYBRIDS.OUT
    call delete_core_states
    call delete_product_basis
    call exit_hybrids
    nullify(input%gw)
    call rereadinput

    Return
End Subroutine
!EOC
