
subroutine parse_gwinput()

    use modinput
    use modmain
    use modgw
    use mod_coulomb_potential, only: vccut
    use modmpi
    use mod_hybrids, only: hybridhf
    implicit none
    integer :: idum
    real(8) :: rdum

    if (associated(input%groundstate%spin)) then
        if (rank==0) call boxmsg(fgw,'!',"WARNING! GW for magnetic materials is not yet implemented!")
        if (ldapu /= 0) then
            if (rank==0) call boxmsg(fgw,'!','Spin-polarized GW@LDA+U version is not implemented!')
            call barrier()
            stop
        end if
    end if

    if (rank==0) call boxmsg(fgw,'=',"GW input parameters")

!-------------------------------------------------------------------------------
! Debugging mode
!-------------------------------------------------------------------------------
    if (input%gw%debug) then
        if (rank==0) call boxmsg(fgw,'-',"Debug mode is activated!")
        if (rank>0) then
            if (rank==0) then
                write(fgw,*) 'WARNING(parse_gwinput):'
                write(fgw,*) '    Debug option is not designed to be used in parallel ...'
                write(fgw,*)
            end if
        end if
        input%gw%debug = input%gw%debug.and.(rank==0)
    end if

!-------------------------------------------------------------------------------
! Task name parser
!-------------------------------------------------------------------------------
    input%gw%taskname = trim(input%gw%taskname)
    if (rank==0) then
      write(fgw,*)
      write(fgw,*) 'GW taskname:'
      write(fgw,*)
    end if
    select case(input%gw%taskname)
        case('skip')
        case('g0w0')
            if (rank==0) write(fgw,*) '  g0w0 - G0W0 run'
        case('g0w0-x')
            if (rank==0) write(fgw,*) '  g0w0-x - Exchange only G0W0 run'
        case('cohsex')
            if (rank==0) write(fgw,*) '  cohsex - Coulomb hole plus screened exchange approximation (COHSEX)'
        case('band','band2')
            if (rank==0) write(fgw,*) '  band - Calculate QP bandstructure'
        case('dos')
            if (rank==0) write(fgw,*) '  dos - Calculate QP DOS'
        case('emac')
            if (rank==0) write(fgw,*) '  emac - Calculate the DFT frequency-dependent macroscopic dielectric function for q=0'
        case('vxc')
            if (rank==0) write(fgw,*) '  vxc  - Calculate the matrix elements of the DFT exchange-correlation potential'
        case('pmat')
            if (rank==0) write(fgw,*) '  pmat  - Calculate the matrix elements of the momentum operator'
        case('evalqp')
            if (rank==0) write(fgw,*) '  evalqp - Solve the quasiparticle equation'
        ! case('eps_r')
        !     if (rank==0) write(fgw,*) '  eps_r - Test only option'
        ! case('chi0_r')
        !     if (rank==0) write(fgw,*) '  chi0_r - Calculate polarizability matrix elements for all q-points and convert them to real-space'
        ! case('chi0_q')
        !     if (rank==0) write(fgw,*) '  chi0_q - Calculate polarizability matrix elements for a given q-point'
        case('kqgen')
            if (rank==0) write(fgw,*) '  kqgen - (testing option) Test generation of k/q-point grids'
        case('lapw')
            if (rank==0) write(fgw,*) '  lapw - (testing option) Calculate LAPW basis functions for plotting'
        case('evec')
            if (rank==0) write(fgw,*) '  evec - (testing option) Calculate DFT eigenvectors for plotting'
        case('prod')
            if (rank==0) write(fgw,*) '  prod - (testing option) Calculate products of eigenvectors for plotting'
        case('mixf')
            if (rank==0) write(fgw,*) '  mixf - (testing option) Calculate products of eigenvectors and'
            if (rank==0) write(fgw,*) '         their expansion in the mixed basis for plotting'
        case('comp')
            if (rank==0) write(*,*)   '  comp - (testing option) Test completeness of the mixed basis'
        case('test_aaa')
            if (rank==0) write(fgw,*) '  Test AAA interpolation'
        case('band_specfunc')
            if (rank==0) write(fgw,*) '  Compute spectral function along the k-path'
        case('sv')
            if (rank==0) write(fgw,*) '  Apply second variation procedure'
        case default
            if (rank==0) write(*,*) 'ERROR(parse_gwinput): Wrong task name!'
            if (rank==0) write(*,*)
            if (rank==0) write(*,*) '  Specified value: taskname = ', trim(input%gw%taskname)
            if (rank==0) write(*,*)
            call barrier()
            stop
    end select
    if (rank==0) call linmsg(fgw,'-','')

!-------------------------------------------------------------------------------
! Frequency integration parameters
!-------------------------------------------------------------------------------
    if (.not.associated(input%gw%freqgrid)) &
     &  input%gw%freqgrid => getstructfreqgrid(emptynode)

    if (input%gw%taskname=='g0w0' .or. &
    &   input%gw%taskname=='emac') then
      ! no frequencies is required
      if (rank==0) write(fgw,*) 'Frequency integration parameters:'
      if (rank==0) write(fgw,*) 'Number of frequencies: ', input%gw%freqgrid%nomeg
      if (rank==0) write(fgw,*) 'Cutoff frequency: ', input%gw%freqgrid%freqmax
      if (rank==0) write(fgw,*) 'Grid type:'
      select case (input%gw%freqgrid%fgrid)
        case('eqdist')
          if (rank==0) write(fgw,*) '  eqdist - Equaly spaced mesh (for tests purposes only)'
        case('gaulag')
          if (rank==0) write(fgw,*) '  gaulag - Grid for Gauss-Laguerre quadrature'
        case('gauleg')
          if (rank==0) write(fgw,*) '  gauleg - Grid for Gauss-Legendre quadrature, from 0 to freqmax'
        case('gauleg2')
          if (rank==0) write(fgw,*) '  gauleg2 - Double Gauss-Legendre grid: [0, freqmax] + [freqmax, infty]'
        case('clencurt2')
          if (rank==0) write(fgw,*) '  clencurt2 - Semi-infinite Clenshaw-Curtis grid'
      end select
      if (rank==0) write(fgw,*) 'Convolution method:'
      select case (input%gw%freqgrid%fconv)
        case('nofreq')
          if (rank==0) write(fgw,*) '  nofreq : no frequency dependence of the weights'
        case('refreq')
          if (rank==0) write(fgw,*) '  refreq : weights calculated for real frequencies'
        case('imfreq')
          if (rank==0) write(fgw,*) '  imfreq : weights calculated for imaginary frequencies'
        case default
          if (rank==0) write(*,*) 'ERROR(parse_gwinput): Unknown frequency convolution method!'
          if (rank==0) write(*,*) '  Currently supported options are:'
          if (rank==0) write(*,*) '  nofreq : no frequency dependence of the weights'
          if (rank==0) write(*,*) '  refreq : weights calculated for real frequencies'
          if (rank==0) write(*,*) '  imfreq : weights calculated for imaginary frequencies'
          stop
      end select
      if (rank==0) call linmsg(fgw,'-','')
    end if

!-------------------------------------------------------------------------------
! Analytic continuation parameters
!-------------------------------------------------------------------------------
    if (.not.associated(input%gw%selfenergy)) &
    &  input%gw%selfenergy => getstructselfenergy(emptynode)
    if (rank==0) write(fgw,*) 'Correlation self-energy parameters:'
    if (input%gw%selfenergy%nempty>0) then
      if (rank==0) write(fgw,*) 'Number of empty states:', input%gw%selfenergy%nempty
    end if

    if (rank==0) write(fgw,*) 'Solution of the QP equation:'
    select case (input%gw%selfenergy%eqpsolver)
        case(0)
            if (rank==0) write(fgw,*) "  0 - perturbative solution"
        case(1)
            if (rank==0) write(fgw,*) "  1 - Z=1 solution"
        case(2)
            if (rank==0) write(fgw,*) "  2 - iterative solution"
        case default
            if (rank==0) write(*,*) 'ERROR(parse_gwinput): Illegal value for input%gw%SelfEnergy%eqpsolver'
            stop
    end select

    if (rank==0) write(fgw,*) 'Energy alignment:'
    select case (input%gw%selfenergy%eshift)
        case(0)
            if (rank==0) write(fgw,*) "  0 - no alignment"
        case(1)
            if (rank==0) write(fgw,*) "  1 - self-consistency at the Fermi level (iterative)"
        case(2)
            if (rank==0) write(fgw,*) "  2 - self-consistency at the Fermi level (perturbative)"
        case default
            if (rank==0) write(*,*) 'ERROR(parse_gwinput): Illegal value for input%gw%SelfEnergy%eshift'
            stop
    end select

    if (rank==0) write(fgw,*) 'Analytic continuation method:'
    select case (trim(input%gw%selfenergy%actype))
        case('pade','Pade','PADE')
            if (rank==0) then
                write(fgw, *) " PADE - Thiele's reciprocal difference method"
                write(fgw, *) " Citation: H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977)"
            endif     
        case('aaa','AAA')
            if (rank==0) then
                write(fgw, *) " AAA - The AAA algorithm for rational approximation"
                write(fgw, *) " Citation: Y. Nakatsukasa, O. Sete, L. N. Trefethen, SIAM J. Sci. Comp. 40 (2018), A1494-A1522"
            endif 
        case default
            if (rank==0) write(*, *) 'ERROR(parse_gwinput): Illegal value for input%gw%SelfEnergy%actype'
            stop
    end select

    if (rank==0) write(fgw, *) 'Scheme to treat singularities:'
    select case (trim(input%gw%selfenergy%singularity))
      case('none')
        if (rank==0) then 
            write(fgw, *) ' No scheme is used (test purpose only)'
        endif 
      case('avg')
        if (rank==0) then
            write(fgw, *) ' Replace the singular term by the corresponding spherical average over small volume around Gamma point.'
        endif     
      case('mpb')
        if (rank==0) then
            write(fgw, *) ' Auxiliary function method "mpb"'
            write(fgw, *) ' Citation: S. Massidda, M. Posternak, and A. Baldereschi, PRB 48, 5058 (1993)'
        endif
      case('crg')
        if (rank==0) then
            write(fgw, *) ' Auxiliary function method "crg"'
            write(fgw, *) ' Citation: P. Carrier, S. Rohra, and A. Goerling, PRB 75, 205126 (2007)'
        endif 
      case('rim')
        if (rank==0) then
            write(fgw, *) 'RIM (experimental):'
            write(fgw, *) 'Citation: Yambo code'
        endif 
      case default
        write(*,*) 'ERROR(parse_gwinput): Unknown singularity treatment scheme!'
        stop
    end select
    if (rank==0) call linmsg(fgw,'-','')

!-------------------------------------------------------------------------------
! Product basis parameters
!-------------------------------------------------------------------------------
    if (.not.associated(input%gw%mixbasis)) &
    &  input%gw%MixBasis => getstructmixbasis(emptynode)
    if (rank==0) write(fgw,*) 'Mixed product basis parameters:'
    if (input%gw%MixBasis%lmaxmb<0) then
        if (rank==0) write(*,*) 'ERROR(parser_gwinput): Illegal value of input%gw%MixBasis%lmaxmb'
        call barrier()
        stop
    end if
    if (rank==0) write(fgw,*) '  MT part:'
    if (rank==0) write(fgw,*) '    Angular momentum cutoff: ', input%gw%MixBasis%lmaxmb
    if (rank==0) write(fgw,*) '    Linear dependence tolerance factor: ', input%gw%MixBasis%epsmb
    if (rank==0) write(fgw,*) '  Interstitial:'
    if (rank==0) write(fgw,*) '    Plane wave cutoff (in units of Gkmax): ', input%gw%MixBasis%gmb
    if (rank==0) call linmsg(fgw,'-','')

!-------------------------------------------------------------------------------
! Bare Coulomb potential parameters
!-------------------------------------------------------------------------------
    if (.not.associated(input%gw%barecoul)) input%gw%barecoul => getstructbarecoul(emptynode)
    if (rank==0) then 
        write(fgw, *) 'Bare Coulomb potential parameters:'
        write(fgw, '(X, A, X, F17.14)') '  Plane wave cutoff (in units of Gkmax*input%gw%MixBasis%gmb): ', input%gw%barecoul%pwm
        write(fgw, *) '  Error tolerance for structure constants: ', input%gw%barecoul%stctol
        write(fgw, *) '  Tolerance factor to reduce the MB size based on'
        write(fgw, *) '  the eigenvectors of the bare Coulomb potential: ', input%gw%barecoul%barcevtol
    endif

    select case (trim(input%gw%barecoul%cutofftype))
      case('none')
        vccut = .false.
      case('0d')
        vccut = .true.
        if (rank==0) write(fgw,*) '  Spherical (0d) cutoff is applied'
      case('1d')
        vccut = .true.
        if (rank==0) write(fgw,*) '  Wired (1d) cutoff is applied (1d periodicity along z-axis)'
      case('2d')
        vccut = .true.
        if (rank==0) write(fgw,*) '  Slab (2d) cutoff is applied (vacuum along z-axis)'
      case default
        if (rank==0) write(*,*) 'ERROR(parse_gwinput): Illegal value for input%gw%barecoul%cutofftype'
        if (rank==0) write(*,*) '  Currently supported options are:'
        if (rank==0) write(*,*) '  none - No cutoff (default)'
        if (rank==0) write(*,*) '  0d   - Spherical (0d) cutoff'
        if (rank==0) write(*,*) '  1d   - Wired (1d) cutoff (periodicity along z-axis)'
        if (rank==0) write(*,*) '  2d   - Slab geometry (vacuum along z-axis)'
        call barrier()
        stop
    end select
    if (vccut) then
        ! Coulomb potential truncation techniques are implemented only for the PW basis
        input%gw%barecoul%basis = "pw"
        ! TODO(Alex/Ronaldo) Issue #142. Default should be taken from schema
        input%gw%barecoul%pwm = 4.d0
    end if
    if (rank==0) call linmsg(fgw,'-','')

    ! TODO(Alex/Ronaldo) Issue #142 . Implement the fixes required to get GW in PW basis consistent
    ! with exciting Nitrogen.
    if (input%gw%barecoul%basis == "pw") then
        call terminate_mpi_env(mpiglobal, message="GW in a plane-wave basis is currently unsupported.")
    end if

!-------------------------------------------------------------------------------
! Parameters for averaging the dielectric function
!-------------------------------------------------------------------------------
    if (.not.associated(input%gw%scrcoul)) &
    &  input%gw%scrcoul => getstructscrcoul(emptynode)
    if (rank==0) write(fgw,*) 'Screened Coulomb potential:'
    select case(trim(input%gw%scrcoul%scrtype))
      case('rpa','RPA')
        if (rank==0) write(fgw,*) '  Full-frequency Random-Phase Approximation'
      case('ppm','PPM')
        if (rank==0) write(fgw,*) '  Godby-Needs plasmon-pole model'
        if (rank==0) write(fgw,*) '  Plasmon frequency: ', input%gw%scrcoul%omegap
      case default
        if (rank==0) write(*,*) 'ERROR(parse_gwinput): Illegal value for input%gw%scrcoul%scrtype'
        if (rank==0) write(*,*) '  The supported options are:'
        if (rank==0) write(*,*) '  rpa - Full-frequency Random-Phase Approximation (default)'
        if (rank==0) write(*,*) '  ppm - Godby-Needs plasmon-pole model'
        call barrier()
        stop
    end select
    if (rank==0) call linmsg(fgw,'-','')

!-------------------------------------------------------------------------------
! Core electrons treatment
!-------------------------------------------------------------------------------
    if (rank==0) write(fgw,*) 'Core electrons treatment:'
    select case (input%gw%coreflag)
        case('all')
            if (rank==0) write(fgw,*) '  all - Core states are included in all calculations'
        case('xal')
            if (rank==0) write(fgw,*) '  xal - Core states are included in exchange but not in correlation'
        case('val')
            if (rank==0) write(fgw,*) '  val - Core states are excluded in all calculations, but kept in &
           & the construction of mixed basis'
        case('vab')
            if (rank==0) write(fgw,*) '  vab - Core states are excluded completely'
        case default
            if (rank==0) write(*,*) 'ERROR(parse_gwinput): Unknown type of core electrons treatment!'
            if (rank==0) write(*,*) '  Currently supported options are:'
            if (rank==0) write(*,*) '    all - Core states are included in all calculations'
            if (rank==0) write(*,*) '    xal - Core states are included in exchange but not in correlation'
            if (rank==0) write(*,*) '    val - Core states are excluded in all calculations, but kept in &
           & the construction of mixed basis'
            if (rank==0) write(*,*) '    vab - Core states are excluded completely'
            call barrier()
            stop
    end select
    if (rank==0) call linmsg(fgw,'-','')

!-------------------------------------------------------------------------------
! Special treatment in case of hybrid functionals
!-------------------------------------------------------------------------------
    hybridhf = .false.
    if (xctype(1) >= 400) then
        hybridhf = .true.
        if (.not. associated(input%groundstate%Hybrid)) &
            input%groundstate%Hybrid => getstructHybrid(emptynode)
        ! use parameters from groundstate/hybrids for consistency
        input%gw%nempty = input%groundstate%nempty
        input%gw%ngridq = input%groundstate%ngridk
        input%gw%vqloff = input%groundstate%vkloff
    end if

    !-------------------------------------------------------------------------------
    ! Band range where GW corrections are applied
    !-------------------------------------------------------------------------------
    if (isspinorb()) then
        input%gw%ibgw = 1
        input%gw%nbgw = int(chgval/2.d0) + input%gw%nempty + 1
    end if
    ibgw = input%gw%ibgw
    nbgw = input%gw%nbgw
    if (nbgw < 1) nbgw = input%gw%nempty
    if (ibgw >= nbgw) then
        if (rank==0) write(*,*) 'ERROR(parse_gwinput): Illegal values for ibgw ot nbgw!'
        if (rank==0) write(*,*) '    ibgw = ', ibgw, '   nbgw = ', nbgw
        call barrier()
        stop
    end if
    if (rank==0) write(fgw,'(a,2i7)') ' Interval of quasiparticle states (ibgw, nbgw): ', ibgw, nbgw
    if (rank==0) write(fgw,*)

!-------------------------------------------------------------------------------
! Number of the empty bands used in GW code
!-------------------------------------------------------------------------------
    if (input%gw%nempty < 1) then
        if (rank==0) write(fgw,*) 'WARNING(parse_gwinput): Number of empty states is not specified!'
        if (rank==0) write(fgw,*) '  This parameter must be carefully chosen based on the convergence tests'
        if (rank==0) write(fgw,*)
        input%gw%nempty = input%groundstate%nempty
    end if
    ! overwrite the GS value to be able to run scf_cycle()
    input%groundstate%nempty = max(input%gw%nempty,input%gw%selfenergy%nempty)
    if (rank==0) write(fgw,*)'Number of empty states (GW): ', input%gw%nempty
    if (rank==0) write(fgw,*)

!-------------------------------------------------------------------------------
! k/q point grids
!-------------------------------------------------------------------------------
    idum = input%gw%ngridq(1)*input%gw%ngridq(1)+ &
    &      input%gw%ngridq(2)*input%gw%ngridq(2)+ &
    &      input%gw%ngridq(3)*input%gw%ngridq(3)
    if (idum == 0) then
        if (rank==0) write(fgw,*) 'WARNING(parse_gwinput): Number of k/q-points is not specified!'
        if (rank==0) write(fgw,*) '  This parameter has a crucial influence on the results and'
        if (rank==0) write(fgw,*) '  must be carefully chosen based on the convergence tests.'
        if (rank==0) write(fgw,*) '  Too large values could make GW calculations very time consuming.'
        if (rank==0) write(fgw,*)
        input%gw%ngridq = input%groundstate%ngridk
    else
        if ((input%gw%ngridq(1)<=0).or. &
        &   (input%gw%ngridq(2)<=0).or. &
        &   (input%gw%ngridq(3)<=0)) then
            if (rank==0) write(fgw,*) 'ERROR(parse_gwinput): Illegal value for k/q-points grid!'
            call barrier()
            stop
        end if
    end if
    if (rank==0) write(fgw,*) 'k/q-points grid: ', input%gw%ngridq
    input%groundstate%ngridk = input%gw%ngridq

    rdum = input%gw%vqloff(1)**2 + &
           input%gw%vqloff(2)**2 + &
           input%gw%vqloff(3)**2
    if (rdum > 1.d-8) then
        if (rank==0) write(fgw,*)'k/q-shift: ', input%gw%vqloff
        input%groundstate%vkloff = input%gw%vqloff
    end if

    !-------------------------------------------------------------------------------
    ! Matrix block size
    !-------------------------------------------------------------------------------
    mblksiz = input%gw%mblksiz
    if (mblksiz > 0) then
        if (rank==0) write(fgw,*)
        if (rank==0) write(fgw,*) ' Matrix block size: ', mblksiz
    else
        mblksiz = 1000000 ! just a big number to account for all available states
    end if

    if (rank==0) call flushifc(fgw)

end subroutine
