!BOP
!
! !ROUTINE: readingw
!
! !INTERFACE:
subroutine parse_gwinput
!
! !DESCRIPTION:
!
! This subroutine check important GW input parameters

! !USES:
 
    use modinput
    use modmain
    use modgw
    use modmpi
    implicit none
 
! !LOCAL VARIABLES:
    integer :: idum
    real(8) :: rdum
      
!EOP
!BOC

    if (rank==0) call boxmsg(fgw,'*',"GW input parameters")
    
!-------------------------------------------------------------------------------
! Debugging mode
!-------------------------------------------------------------------------------
    if (input%gw%debug) then
        if (rank==0) write(fgw,*) 'The code run in debugging mode'
        if (rank>0) then
            if (rank==0) write(fgw,*) 'WARNING(parse_gwinput): Debug option is not supposed to &
           & be used in parallel ...'
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
        case('g0w0_x')
            if (rank==0) write(fgw,*) '  g0w0_x - Exchange only G0W0 run'
        case('cohsex')
            if (rank==0) write(fgw,*) '  cohsex - Coulomb hole plus screened exchange approximation (COHSEX)'
        case('gw0')
            if (rank==0) write(fgw,*) '  gw0  - GW0 self-consistent run'
        case('band','band2')
            if (rank==0) write(fgw,*) '  band - Calculate QP bandstructure'
        case('dos')
            if (rank==0) write(fgw,*) '  dos - Calculate QP DOS'
        case('emac')
            if (rank==0) write(fgw,*) '  emac - Calculate the DFT frequency-dependent macroscopic dielectric function for q=0'
        case('emac_q')
            if (rank==0) write(fgw,*) '  emac_q - Calculate the DFT q-dependent macroscopic dielectric function for \omega=0'
        case('vxc')
            if (rank==0) write(fgw,*) '  vxc  - Calculate the matrix elements of the DFT exchange-correlation potential'
        case('pmat')
            if (rank==0) write(fgw,*) '  pmat  - Calculate the matrix elements of the momentum operator'    
        case('acon')
            if (rank==0) write(fgw,*) '  acon - Perform only the analytic continuation of &
            &the correlation self energy and recalculate QP energies'
        case('epsilon')
            if (rank==0) write(fgw,*) '  epsilon - Calculate dielectric function matrix elements'
        case('eps_r')
            if (rank==0) write(fgw,*) '  eps_r - Test only option'
        case('chi0_r')
            if (rank==0) write(fgw,*) '  chi0_r - Calculate polarizability matrix elements for all q-points and convert to real-space'
        case('chi0_q')
            if (rank==0) write(fgw,*) '  chi0_q - Calculate polarizability matrix elements for a given q-point'
!       case('epsev')
!           if (rank==0) write(fgw,*) '  epsev - Calculate eigenvalues of the dielectric matrix' 
!       case('wev')
!           if (rank==0) write(fgw,*) '  wev  - Calculate eigenvalues of the screened Coulomb potential'
!       case('epsgw')
!           if (rank==0) write(fgw,*) '  epsgw - Calculate the GW macroscopic dielectric function'
        case('kqgen')
            if (rank==0) write(fgw,*) '  kqgen - (testing option) Test generation of k/q-point grids'
        case('bzintw')
            if (rank==0) write(fgw,*) '  bzintw - (testing option) Test generation of k- k/q-dependent BZ integration weights'
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
            if (rank==0) write(fgw,*) '  comp - (testing option) Test completeness of the mixed basis'
        case('coul')
            if (rank==0) write(fgw,*) '  coul - (testing option) Test bare Coulomb potential'
        case('sepl')
            if (rank==0) write(fgw,*) '  sepl - (testing option) Plot Selfenergy as a function of frequency'
        case('rotmat')
            if (rank==0) write(fgw,*) '  rotmat - (testing option) Calculate and check the MB rotation matrices (symmetry feature)'
        case('wannier')
            if (rank==0) write(fgw,*) '  wannoer - (testing option) Wannier-interpolate QP-energies'
        
        case default
            if (rank==0) write(*,*) 'ERROR(parse_gwinput): Wrong task name!'
            if (rank==0) write(*,*)
            if (rank==0) write(*,*) 'Specified value: taskname = ', trim(input%gw%taskname)
            if (rank==0) write(*,*)
            if (rank==0) write(*,*) 'Currently supported options are'
            if (rank==0) write(*,*) '  skip - Skip GW part execution'
            if (rank==0) write(*,*) '  g0w0   - Perform G0W0 run'
            if (rank==0) write(*,*) '  gw0   - Perform GW0 self-consistent run'
            if (rank==0) write(*,*) '  g0w0_x - Exchange only G0W0 run'
            if (rank==0) write(*,*) '  cohsex - Coulomb hole plus screened exchange approximation (COHSEX)'
            if (rank==0) write(*,*) '  band - Calculate QP bandstructure'
            if (rank==0) write(*,*) '  emac - Calculate the DFT macroscopic dielectric function'
            if (rank==0) write(*,*) '  vxc  - Calculate the matrix elements of the DFT exchange-correlation potential'
            if (rank==0) write(*,*) '  pmat - Calculate the matrix elements of the momentum operator'
            if (rank==0) write(*,*) '  acon - Perform only the analytic continuation of the correlation self energy and &
            &  recalculate QP energies'
            if (rank==0) write(*,*) '  epsev - Calculate eigenvalues of the dielectric matrix' 
            !if (rank==0) write(fgw,*) '  epgw - Calculate the GW macroscopic dielectric function'
            !if (rank==0) write(fgw,*) '  wev  - Calculate eigenvalues of the screened Coulomb potential'
            if (rank==0) write(*,*) '  lapw - (test option) Calculate LAPW basis functions for plotting'
            if (rank==0) write(*,*) '  evec - (test option) Calculate DFT eigenvectors for plotting'
            if (rank==0) write(*,*) '  prod - (test option) Calculate products of eigenvectors for plotting'
            if (rank==0) write(*,*) '  mixf - (test option) Calculate products of eigenvectors and &
            &  their expansion in the mixed basis for plotting'
            if (rank==0) write(*,*) '  comp - (test option) Test completeness of the mixed basis'
            if (rank==0) write(*,*) '  coul - (test option) Test bare Coulomb potential'
            if (rank==0) write(*,*) '  sepl - (test option) Plot Selfenergy as a function of frequency'
            if (rank==0) write(*,*) '  rotmat - (test option) Calculate and check the MB rotation matrices (symmetry feature)'
            if (rank==0) write(*,*) '  kqgen - (test option) Test generation of k/q-point grids'
            if (rank==0) write(*,*)
            stop
    end select  
    if (rank==0) call linmsg(fgw,'-','')

!-------------------------------------------------------------------------------
! Frequency integration parameters    
!-------------------------------------------------------------------------------
    if (.not.associated(input%gw%freqgrid)) &
     &  input%gw%freqgrid => getstructfreqgrid(emptynode)
     
    if (input%gw%taskname=='g0w0' .or. &
    &   input%gw%taskname=='gw0'  .or. &
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
        case('gaule2')
          if (rank==0) write(fgw,*) '  gaule2 - Grid for double Gauss-Legendre quadrature,'
          if (rank==0) write(fgw,*) '           from 0 to freqmax and from freqmax to infinity'
        case('gauleg')
          if (rank==0) write(fgw,*) '  gauleg - Grid for Gauss-Legendre quadrature, from 0 to freqmax'
        case default
          if (rank==0) write(*,*) 'ERROR(parse_gwinput): Unknown frequency grid type!'
          stop
      end select
      if (rank==0) write(fgw,*) 'Convolution method:'
      select case (input%gw%freqgrid%fconv)
        case('nofreq')
          if (rank==0) write(fgw,*) '  nofreq : no frequecy dependence of the weights'
        case('refreq')
          if (rank==0) write(fgw,*) '  refreq : weights calculated for real frequecies'
        case('imfreq')
          if (rank==0) write(fgw,*) '  imfreq : weights calculated for imaginary frequecies'
        case default
          if (rank==0) write(*,*) 'ERROR(parse_gwinput): Unknown frequency convolution method!'
          if (rank==0) write(*,*) '  Currently supported options are:'
          if (rank==0) write(*,*) '  nofreq : no frequecy dependence of the weights'
          if (rank==0) write(*,*) '  refreq : weights calculated for real frequecies'
          if (rank==0) write(*,*) '  imfreq : weights calculated for imaginary frequecies'
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
    select case (input%gw%selfenergy%iopes)
        case(-2)
            if (rank==0) write(fgw,*) " -2 - perturbative G0W0 without renormalization Z (testing)"
        case(0)
            if (rank==0) write(fgw,*) "  0 - perturbative G0W0 without energy shift"
        case(1)
            if (rank==0) write(fgw,*) "  1 - perturbative G0W0 with energy shift"
        case(2)
            if (rank==0) write(fgw,*) "  2 - iterative G0W0 with energy shift"
        case(3)
            if (rank==0) write(fgw,*) "  3 - iterative G0W0 without energy shift"
        case default
            if (rank==0) write(*,*) 'ERROR(parse_gwinput): Illegal value for input%gw%SelfEnergy%iopes'
            if (rank==0) write(*,*) '  Currently supported options are:'
            if (rank==0) write(*,*) '  0 - perturbative G0W0 without energy shift'
            if (rank==0) write(*,*) '  1 - perturbative G0W0 with energy shift'
            if (rank==0) write(*,*) '  2 - iterative G0W0 with energy shift'
            if (rank==0) write(*,*) '  3 - iterative G0W0 without energy shift'
            stop        
    end select
    if (rank==0) write(fgw,*) 'Analytic continuation method:'
    select case (trim(input%gw%selfenergy%actype))
        case('pade','Pade','PADE')
            iopac = 0
            if (rank==0) write(fgw,*) " pade - Thiele's reciprocal difference method &
            &(by H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179 (1977))"
        case('mpf','MPF')
            iopac = 1
            if (rank==0) write(fgw,*) " mpf - Multi-pole fitting (by Rojas, Godby and Needs PRL 74, 1827 (1995))"
        case default
            if (rank==0) write(*,*) 'ERROR(parse_gwinput): Illegal value for input%gw%SelfEnergy%actype'
            if (rank==0) write(*,*) '  Currently supported options are:'
            if (rank==0) write(*,*) "  pade - Thiele's reciprocal difference method &
            &(by Vidberg and Serence, J. Low Temp. Phys. 29, 179 (1977))"
            if (rank==0) write(*,*) "  mpf  - Multi-pole fitting (by Rojas, Godby and Needs (PRL 74, 1827 (1995))"
    end select
    if (input%gw%selfenergy%npol==0) then
        select case (trim(input%gw%selfenergy%actype))
        case('pade','Pade','PADE')
            input%gw%selfenergy%npol = input%gw%freqgrid%nomeg/2
        case('mpf','MPF')
            input%gw%selfenergy%npol = 2
        end select
    else
        if (input%gw%selfenergy%npol>input%gw%freqgrid%nomeg) then
            if (rank==0) then
                write(*,*) 'ERROR(parse_gwinput): the number of poles cannot exceed the number of frequencies!'
                stop
            end if
        end if
    endif
    if (rank==0) write(fgw,*) 'Number of poles used in analytic continuation:', input%gw%selfenergy%npol
    if (rank==0) write(fgw,*) 'Scheme to treat singularities:'
    select case (trim(input%gw%selfenergy%singularity))
      case('none')
        if (rank==0) write(fgw,*) 'No scheme is used (test purpose only)'
      case('mpb')
        if (rank==0) write(fgw,*) 'Auxiliary function method by &
        &S. Massidda, M. Posternak, and A. Baldereschi, PRB 48, 5058 (1993)'
      case('crg')  
        if (rank==0) write(fgw,*) 'Auxiliary function method by &
        &P. Carrier, S. Rohra, and A. Goerling, PRB 75, 205126 (2007)'
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
    if (.not.associated(input%gw%barecoul)) &
    &  input%gw%barecoul => getstructbarecoul(emptynode)
    if (rank==0) write(fgw,*) 'Bare Coulomb potential parameters:'
    if (rank==0) write(fgw,*) '  Plane wave cutoff (in units of Gkmax*input%gw%MixBasis%gmb): ', &
    &  input%gw%barecoul%pwm
    if (rank==0) write(fgw,*) '  Error tolerance for structure constants: ', &
    &  input%gw%barecoul%stctol
    if (rank==0) write(fgw,*) '  Tolerance factor to reduce the MB size based on'
    if (rank==0) write(fgw,*) '  the eigenvectors of the bare Coulomb potential: ', input%gw%barecoul%barcevtol
    select case (trim(input%gw%barecoul%cutofftype))
      case('none')
        vccut = .false.
      case('0d')
        vccut = .true.
        if (rank==0) write(fgw,*) '  Spherical (0d) cutoff is applied'
      case('2d')
        vccut = .true.
        if (rank==0) write(fgw,*) '  Slab (2d) cutoff is applied (vacuum along z-axis)'
      case default
        if (rank==0) write(*,*) 'ERROR(parse_gwinput): Illegal value for input%gw%barecoul%cutofftype'
        if (rank==0) write(*,*) '  Currently supported options are:'
        if (rank==0) write(*,*) '  none - No cutoff (default)'
        if (rank==0) write(*,*) '  0d   - Spherical (0d) cutoff'
        if (rank==0) write(*,*) '  2d   - Slab geometry (vacuum along z-axis)'
        stop
    end select
    if (rank==0) call linmsg(fgw,'-','')
    
!-------------------------------------------------------------------------------
! Parameters for averaging the dielectric function
!-------------------------------------------------------------------------------     
    if (.not.associated(input%gw%scrcoul)) &
    &  input%gw%scrcoul => getstructscrcoul(emptynode)
    if (rank==0) write(fgw,*) 'Screened Coulomb potential parameters:'
    if (rank==0) write(fgw,*) '  Type of averaging: ', trim(input%gw%scrcoul%sciavtype)
    select case (trim(input%gw%scrcoul%sciavtype))
      case('isotropic')
        if (rank==0) write(fgw,'(a,3f8.4)') '  Averaging direction: ', input%gw%scrcoul%q0eps
      case('anisotropic')
        if (rank==0) write(fgw,*) '  Angular momentum cutoff: ', input%gw%scrcoul%lmaxdielt
        if (rank==0) write(fgw,*) '  Number of points used for the Lebedev-Laikov grid: ', input%gw%scrcoul%nleblaik
        if (rank==0) write(fgw,*) '  Averaging of the body of the dielectric function: ', input%gw%scrcoul%sciavbd
      case default
        if (rank==0) write(*,*) 'ERROR(parse_gwinput): Illegal value for input%gw%scrcoul%sciavtype'
        if (rank==0) write(*,*) '  Currently supported options are:'
        if (rank==0) write(*,*) '  isotropic - Simple averaging along the direction given in q0eps'
        if (rank==0) write(*,*) '  anisotropic - Anisotropic screening following Freysold et al.'
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
            stop
    end select
    if (rank==0) call linmsg(fgw,'-','')    
    
!-------------------------------------------------------------------------------
! Number of the empty bands used in GW code
!-------------------------------------------------------------------------------
    if (input%gw%nempty<1) then
        if (rank==0) write(fgw,*) 'WARNING(parse_gwinput): Number of empty states is not specified!'
        if (rank==0) write(fgw,*) '  This parameter must be carefully chosen based on the convergence tests'
        if (rank==0) write(fgw,*) '  Too large values make GW calculations very time consuming'
        if (rank==0) write(fgw,*)
        if (rank==0) write(fgw,*) '  Used default (small) value for input%gw%nempty'
        input%gw%nempty = 10
    end if
    input%groundstate%nempty = max(input%gw%nempty,input%gw%selfenergy%nempty)
    if (rank==0) write(fgw,*)'Number of empty states (GW): ', input%gw%nempty
        
!-------------------------------------------------------------------------------
! Band range where GW corrections are applied
!-------------------------------------------------------------------------------

! lower QP band index
    ibgw = input%gw%ibgw

! upper QP band index
    nbgw = input%gw%nbgw
    if (nbgw < 1) nbgw = input%groundstate%nempty

    if ( ibgw >= nbgw ) then
        if (rank==0) write(*,*) 'ERROR(parse_gwinput): Illegal values for ibgw ot nbgw!'
        if (rank==0) write(*,*) '    ibgw = ', ibgw, '   nbgw = ', nbgw
        stop
    end if

    if (rank==0) write(fgw,'(a,2i7)') ' Quasiparticle band range: ', ibgw, nbgw
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
        if (rank==0) write(fgw,*) '  Too large values make GW calculations very time consuming.'
        if (rank==0) write(fgw,*)
        if (rank==0) write(fgw,*) '  Set the default value for input%gw%ngridq:' 
        input%gw%ngridq = (/2, 2, 2/)
    else
        if ((input%gw%ngridq(1)<=0).or. &
        &   (input%gw%ngridq(2)<=0).or. &
        &   (input%gw%ngridq(3)<=0)) then
            if (rank==0) write(fgw,*) 'ERROR(parse_gwinput): Illegal value for k/q-points grid!'
            stop
        end if
    end if
    if (rank==0) write(fgw,*) '  k/q-points grid: ', input%gw%ngridq
    input%groundstate%ngridk = input%gw%ngridq

    rdum = input%gw%vqloff(1)**2+ &
    &      input%gw%vqloff(2)**2+ &
    &      input%gw%vqloff(3)**2
    if (rdum > 1.d-8) then
        if (rank==0) write(fgw,*)'Attention! k/q-point shift is specified!'
        if (rank==0) write(fgw,*)'k/q-shift: ', input%gw%vqloff
        input%groundstate%vkloff = input%gw%vqloff
    end if
    
    call flushifc(fgw)
      
    return
end subroutine
!EOC      
