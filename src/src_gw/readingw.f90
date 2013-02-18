!BOP
!
! !ROUTINE: readingw
!
! !INTERFACE:
      subroutine readingw()
      
!
! !DESCRIPTION:
!
! This subroutine reads the main input file for the gw calculation

! !USES:
 
      use modinput
      use modmain
      use modgw
 
! !LOCAL VARIABLES:      

      implicit none
      
      character(8) :: test
      character(6) :: fdep
      character(6) :: fgrid
      character(3) :: corflag
      real(8)      :: len
      
! !EXTERNAL ROUTINES: 

!EOP
!BOC
      call boxmsg(fgw,'*',"GW input parameters")
!
!     Define the task to calculate
!      
      test=input%gw%taskname
      write(fgw,*)
      call linmsg(fgw,'-','')
      write(fgw,*)  ' TASK:', test
      call linmsg(fgw,'-','')
      write(fgw,*)

      select case (test)
!
        case('GW','gw')
          testid=0
!      
        case('LAPW','lapw')
          testid=1
!      
        case('EVEC','evec')
          testid=2
!
        case('PROD','prod')
          testid=3
!
        case('MIXF','mixf')
          testid=4
!
        case('COMP','comp')
          testid=5
!
        case('COUL','coul')
          testid=6
!
        case('EMAC','emac')
          testid=7
!
        case('ACON','acon')
          testid=9
!
        case('DFEV','epsev')
          testid=10
!        
        case('WEV','wev')
          testid=11
!          
        case('SEPL','sepl')
          testid=12
!
        case('VXC','vxc')
          testid=13
!
        case('EX','ex')
          testid=14
!
!        case('BAND','band')
!          testid=16
!          
        case('ROTMAT','rotmat')
          testid=18   
!               
        case default
          write(fgw,*)'ERROR: Wrong test option!! Valid options are:'  
          write(fgw,*)'lapw - Calculate LAPW basis functions for plotting'
          write(fgw,*)'evec - Calculate DFT eigenvectors for plotting'
          write(fgw,*)'prod - Calculate products of eigenvectors for plotting'
          write(fgw,*)'mixf - Calculate products of eigenvectors and &
         &      their expansion in the mixed basis for plotting'
          write(fgw,*)'comp - Test completeness of the mixed basis'
          write(fgw,*)'coul - Test bare Coulomb potential'
          write(fgw,*)'emac - Calculate the DFT macroscopic dielectric function'
          write(fgw,*)'epsev - Calculate eigenvalues of the dielectric matrix' 
          !write(fgw,*)'epgw - Calculate the GW macroscopic dielectric function'
          !write(fgw,*)'wev  - Calculate eigenvalues of the screened Coulomb potential'
          write(fgw,*)'acon - Perform only the analytic continuation'
          write(fgw,*)'sepl - Plot Selfenergy as a function of frequency'
          write(fgw,*)'vxc  - Calculate the matrix elements of the DFT exchange-correlation potential'
          write(fgw,*)'ex   - Hartree-Fock calculation (exchange only)'
          !write(fgw,*)'band - Calculate the bandstructure'
          write(fgw,*)'rotmat - Calculate and check the MB rotation matrices (symmetry feature)'
          write(fgw,*)'gw   - Performs one complete GW cycle'
          write(fgw,*)
      end select  
      call linmsg(fgw,'*',"end readingw")
!
!     Read the BZ convolution method
!           
      select case (input%gw%bzconv)
        case ('sum','SUM')
          convflg=0
        case ('tetra','TETRA')
          convflg=1
        case default
          write(fgw,*) 'Warning: Wrong BZ convolution option!! Valid options are:'
          write(fgw,*) 'tetra: Use the linearized tetrahedrom method'
          write(fgw,*) 'sum: Simple sum over k-points'
          write(fgw,*) 'Taking default value: tetra'
          convflg=1
      end select
!
!     Current implementation is only using imaginary frequencies (see FHI-gap for more info)
!
      fdep='imfreq'
      fflg=3

      write(fgw,*)'BZ Convolution method:'
      write(fgw,*) trim(input%gw%bzconv), '   ', fdep

      call linmsg(fgw,'-','')
!      
!     Read the data for the frequecy grid
!      
      if (associated(input%gw%FreqGrid)) then
        fgrid=input%gw%FreqGrid%fgrid
        nomeg=input%gw%FreqGrid%nomeg
        freqmax=input%gw%FreqGrid%freqmax
      else
        fgrid='gaule2'
        nomeg=16
        freqmax=0.42
      end if
      select case (fgrid)
        case('eqdist','EQDIST')
         wflag=1
        case('gaulag','GAULAG')
         wflag=2
        case('gaule2','GAULE2')
         wflag=3
        case('gauleg','GAULEG')
         wflag=4
        case default
          write(fgw,*) 'Warning: Wrong frequency grid option!! Valid options are:'
          write(fgw,*) 'eqdist - equidistant frequencies from 0 &
     & to freqmax (no integration weights)'
          write(fgw,*) 'gaulag - grid for Gauss-Laguerre &
     & quadrature from 0 to infinity'
          write(fgw,*) 'gauleg - grid for Gauss-Legendre &
     & quadrature from 0 to freqmax'
          write(fgw,*) 'gaule2 - grid for Gauss-Legendre &
     & quadrature from 0 to freqmax and from freqmax to infinity'
          write(fgw,*)'Taking default value: gaule2'
          fgrid='gaule2'
          wflag=3
      end select
      write(fgw,*) 'Frequency grid: fgrid, nomeg, freqmax'
      write(fgw,*) fgrid, nomeg, freqmax
      
      call linmsg(fgw,'-','')

!     Options for the analytic continuation
      if(associated(input%gw%SelfEnergy))then
         iopes=input%gw%SelfEnergy%iopes
         iopac=input%gw%SelfEnergy%iopac
         npol=input%gw%SelfEnergy%npol
      else
         ! default values
         npol=2
         iopes=0
         iopac=1
      end if

      write(fgw,*) '- Option for calculating selfenergy (iopes): ', iopes
      select case (iopes)
        case(0)
          write(fgw,*) "  -- perturbative G0W0 without energy shift"
        case(1)
          write(fgw,*) "  -- perturbative G0W0 with energy shift"
        case(2)
          write(fgw,*) "  -- iterative G0W0 with energy shift"
        case(3)
          write(fgw,*) "  -- iterative G0W0 without energy shift"
        case default
          call errmsg(1,'readingw','Wrong value for iopes')
      end select

      write(fgw,*) '- Type of analytic continuation (iopac): ', iopac
      select case (iopac)
        case(1) 
          write(fgw,*) "  -- RGN method(Rojas, Godby and Needs)"
        case(2)
          write(fgw,*) "  -- Pade's approximation"
        case default
          call errmsg(1,'readingw','Wrong value for iopac')
      end select

      if(npol.eq.0) then
        write(fgw,*) "The input npol == 0: use default npol setup"
        if(iopac.eq.1) then
          npol=2
        else
          npol=nomeg/2
        endif
      endif
      write(fgw,*) '- Nr. of poles used in analytic continuation: ', npol

      call linmsg(fgw,'-','')
!
!     Number of the empty bands used in GW code
!
      if ((input%gw%nempty<1).or.(input%gw%nempty>nmatmax)) then
        input%gw%nempty=input%groundstate%nempty
      else
        input%groundstate%nempty=input%gw%nempty
        nstfv=int(chgval/2.d0)+input%gw%nempty+1
        nstsv=nstfv*nspinor
      end if
      write(fgw,*)'Number of empty states (groundstate): ', input%groundstate%nempty
      write(fgw,*)'Number of empty states (GW): ', input%gw%nempty
      write(fgw,*)
!
!     new number of first-variational states (to be used in GW)
!
!     Band interval where GW corections are applied
!
      ibgw=input%gw%ibgw
      if ((ibgw<1).or.(ibgw>nstfv)) then
        ibgw=1
      end if
      nbgw=input%gw%nbgw
      if ((nbgw<1).or.(nbgw>nstfv)) then
        ! use just 30 empty states for calculating QP states
        nbgw=min(nstfv,int(chgval/2.d0)+30)
      end if
      write(fgw,'(a,2i7)') ' GW output band range: ', ibgw, nbgw
      write(fgw,*)
!
!     Read the options for the mixed basis functions
!
      write(fgw,*) 'Mixed basis parameters:'
      write(fgw,*) '- Interstitial:'
      write(fgw,*) '  -- maximum |G| of IPW in gmaxvr units (gmb):', input%gw%MixBasis%gmb
      write(fgw,*) '- MT-Spheres:'
      write(fgw,*) '  -- l_max (lmaxmb): ', input%gw%MixBasis%lmaxmb
      write(fgw,*) '  -- linear dependence tolerance (epsmb): ', input%gw%MixBasis%epsmb
      call linmsg(fgw,'-','')
!      
!     Read the parameters for the Bare coulomb potential
!      
      if(associated(input%gw%BareCoul)) then
        pwm = input%gw%BareCoul%pwm
        stctol = input%gw%BareCoul%stctol
        barcevtol=input%gw%BareCoul%barcevtol
      else
        pwm=2.0
        stctol=1.0d-10
        barcevtol=-1.0d-10
      endif 
      write(fgw,*) 'Bare Coulomb parameters:'
      write(fgw,*) 'Maximum |G| in gmaxvr*gmb units:', pwm
      write(fgw,*) 'Error tolerance for struct. const.:', stctol
      write(fgw,*) 'Tolerance to choose basis functions from bare Coulomb matrix eigenvectors: ', barcevtol

      call linmsg(fgw,'-','')
!
!     Read pmat (the momentum matrix elements from file)
!
      write(fgw,*) 'Read options for the momentum matrix elements:'
      write(fgw,*) 'read pmat = ', input%gw%rpmat

      call linmsg(fgw,'-','')
!
!     Read epsilon (the dielectric function matrix elements is read from file)
!
      write(fgw,*) 'Read options for the dielectric function matrix elements:'
      write(fgw,*) 'read reps = ', input%gw%reps
            
      call linmsg(fgw,'-','')    
!
      select case (input%gw%coreflag)
        case('all','ALL')
          iopcore=0
          write(fgw,*)' all: All electron calculation'
        case('xal','XAL')
          iopcore=1
          write(fgw,*)' xal: all electron for exchange, valence only for correlation'
        case('val')
          iopcore=2
          write(fgw,*)' val: Valence electrons only'
        case('vab')
          iopcore=3
          write(fgw,*)' vab: Valence-only without core states in mixbasis'
       end select
!
      call linmsg(fgw,'-','')
!
!     To take into account anisotropy of \epsilon for q->0
!      
      !q0_eps=(/1.d0,1.d0,1.d0/)/sqrt(3.0d0)
      q0_eps=input%gw%q0eps
     
!
!     K/Q point grids
!     
      if ((input%gw%ngridq(1).gt.0).and. &
     &    (input%gw%ngridq(2).gt.0).and. &
          (input%gw%ngridq(3).gt.0)) then
        if ((input%gw%ngridq(1).ne.input%groundstate%ngridk(1)).or. &
       &    (input%gw%ngridq(2).ne.input%groundstate%ngridk(2)).or. &
            (input%gw%ngridq(3).ne.input%groundstate%ngridk(3))) then
          write(fgw,*)
          write(fgw,*)'GW calculations performed on the k/q-grid:'
          write(fgw,*)'ngridq = ', input%gw%ngridq
          input%groundstate%ngridk=input%gw%ngridq
        end if
      end if

      len=input%gw%vqloff(1)**2+ &
     &    input%gw%vqloff(2)**2+ &
     &    input%gw%vqloff(3)**2
      if (len.gt.1.0d-10) then
        write(fgw,*)
        write(fgw,*)'Attention! k/q-point shift is specified!'
        write(fgw,*)
        write(fgw,*)'k/q-shift: vqloff = ', input%gw%vqloff
        input%groundstate%vkloff=input%gw%vqloff
      end if
      
      return
      end subroutine readingw
      
!EOC      
