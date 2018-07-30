
subroutine gw_main()

    use modinput
    use modmain
    use modgw
    use modmpi
    use mod_mpi_gw
    use m_getunit
    use mod_hdf5
    
    implicit none
    character(80) :: fname
    real(8) :: tstart, tend

    !--------------------------------------------
    ! Skip initializing and running any GW task
    !--------------------------------------------
    if (trim(input%gw%taskname)=='skip') return
    
    ! Initialize timing and memory usage variables
    call timesec(tstart)
    call init_timing    


    !---------------------
    ! Main GW output file
    !---------------------
    if (rank==0) then
        call getunit(fgw)
        open(fgw,File='GW_INFO.OUT')
        if (input%gw%debug) then
            call getunit(fdebug)
            open(fdebug,File='debug.info',Action='Write')
        end if
        call boxmsg(fgw,'=','Main GW output file')
        call flushifc(fgw)
    end if
    
    !----------------------------------------------
    ! initialize GW MPI environment
    !----------------------------------------------
    call init_mpi_gw
    
    !-----------------------------------------------------------
    ! Parse and check the validity of some GW input parameters
    !-----------------------------------------------------------
    call parse_gwinput
    
    !----------------------------------------------
    ! Store all important results to the hdf5 file
    !----------------------------------------------
#ifdef _HDF5_
    call hdf5_initialize()
    fgwh5 = "gw_output.h5"
    if (rank==0) then
      select case (input%gw%taskname)
      
        case('acon','band', 'sepl')
          continue
        
        case default
          call hdf5_create_file(fgwh5)
          call hdf5_create_group(fgwh5,"/","parameters")
          call hdf5_create_group(fgwh5,"/","kpoints")
          call write_gw_parameters_hdf5
          
      end select
    end if
#endif
    
    !----------------
    ! Task selector
    !----------------
    select case(input%gw%taskname)
    
        ! GW calculations
        case('g0w0','g0w0_x','gw0','cohsex')
            call task_gw
            
        ! Calculate the QP band structure
        case('band')
            if (rank==0) call task_band
            
        ! Calculate QP DOS
        case('dos')
            if (rank==0) call task_dos
            
        ! Calculate the macroscopic dielectric function
        case('emac')
            call task_emac
        
        ! test option: q-dependent \epsilon along a k-path
        ! case('emac_q')
        !     call task_emac_q

        ! Calculate diagonal matrix elements of the exchange-correlation potential
        case('vxc')
            call init_gw
            call calcvxcnn
            
        ! Calculate matrix elements of the momentum operator
        case('pmat')
            call init_gw
            call calcpmatgw
            
        ! Perform analytic continuation of the correlation self-energy and 
        ! calculate QP energies
        case('acon')
            if (rank==0) call task_analytic_continuation
        
        ! Calculate and store the (q,\omega)-dependent dielectric function
        case('epsilon')
            call task_epsilon

        ! Calculate and store the (q,\omega)-dependent dielectric function
        case('chi0_r')
            call task_chi0_r
            
        ! Calculate and store the (q,\omega)-dependent dielectric function
        case('chi0_q')
            call task_chi0_q

        ! Calculate and store the (q,\omega)-dependent dielectric function
        case('eps_r')
            call task_eps_r
                        
        ! Calculate the eigenvalues the LDA dielectric function and its inverse
        ! case('epsev')
        !    call task_epsev

        ! Calculate the eigenvalues the GW dielectric function and its inverse
        ! case('epsgw') 
        !    call task_epsgw

        ! Calculate the eigenvalues of the screened coulomb potential
        ! case('wev') 
        !    call task_wev

        ! (testing option) Check generation of k-, q-, k+G, q+G, etc. sets 
        case('kqgen')
            if (rank==0) call test_kqpts

        ! (testing option) Test generation of k- k/q-dependent BZ integration weights
        case('bzintw')
            if (rank==0) call test_bzintw

        ! (testing option) Calculate LAPW basis functions for plotting
        case('lapw')
            call init_gw
            if (rank==0) call plot_lapw

        ! (testing option) Calculate LAPW eigenvectors for plotting (test option)
        case('evec')
            call init_gw
            if (rank==0) call plot_evec

        ! (testing option) Calculate LAPW eigenvectors products for plotting (test option)
        case('prod') 
            call init_gw
            if (rank==0) call test_prodfun

        ! (testing option) Calculate eigenvectors products compared 
        ! with mix basis expansion for plotting
        case('mixf')
            call init_gw
            if (rank==0) call test_mixfun
            
        ! (testing option) Integrate eigenvector products directly and 
        ! as a sum of the Minm matrix elements
        case('comp')
            call init_gw
            if (rank==0) call test_mixcomp

        ! (testing option) Test the bare coulomb matrix for various q-points
        ! case('coul')
        !    if (rank==0) call test_coulpot

        ! (testing option) Plot the self-energy
        case('sepl') 
            if (rank==0) call plot_selfenergy

        case('wannier')
            !call init_gw
            if( associated( input%properties%wannier)) then
              input%properties%wannier%input = "gw"
              !input%properties%wannier%fst = ibgw
              !input%properties%wannier%lst = nbgw
              call wannierlauncher
              !call task_eph ()
              call task_band 
            end if

        ! (testing option) Check the rotational matrix for MB functions
        ! case('rotmat')
            ! if (rank==0) call test_mbrotmat

    end select
    
    ! output timing info
    call timesec(tend)
    time_total = time_total+tend-tstart
    call print_timing
      
    ! close GW output file
    if (rank==0) then
      close(fgw)
      if (input%gw%debug) close(fdebug)
    end if
    
#ifdef _HDF5_
    call hdf5_finalize()
#endif

    return
end subroutine
