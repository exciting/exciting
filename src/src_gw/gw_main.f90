
subroutine gw_main()

    use modinput
    use modmain
    use modgw
    use modmpi
    use mod_mpi_gw
    use m_getunit
    use mod_hdf5
    use mod_aaa_approximant

    implicit none
    real(8) :: tstart, tend

    !--------------------------------------------
    ! Skip initializing and running any GW task
    !--------------------------------------------
    if (trim(input%gw%taskname)=='skip') return

    ! Initialize timing and memory usage variables
    call timesec(tstart)
    call init_timing()

    !---------------------
    ! Main GW output file
    !---------------------
    if (rank == 0) then
        call getunit(fgw)
        open(fgw, File='GW_INFO.OUT')
        ! open(fgw, File='GW_INFO.OUT', Access='Append')
        if (input%gw%debug) then
            call getunit(fdebug)
            open(fdebug, File='debug.info', Action='Write')
        end if
    end if

    !----------------------------------------------
    ! initialize GW MPI environment
    !----------------------------------------------
    call init_mpi_gw

    !-----------------------------------------------------------
    ! Parse and check the validity of some GW input parameters
    !-----------------------------------------------------------
    call init0() ! to allow access to some global data such as charges, species, etc.
    call parse_gwinput()

    !----------------
    ! Task selector
    !----------------
    select case(input%gw%taskname)

        ! GW calculations
        case('g0w0','g0w0-x','cohsex')
            call task_gw()

        ! Calculate the QP band structure
        case('band')
            if (rank==0) call task_band()

        ! Calculate QP DOS
        case('dos')
            if (rank==0) call task_dos()

        ! Calculate the macroscopic dielectric function
        case('emac')
            call task_emac()

        ! Calculate diagonal matrix elements of the exchange-correlation potential
        case('vxc')
            call init_gw()
            call calcvxcnn()

        ! Calculate matrix elements of the momentum operator
        case('pmat')
            call init_gw()
            call calcpmatgw()

        ! Perform analytic continuation of the correlation self-energy and calculate QP energies
        case('evalqp')
            if (rank==0) call task_evalqp()

        ! Visualize the spectral function along the bandstructure path
        case('band_specfunc')
            call task_band_specfunc()

        ! Apply second-variational procedure to GW results to include SOC
        case('sv')
            input%gw%skipgnd = .True.
            call init_gw()
            if (rank==0) call task_second_variation()

        ! (for testing only) Check generation of k-, q-, k+G, q+G, etc. sets
        case('kqgen')
            if (rank==0) call test_kqpts()

        ! (for testing only) Calculate LAPW basis functions for plotting
        case('lapw')
            call init_gw()
            if (rank==0) call plot_lapw()

        ! (for testing only) Calculate LAPW eigenvectors for plotting (test option)
        case('evec')
            call init_gw()
            if (rank==0) call plot_evec()

        ! (for testing only) Calculate LAPW eigenvectors products for plotting (test option)
        case('prod')
            call init_gw()
            if (rank==0) call test_prodfun()

        ! (for testing only) Calculate eigenvectors products and compare them with mixed-product basis expansion
        case('mixf')
            call init_gw()
            if (rank==0) call test_mixfun()

        ! (for testing only) Integrate eigenvector products directly and
        ! as a sum of the Minm matrix elements
        case('comp')
            call init_gw()
            if (rank==0) call test_mixcomp()

        ! Calculate and store the (q,\omega)-dependent dielectric function
            ! case('epsilon')
        !     call task_epsilon()

            ! Compute and output q-dependent \epsilon_00 along a k-path
            ! case('emac_q')
            !     call task_emac_q

        ! Compute and output the polarizability in the real space
        ! case('chi0_r')
        !     call task_chi0_r

        ! Compute and output the polarizability in the reciprocal space
            ! case('chi0_q')
        !     call task_chi0_q

        ! Compute and output the dielectric function in the real space
            ! case('eps_r')
        !     call task_eps_r

        ! (for testing only) Test the AAA analytical continuation technique
        ! case('test_aaa')
        !     if (rank==0) call test_aaa_1()
        !     if (rank==0) call test_aaa_2()

    end select

    ! output timing info
    call timesec(tend)
    time_total = time_total+tend-tstart

    ! output some timing info
    call print_timing()

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
