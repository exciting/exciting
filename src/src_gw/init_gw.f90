
subroutine init_gw()

    use modinput
    use modmain
    use modgw
    use mod_gaunt_coefficients
    use modmpi
    use modxs, only: isreadstate0
    use mod_hybrids, only: hybridhf
    use m_filedel
    use mod_hdf5
    use gw_scf, only: set_gs_solver_threads, thread_consistent_scf

    implicit none
    logical :: reducek
    integer :: lmax, ik
    real(8) :: t0, t1, tstart, tend

    !> Threads to use for GW call to ground state SCF
    integer :: gs_solver_threads

    call timesec(tstart)

    ! Importantly, current implementation uses exclusively
    ! "Extended linear tetrahedron method for the calculation of q-dependent
    !  dynamical response functions", to be published in Comp. Phys. Commun. (2010)
    input%groundstate%stypenumber = -1
    task = 1

    ! It is important at this stage to switch off the update of the effective potential
    ! (critical for OEP/HYBRIDS related methods)
    input%groundstate%xctypenumber = 1
    xctype(1) = 1

    ! initialize global exciting variables with GW parameters (see parse_gwinput.f90)
    call init0()

    call timesec(t0)

    if (hybridhf) then

      ! In this case we are forced to use the same parameters as for the underlying
      ! hybrid functional, i.e., nempty and ngridk.
      ! Remark: the wavefunction for a general k-point will be obtained applying
      ! a rotation algorithm implemented in getevecfv.f90.
      ! Unfortunately, there are known some small artifacts caused by the rotation.
      call init1()

      isreadstate0 = .false.
      filext = '_PBE.OUT'
      call readstate()
      filext = '.OUT'

    else

      ! to preserve the WF phase one recompute the eigenfunctions for the entire BZ
      reducek = input%groundstate%reducek
      input%groundstate%reducek = .false.
      call init1()
      input%groundstate%reducek = reducek

      isreadstate0 = .true.
      filext = "_GW.OUT"

      if (.not. input%gw%skipgnd) then
        input%groundstate%maxscl = 1
        gs_solver_threads = set_gs_solver_threads(mpiglobal)
        call thread_consistent_scf(gs_solver_threads)

        if (rank == 0) then
          ! safely remove unnecessary files
          call filedel('EIGVAL'//trim(filext))
          call filedel('LINENGY'//trim(filext))
          call filedel('EVALCORE'//trim(filext))
          call filedel('OCCSV'//trim(filext))
          call filedel('EFERMI'//trim(filext))
          call filedel('BROYDEN.OUT')
          call writefermi
        end if
        call barrier
      end if

      call readstate()

    end if

    call timesec(t1)
    time_initscf = time_initscf+t1-t0

    ! Generate the k/q-point meshes
    call timesec(t0)
    call init_kqpoint_set()
    call timesec(t1)
    time_initkpt = time_initkpt+t1-t0

    ! Intialize auxiliary arrays used further for convenience
    call init_misc_gw()

    ! Mixed basis initialization
    call timesec(t0)
    call init_product_basis()
    ! Print the product-basis info
    if (rank==0) then
      call boxmsg(fgw,'-',"Mixed product WF info")
      write(fgw,*) ' Maximal number of MT wavefunctions per atom: ', lmixmax
      write(fgw,*) ' Total number of MT wavefunctions:            ', locmatsiz
      write(fgw,*) ' Maximal number of PW wavefunctions:          ', Gqset%ngkmax
      write(fgw,*) ' Total number of mixed-product wavefunctions: ', matsizmax
      write(fgw,*)
    end if
    call timesec(t1)
    time_initmb = time_initmb+t1-t0

    ! Frequency grid initialization
    call timesec(t0)

    if (input%gw%taskname=='g0w0' .or. &
        input%gw%taskname=='emac') then
      call generate_freqgrid(freq, &
      &                      input%gw%freqgrid%fgrid, &
      &                      input%gw%freqgrid%fconv, &
      &                      input%gw%freqgrid%nomeg, &
      &                      input%gw%freqgrid%freqmin, &
      &                      input%gw%freqgrid%freqmax)
      if (rank==0) call print_freqgrid(freq,fgw)
    else
      ! frequency independent method
      freq%nomeg = 1
      freq%freqmin = 0.d0
      freq%freqmax = 0.d0
      freq%fconv = 'imfreq'
      allocate(freq%freqs(freq%nomeg))
      freq%freqs(1) = 1.d-3
      allocate(freq%womeg(freq%nomeg))
      freq%womeg(1) = 1.d0
    end if

    call timesec(t1)
    time_initfreq = time_initfreq+t1-t0

    ! Gaunt coefficients
    lmax = max(input%groundstate%lmaxapw+1, &
    &          2*(input%gw%mixbasis%lmaxmb+1))
    call calcgauntcoef(lmax)

    ! Get Kohn-Sham eigenvalues
    call timesec(t0)
    call init_dft_eigenvalues()
    call timesec(t1)
    time_initeval = time_initeval+t1-t0

    ! timing
    call timesec(tend)
    time_initgw = time_initgw+tend-tstart

    return
end subroutine init_gw
