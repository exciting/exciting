
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

    implicit none
    logical :: reducek_
    integer :: lmax, ik
    real(8) :: t0, t1, tstart, tend
    integer :: stype_

    call timesec(tstart)
    
    ! Importantly, current implementation uses exclusively 
    ! "Extended linear tetrahedron method for the calculation of q-dependent
    !  dynamical response functions", to be published in Comp. Phys. Commun. (2010)
    stype_ = input%groundstate%stypenumber
    input%groundstate%stypenumber = -1
    task = 1

    ! It is important at this stage to switch off the update of the effective potential
    ! (critical for OEP/HYBRIDS related methods)
    input%groundstate%xctypenumber = 1
    xctype(1) = 1

    ! initialize global exciting variables with GW parameters (see parse_gwinput.f90)
    call init0
    call init1    

    if (hybridhf) then
      
      isreadstate0 = .false.
      filext = '_PBE.OUT'
      call readstate()
      filext = '.OUT'

    else
      
      isreadstate0 = .true.
      filext = "_GW.OUT"
      
      call timesec(t0)
      if (.not. input%gw%skipgnd) then
        input%groundstate%maxscl = 1
        call scf_cycle(-2)
        if (rank == 0) then
          ! safely remove unnecessary files
          call filedel('EIGVAL'//trim(filext))
          call filedel('LINENGY'//trim(filext))
          call filedel('EVALCORE'//trim(filext))
          call filedel('OCCSV'//trim(filext))
          call filedel('BROYDEN.OUT')
          call writefermi
          ! call writeeval
        end if
        call barrier
      end if
      
      call readstate()

    end if
    
    ! if BSE is used just after GW, libbzint
    ! may create the segmentation faults when trying to use a general
    ! k-point shift vkloff
    input%groundstate%stypenumber = stype_
        
    call timesec(t1)
    time_initscf = time_initscf+t1-t0
    
    ! Generate the k/q-point meshes
    call timesec(t0)
    call init_kqpoint_set
    call timesec(t1)
    time_initkpt = time_initkpt+t1-t0

    ! Intialize auxiliary arrays used further for convenience
    call init_misc_gw
      
    ! Mixed basis initialization
    call timesec(t0)
    call init_product_basis
    call timesec(t1)
    time_initmb = time_initmb+t1-t0
        
    ! Get Kohn-Sham eigenvalues
    call timesec(t0)
    call init_dft_eigenvalues
    call timesec(t1)
    time_initeval = time_initeval+t1-t0
    
    ! Frequency grid initialization
    call timesec(t0)

    if (input%gw%taskname=='g0w0' .or. &
        input%gw%taskname=='gw0'  .or. &
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
    
    ! timing
    call timesec(tend)
    time_initgw = time_initgw+tend-tstart

    return
end subroutine init_gw
