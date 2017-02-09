
subroutine init_gw()

    use modinput
    use modmain
    use modgw
    use mod_gaunt_coefficients
    use modmpi
    use modxs,      only: isreadstate0
    use m_filedel
    use mod_hdf5

    implicit none
    character(256) :: filext_save
    logical :: reducek_
    integer :: lmax
    real(8) :: t0, t1, tstart, tend
    integer :: stype_

    call timesec(tstart)
    
    !spinpol = associated(input%groundstate%spin)
    !if (spinpol) then
    !    write(*,*) 'GW EMERGENCY STOP!'
    !    stop 'Spin polarization is not yet implemented'
    !end if
    
    ! Importantly, current implementation uses exclusively 
    ! "Extended linear tetrahedron method for the calculation of q-dependent
    !  dynamical response functions", to be published in Comp. Phys. Commun. (2010)
    stype_ = input%groundstate%stypenumber
    input%groundstate%stypenumber = -1

    ! Extra call of the groundstate part to generate the data consisted with 
    ! the specified GW parameters
    call timesec(t0)
    if (.not.input%gw%skipgnd) then
        
        ! One needs to regenerate eigenvectors for the complete (non-reduced) k-point set
        reducek_ = input%groundstate%reducek
        input%groundstate%reducek = .false.
        
        ! Resume scf KS calculations and diagonalize one more time the Hamiltonian
        ! for new set of k- and nempty parameters.
        ! It is important at this stage to switch off the update of the effective potential
        ! (critical for OEP/HYBRIDS related methods)
        input%groundstate%xctypenumber = 1
        xctype(1) = 1
        
        ! initialize global exciting variables with (new) GW definitions
        call init0
        call init1
        
        task = 1
        input%groundstate%maxscl = 1
        filext_save = trim(filext)
        filext = "_GW.OUT"
        isreadstate0 = .true.
        call scf_cycle(-2)
        
        !-----------------------------------------------------
        ! rearrange the way to store eigenvalues and -vectors
        !-----------------------------------------------------
        if (rank==0) then
          call write_dft_orbitals
          ! safely remove unnecessary files
          call filedel('LINENGY'//trim(filext))
          call filedel('EVALFV'//trim(filext))
          call filedel('EVALSV'//trim(filext))
          call filedel('EVECFV'//trim(filext))
          call filedel('EVECSV'//trim(filext))
          call filedel('EVALCORE'//trim(filext))
          call filedel('OCCSV'//trim(filext))
          call filedel('BROYDEN.OUT')
          !call writeeval
          !call writefermi
        end if
        call barrier

        ! restore the initial value
        input%groundstate%reducek = reducek_
        filext = trim(filext_save)
        
    end if
    
    ! reinitialize global exciting variables
    call init0
    call init1
    
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

    ! Get the self-consistent effective + exchange-correlation potential
    call timesec(t0)
    call readstate
    call timesec(t1)
    time_io = time_io+t1-t0
      
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
    &   input%gw%taskname=='gw0'  .or. &
    &   input%gw%taskname=='emac') then
      call generate_freqgrid(freq, &
      &                      input%gw%freqgrid%fgrid, &
      &                      input%gw%freqgrid%fconv, &
      &                      input%gw%freqgrid%nomeg, &
      &                      input%gw%freqgrid%freqmax)
      if (rank==0) call print_freqgrid(freq,fgw)
#ifdef _HDF5_      
      if (rank==0) then
        call hdf5_write(fgwh5,"/parameters/freqgrid","freqs", &
        &               freq%freqs(1),(/freq%nomeg/))
        call hdf5_write(fgwh5,"/parameters/freqgrid","womeg", &
        &               freq%womeg(1),(/freq%nomeg/))
      end if
#endif

    else
      ! frequency independent method
      freq%nomeg = 1
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
    
    !------------------------------
    ! print the memory usage info
    !------------------------------
    call print_memory_usage

    ! timing
    call timesec(tend)
    time_initgw = time_initgw+tend-tstart

    return
    
contains

    subroutine print_memory_usage
        implicit none
        real(8) :: msize_tot
        if (rank==0) then
          call boxmsg(fgw,'-','Peak memory estimate (Mb, per process):')
          msize_tot = 0.d0
          !---------------------------------------
          ! kset
          msize = sizeof(kset)*b2mb
          !write(fgw,'(" kset:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! Gset
          msize = sizeof(Gset)*b2mb
          !write(fgw,'(" Gset:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! Gkset
          msize = sizeof(Gkset)*b2mb
          !write(fgw,'(" Gkset:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! kqset
          msize = sizeof(kqset)*b2mb
          !write(fgw,'(" kqset:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! Gqset
          msize = sizeof(Gqset)*b2mb
          !write(fgw,'(" Gqset:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          !---------------------------------------
          msize = sizeof(freq)*b2mb
          !write(fgw,'(" freq:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          !---------------------------------------
          ! umix
          msize = sizeof(umix)*b2mb
          !write(fgw,'(" umix:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! rtl
          msize = sizeof(rtl)*b2mb
          !write(fgw,'(" rtl:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! rrint
          msize = sizeof(rrint)*b2mb
          !write(fgw,'(" rrint:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! bradketc
          msize = sizeof(bradketc)*b2mb
          !write(fgw,'(" bradketc:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! bradketa
          msize = sizeof(bradketa)*b2mb
          !write(fgw,'(" bradketa:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! bradketlo
          msize = sizeof(bradketlo)*b2mb
          !write(fgw,'(" bradketlo:",T40,f8.2)') msize
          msize_tot = msize_tot+msize 
          !---------------------------------------
          msize = sizeof(gauntcoef)*b2mb
          !write(fgw,'(" gauntcoef:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          write(fgw,'(" Global data:",T40,f8.2)') msize
          !---------------------------------------
          ! vcxnn
          msize = nstsv*kset%nkpt*b2mb*16
          write(fgw,'(" <n|Vxc|n>:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! kintw = kiw+kwfer+ciw
          msize = (nstsv*kqset%nkpt+nstsv*kqset%nkpt+ncmax*natmtot)*b2mb*8
          write(fgw,'(" BZ integration weights:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! mpwipw
          msize = maxval(Gqbarc%ngk(1,:))*maxval(Gqset%ngk(1,:))*b2mb*16
          write(fgw,'(" PW-IPW overlap matrix:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! barc = barcev+vmat
          msize = (matsizmax*matsizmax)*b2mb*16+matsizmax*b2mb*8
          write(fgw,'(" Coulomb potential:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          msize_tot = msize_tot+msize
          ! KS eigenvectors = A*C + C + c.c.
          msize = 2*(nstsv*apwordmax*lmmaxapw*natmtot*nspnfv + &
          &          nmatmax*nstsv*nspnfv)*b2mb*16
          write(fgw,'(" KS eigenvectors:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! M^i_{nm}(k,q)
          msize = matsizmax*(nbandsgw+ncg)*nstsv*b2mb*16
          write(fgw,'(" Minm:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          ! G0W0 specific case
          if (trim(input%gw%taskname).ne.'g0w0_x') then
            ! pmatvv+pmatcv
            msize = (3*nstsv*nstsv + 3*ncg*nstsv)*b2mb*16
            write(fgw,'(" Momentum matrix elements:",T40,f8.2)') msize
            msize_tot = msize_tot+msize
            ! qdepwtet = kcw+unw
            msize = (nstsv*nstsv*freq%nomeg*kqset%nkpt + &
            &        natmtot*ncmax*nstsv*freq%nomeg*kqset%nkpt)*b2mb*16
            write(fgw,'(" (q,omega)-BZ integration weights:",T40,f8.2)') msize
            ! epsilon+wings
            msize = (matsizmax*matsizmax*freq%nomeg + &
            &        matsizmax*freq%nomeg*3)*b2mb*16
            write(fgw,'(" Dielectric function:",T40,f8.2)') msize
            msize_tot = msize_tot+msize
            ! M*W*M
            msize = nbandsgw*nstsv*freq%nomeg*b2mb*16
            write(fgw,'(" M*W*M:",T40,f8.2)') msize
            msize_tot = msize_tot+msize
          end if
          ! Self-energy = evalks + evalqp + \Sigma_x + \Sigma_c
          msize = 3*nbandsgw*kqset%nkpt*b2mb*16
          if (trim(input%gw%taskname).ne.'g0w0_x') then
            msize = msize+nbandsgw*kqset%nkpt*freq%nomeg*b2mb*16
          end if
          write(fgw,'(" Self-energy:",T40,f8.2)') msize
          msize_tot = msize_tot+msize
          write(fgw,*) '_______________________________________________'
          write(fgw,'(" Total G0W0:",T40,f8.2)') msize_tot
          call linmsg(fgw,'-','')
          call flushifc(fgw)
        end if
        
        return
    end subroutine
    
end subroutine init_gw
