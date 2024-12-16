module rttddft_io
  use asserts, only: assert
  use mod_misc, only: filext, versionname, githash
  use modinput, only: input
  use modmpi, only: rank, procs, barrier
  use mod_mpi_env, only: mpiinfo
#ifdef MPI
  use rttddft_io_parallel, only: read_array, write_array
#else
  use rttddft_io_serial, only: read_array, write_array
#endif
  use precision, only: dp, i32
  use rttddft_Energy, only: TotalEnergy
  use rttddft_CurrentDensity, only: Current_Density_Field
  use rttddft_Polarization, only: Polarization
  use rttddft_timings, only: Print_Timings, Timing_RTTDDFT_and_MD
  use rttddft_VectorField, only: Uniform_Vector_Field, x, y, z
  use rttddft_VectorPotential, only: Vector_Potential_Field
  
  implicit none

  private
  ! procedures
  public :: open_files_jpa, close_files_jpa, write_jpa, &
            open_file_etot, close_file_etot, write_total_energy, &
            open_file_nexc, close_file_nexc, write_nexc, &
            open_file_info, close_file_info, write_file_info, &
            write_file_info_header, write_file_info_fill_line_with_char, &
            open_file_timing, close_file_timing, write_timing, &
            file_pmat_exists, read_pmat, write_pmat, get_filename_pmat, &
            file_pmat_mt_exists, read_pmat_mt, write_pmat_mt, get_filename_pmat_mt, &
            write_wavefunction, write_real_function_xsf, transform_real_function_to_rgrid

  !> Number of the unit to print timings
  integer(i32)                   :: file_time
  !> number of the unit to write the vector potential
  integer(i32)                   :: file_avec
  !> number of the unit to write the polarization field
  integer(i32)                   :: file_pvec
  !> number of the unit to write the current density
  integer(i32)                   :: file_jind
  !> number of the unit to write the number of excited electrons (per unit cell) 
  integer(i32)                   :: file_nexc
  !> number of the unit to write total energy
  integer(i32)                   :: file_etot
  !> number of the unit to write general information about the RT-TDDFT calculation
  integer(i32)                   :: file_info
  !> Format of the timing outputs in RT-TDDFT
  character(len=*), parameter :: format_timing = '(A30,F12.6)'
  !> Format of the outputs: `JIND` and `PVEC`
  character(len=*), parameter :: format_j_p = '(F12.4,3F20.12)'
  !> Format of the output: `AVEC`
  character(len=*), parameter :: format_avec = '(F12.4,6F20.12)'
  !> Default name of the file where the vector potential is printed out
  character(len=*), parameter :: filename_avec = 'AVEC'
  !> Default name of the file where the polarization is printed out
  character(len=*), parameter :: filename_pvec = 'PVEC'
  !> Default name of the file where the current density is printed out
  character(len=*), parameter :: filename_jind = 'JIND'
  !> Default name of the file where the total energy is printed out
  character(len=*), parameter :: filename_etot = 'ETOT_RTTDDFT'
  !> Default name of the file where the number of excited electrons is printed out
  character(len=*), parameter :: filename_nexc = 'NEXC'
  !> Default name of the file with general information about the RT-TDDFT calculation
  character(len=*), parameter :: filename_info = 'RTTDDFT_INFO'
  !> Default name of the file where timigs are printed out
  character(len=*), parameter :: filename_timing = 'TIMING_RTTDDFT'
  !> Default name of the file where `pmat` is printed out
  character(len=*), parameter :: filename_pmat = 'PMATBASIS'
  !> Default name of the file where `pmat_mt` is printed out
  character(len=*), parameter :: filename_pmat_mt = 'PMATMTBASIS'
  

  interface write_timing
    module procedure :: write_timing_initialization
    module procedure :: write_timing_RTTDDFT_steps
  end interface

contains 
  !> (private) add the default extension (usually .OUT) to the base file name
  pure function add_default_extension( file_name )
    !> base file name
    character(len=*), intent(in)  :: file_name
    !> file name with default extension
    character(len=:), allocatable :: add_default_extension
    add_default_extension = trim(file_name)//filext
  end function

  !> (private) generic subroutine to open a file
  subroutine open_file_generic( unit, file_name )
    !> unit number of file to open
    integer, intent(out) :: unit
    !> name of file to open
    character(len=*), intent(in) :: file_name
    
    open( newunit=unit, file=trim(file_name), status='replace' )
  end subroutine

  !> Prints the current density \(\mathbf{J}\), or the polarization 
  !> \(\mathbf{P}\), or the vector potential \(\mathbf{A}\)
  subroutine write_jpa( times, first, second )
    !> Array with the values of time \( t \)
    real(dp), intent(in) :: times(:)
    !> Array with the \( x, y, z \) components of \(\mathbf{J}\) , 
    !> \(\mathbf{P}\) or \(\mathbf{A}\) for each time \( t \)
    class(Uniform_Vector_Field), intent(in) :: first(:)
    !> Same as before, but for the second array - usually \(\mathbf{A}\)
    class(Uniform_Vector_Field), optional :: second(:)

    integer(i32) :: i, n, unit
    logical :: twoArrays

    twoArrays = present( second )
    n = size( times )

    call assert( size(first) == n, 'first array must have size = n')
    if( twoArrays ) call assert( size(second) == n, 'second array must have size = n')
    
    select type( first )
      type is( Vector_Potential_Field )
        unit = file_avec
        call assert( twoArrays, '2nd argument must be passed for the case of Vector_Field')
        call assert( same_type_as( first, second ), '2nd argument must be of type(Vector_Field)')
      type is( Polarization )
        unit = file_pvec
      type is( Current_Density_Field )
        unit = file_jind
      class default
        call assert( .false., 'unrecognized type passed to write_jpa')
    end select

    if( twoArrays ) then
      do i = 1, n
        write( unit, '(F9.3,6F20.12)' ) times(i), first(i)%components(x), second(i)%components(x), &
          & first(i)%components(y), second(i)%components(y), first(i)%components(z), second(i)%components(z)
      end do
    else 
      do i = 1, n
        write( unit, '(F9.3,3F20.12)' ) times(i), first(i)%components(x), first(i)%components(y), first(i)%components(z)
      end do
    end if
  end subroutine 

  !> Open files for writing jind, pvec and avec
  subroutine open_files_jpa
    call open_file_generic( file_jind, add_default_extension(filename_jind) )
    call open_file_generic( file_pvec, add_default_extension(filename_pvec) )
    call open_file_generic( file_avec, add_default_extension(filename_avec) )
  end subroutine

  subroutine close_files_jpa
    close( file_jind )
    close( file_pvec )
    close( file_avec )
  end subroutine

  subroutine open_file_etot
    call open_file_generic( file_etot, add_default_extension(filename_etot) )
  end subroutine

  subroutine close_file_etot
    close( file_etot )
  end subroutine

  !> Print out the total energy
  subroutine write_total_energy( printHeader, nArrayElements, &
    & timeArray, etotArray )
    !> Number of lines to printed = number of elements of the arrays:
    !> `timeArray` and `etotArray`.
    integer, intent(in) :: nArrayElements
    !> Tells if a header must be printed (useful when the file is opened for 
    !> the 1st time)
    logical, intent(in) :: printHeader
    !> Array with the values of time \( t \)
    real(8), intent(in) :: timeArray(nArrayElements)
    !> Array with the energies (total energy, XC, Madelung, etc.)
    type(TotalEnergy), intent(in) :: etotArray(nArrayElements)

    integer :: i

    if ( printHeader ) then
      write(file_etot,'(A9,8A20)') 'Time','ETOT','Madelung','Eigenvalues-Core',&
        & 'Eigenvalues-Valence','Exchange','Correlation','XC-potential',&
        & 'Coulomb pot. energy'
    end if
    do i = 1, nArrayElements
      write(file_etot,'(F9.3,8F20.10)') timeArray(i), &
        & etotArray(i)%total_energy, etotArray(i)%madelung, &
        & etotArray(i)%eigenvalues_core, etotArray(i)%hamiltonian,&
        & etotArray(i)%exchange, etotArray(i)%correlation, &
        & etotArray(i)%integral_vxc_times_density, etotArray(i)%Coulomb
    end do
  end subroutine

  subroutine open_file_nexc
    call open_file_generic( file_nexc, add_default_extension(filename_nexc) )
  end subroutine

  subroutine close_file_nexc
    close( file_nexc )
  end subroutine

  !> Prints the number of excitations
  subroutine write_nexc( printHeader, nArrayElements, &
    & timeArray, nex, ngs, ntot )
    !> Number of lines to printed = number of elements of the arrays:
    !> `nex`, `ngs` and `ntot`.
    integer, intent(in)   :: nArrayElements
    !> printHeader: If we need to print a header (useful when we open the file for the 1st time)
    logical, intent(in)   :: printHeader
    !> timeArray         array with the values of time \( t \)
    real(dp), intent(in)  :: timeArray(nArrayElements)
    !> nex               Number of electrons which were excited
    real(dp), intent(in)  :: nex(nArrayElements)
    !> ngs               Number of electrons on the groundstate
    real(dp), intent(in)  :: ngs(nArrayElements)
    !> ntot              Sum of ngs and nex
    real(dp), intent(in)  :: ntot(nArrayElements)

    integer :: i

    if ( printHeader ) then
      write( file_nexc,'(A9,3A20)' ) 'Time','N.Elec.GS', &
        & 'N.XS', 'Sum'
    end if
    do i = 1, nArrayElements
      write( file_nexc, '(F9.3,3F20.10)' ) timeArray(i), &
        & ngs(i), nex(i), ntot(i)
    end do
  end subroutine

  subroutine open_file_info
    call open_file_generic( file_info, add_default_extension(filename_info) )
  end subroutine

  subroutine close_file_info
    close( file_info )
  end subroutine

  subroutine write_file_info( string, string_format )
    !> string to be printed out
    character(len=*), intent(in) :: string
    !> format for printing out `string`
    character(len=*), intent(in), optional :: string_format

    character(len=:), allocatable :: format_
    if( present(string_format) ) then
      format_ = string_format
    else
      format_='(A)'
    end if
    write( file_info, format_ ) string
  end subroutine

  subroutine write_file_info_fill_line_with_char( ch )
    character, intent(in) :: ch
    call printline(file_info, ch)
  end subroutine

  subroutine write_file_info_header
    character(len=100)      :: string
    character(10)           :: dat, tim

    call write_file_info( 'Real-time TDDFT calculation started' )
    call write_file_info( 'EXCITING '//trim(versionname)//' started' ) 
    if (len(trim(githash)) > 0) call write_file_info('version hash id: '//trim(githash))
#ifdef MPI
    write( string, '(A,I6,A)') 'MPI version using ', procs, ' processor(s)'
    call write_file_info( string )
#endif
    call date_and_time(date=dat, time=tim)
    write( string, '("Date (DD-MM-YYYY) : ", A2, "-", A2, "-", A4)') &
    &  dat (7:8), dat (5:6), dat (1:4)
    call write_file_info( string )
    write( string, '("Time (hh:mm:ss)   : ", A2, ":", A2, ":", A2)') &
    &  tim (1:2), tim (3:4), tim (5:6)
    call write_file_info( string )
    call write_file_info( 'All units are atomic (Hartree, Bohr, etc.)' )
  end subroutine

  subroutine open_file_timing
    call open_file_generic( file_time, add_default_extension(filename_timing) )
  end subroutine

  subroutine close_file_timing
    close( file_time )
  end subroutine

  subroutine write_timing_initialization( timing_init )
    !> Time taken to initialize RT-TDDFT
    real(dp), intent(in) :: timing_init

    call write_nonzero_timing( 'Initialization (sec):', timing_init )
  end subroutine

  !> Subroutine to output the timings into `TIMING_RTTDDFT.OUT`
  subroutine write_timing_RTTDDFT_steps( itNumber, timing, molecular_dynamics )
    !> itNumber: The actual number of the counter that tells how many time steps 
    !> have already been executed
    integer, intent(in) :: itNumber
    !> timing: Array of timings. Each elements contains information
    !>   about how many seconds (timings) were spent in different parts of code
    type(Timing_RTTDDFT_and_MD), intent(in) :: timing(:)
    !> Does timings about MD need to be printed?
    logical, intent(in), optional :: molecular_dynamics
    
    integer  :: ip, shift, n
    logical  :: MD

    n = size( timing )

    shift = itNumber - n
    MD = .False.
    if( present(molecular_dynamics) ) MD = molecular_dynamics

    do ip = 1, n
      associate( t_rttddft => timing(ip)%t_RTTDDFT )
      write( file_time, '(A30,I10)' ) 'Time (sec) spent in iteration:', ip + shift
      call write_nonzero_timing( 'updatewvf:', t_rttddft%wavefunction )
      call write_nonzero_timing( 'updatedens:', t_rttddft%dens%total )
      call write_nonzero_timing( '-- rhovalk and genrhoir:', t_rttddft%dens%rho )
      call write_nonzero_timing( '-- symrf:', t_rttddft%dens%symrf )
      call write_nonzero_timing( '-- rfmtctof:', t_rttddft%dens%rfmtctof )
      call write_nonzero_timing( '-- addrhocr:', t_rttddft%dens%addrhocr )
      call write_nonzero_timing( '-- charge:', t_rttddft%dens%charge )
      call write_nonzero_timing( '-- rhonorm:', t_rttddft%dens%rhonorm )
      call write_nonzero_timing( 'updatepot:', t_rttddft%pot%total )
      call write_nonzero_timing( '-- poteff:', t_rttddft%pot%poteff )
      call write_nonzero_timing( '-- genveffig:', t_rttddft%pot%genveffig )
      call write_nonzero_timing( '-- genmeffig:', t_rttddft%pot%genmeffig )
      call write_nonzero_timing( 'UpdateCurrentDensity:', t_rttddft%current_density )
      call write_nonzero_timing( 'ObtainA:', t_rttddft%vector_potential )
      call write_nonzero_timing( 'updatehamiltonian:', t_rttddft%ham%total )
      call write_nonzero_timing( '-- hmlint:', t_rttddft%ham%hmlint )
      call write_nonzero_timing( '-- other subs:', t_rttddft%ham%rest )
      call write_nonzero_timing( 'All cycles of predcorr:', t_rttddft%pred_corr )
      call write_nonzero_timing( 'Total Energy:', t_rttddft%energy )
      call write_nonzero_timing( 'nexc:', t_rttddft%n_exc )
      call write_nonzero_timing( 'Screenshots:', t_rttddft%screenshot )
      end associate
      if( MD ) then
        associate( t_MD => timing(ip)%t_Ehrenfest )
        call write_nonzero_timing( 'MD:', t_MD%t_MD_step )
        call write_nonzero_timing( '-- 1st part of forces:', t_MD%t_MD_1st )
        call write_nonzero_timing( '-- 2nd part of forces:', t_MD%t_MD_2nd )
        call write_nonzero_timing( '-- sum forces:', t_MD%t_MD_sumforces )
        call write_nonzero_timing( '-- move ions:', t_MD%t_MD_moveions )
        call write_nonzero_timing( '-- update basis:', t_MD%t_MD_updateBasis )
        call write_nonzero_timing( '-- update H, S:', t_MD%hamoverl )
        call write_nonzero_timing( '-- update pmat:', t_MD%pmat )
        end associate
      end if
      write( file_time, format_timing ) 'time per iteration:', timing(ip)%t_iteration
    end do
  end subroutine

  logical function file_pmat_exists()
    inquire( file=trim(add_default_extension(filename_pmat)), exist=file_pmat_exists )
  end function

  function get_filename_pmat() result(name)
    character(len=:), allocatable :: name
    name = add_default_extension( filename_pmat )
  end function

  function get_filename_pmat_mt() result(name)
    character(len=:), allocatable :: name
    name = add_default_extension( filename_pmat_mt )
  end function

  !> Read the momentum matrix elements from file
  subroutine read_pmat( first_kpt, pmat, mpi_env )
    !> Index of the first `k-point` to be considered
    integer,intent(in)        :: first_kpt
    !> Momentum matrix elements
    complex(dp), intent(out)  :: pmat(:, :, :, first_kpt:)
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(inout) :: mpi_env
    
    call read_array(add_default_extension( filename_pmat ), first_kpt, pmat, mpi_env )
  end subroutine

  !> Write the momentum matrix elements to file
  subroutine write_pmat( first_kpt, pmat, mpi_env )
    !> Index of the first `k-point` to be considered in the sum
    integer,intent(in)        :: first_kpt
    !> Momentum matrix elements
    complex(dp), intent(in)   :: pmat(:, :, :, first_kpt:)
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(inout) :: mpi_env
    
    call write_array(add_default_extension( filename_pmat ), first_kpt, pmat, mpi_env )
  end subroutine

  !> Check if file with `pmat_mt` exists
  logical function file_pmat_mt_exists()
    inquire( file=trim(add_default_extension(filename_pmat_mt)), exist=file_pmat_mt_exists )
  end function

  !> Read the muffin-tin part of the momentum matrix (`pmat_mt`) from file
  subroutine read_pmat_mt( first_kpt, pmat_mt, mpi_env )
    !> Index of the first `k-point` to be considered in the sum
    integer,intent(in)        :: first_kpt
    !> Muffin-tin part of the momentum matrix
    complex(dp), intent(out)  :: pmat_mt(:, :, :, :, first_kpt:)
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(inout) :: mpi_env

    call read_array( add_default_extension( filename_pmat_mt ), first_kpt, pmat_mt, mpi_env )
  end subroutine

  !> Write the muffin-tin part of the momentum matrix (`pmat_mt`) to file
  subroutine write_pmat_mt( first_kpt, pmat_mt, mpi_env )
    !> Index of the first `k-point` to be considered in the sum
    integer,intent(in)        :: first_kpt
    !> Muffin-tin part of the momentum matrix
    complex(dp), intent(in)   :: pmat_mt(:, :, :, :, first_kpt:)
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(inout) :: mpi_env
    
    call write_array( add_default_extension( filename_pmat_mt ), first_kpt, pmat_mt, mpi_env )
  end subroutine

  subroutine write_wavefunction( first_kpt, wavefunction )
    !> index of the first `k-point` to be considered in the sum
    integer,intent(in)        :: first_kpt
    !> wavefunction coefficients
    complex(dp), intent(in)   :: wavefunction(:, :, first_kpt:)
    
    integer(i32) :: count, ik, last_kpt

    last_kpt = ubound( wavefunction, 3 )
    do count = 1, procs
      if ( rank == count-1 ) then
        do ik = first_kpt, last_kpt
          call putevecfv( ik, wavefunction(:,:,ik) )
        end do
      end if
      call barrier()
    end do
  end subroutine

    !> Write timing only if it is nonzero (> tol)
  subroutine write_nonzero_timing( description, timing )
    !> Action name
    character(len = *), intent(in) :: description
    !> Action duration
    real(dp), intent(in) :: timing

    real(dp), parameter :: tol = 1.0e-7_dp

    if ( timing > tol ) write( file_time, format_timing ) description, timing

  end subroutine

  !> Write real-valued coordinate-space-defined 3d function in xsf file
  subroutine write_real_function_xsf( grid, iteration, function_rgrid, label )
    use mod_xsf_format, only: write_structure_xsf, write_3d_xsf
    use mod_rgrid, only: rgrid

    implicit none
    !> Pre-generated grid
    type(rgrid), intent(in) :: grid
    !> Iteration number used in the filename
    integer, intent(in) :: iteration
    !> Real-valued function on the grid (grid%npt)
    real(dp), intent(in) :: function_rgrid(:)
    !> User-defined function label (e.g. observable name)
    character(len = *), intent(in) :: label
  
    character(80) :: fname
    
    write( fname, '("-",i5,".xsf")' ) iteration
    fname = trim( label )//fname
    call str_strip( fname )
    call write_structure_xsf( fname )
    call write_3d_xsf( fname, label, grid%boxl(1 : 4, :), grid%ngrid, &
    grid%npt, function_rgrid )

  end subroutine

  !> Convert real-valued coordinate-space-defined 3d function from the 
  !> IR-MT representation to the coordinate representation
  subroutine transform_real_function_to_rgrid( grid, lmax, function_mt, function_ir, function_rgrid )
    use mod_rgrid, only: rgrid

    implicit none
    !> Pre-generated grid
    type(rgrid), intent(in) :: grid
    !> Maximum value of l used for the MT expansions
    integer, intent(in) :: lmax
    !> Real-valued function in MT (lmmaxvr, nrmtmax, natmtot)
    real(dp), intent(in) :: function_mt(:, :, :)
    !> Real-valued function in IR region (ngrtot)
    real(dp), intent(in) :: function_ir(:)
    !> Real-valued function on the grid (grid%npt)
    real(dp), intent(out) :: function_rgrid(:)

    call rfarray( lmax, size( function_mt, 1 ), function_mt, &
    function_ir, grid%npt, grid%vpl, function_rgrid )

  end subroutine 

end module