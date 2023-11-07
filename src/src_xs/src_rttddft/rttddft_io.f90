module rttddft_io
  use asserts, only: assert
  use m_getunit, only: getunit
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
  use rttddft_GlobalVariables, only: Timing_RTTDDFT_and_MD, calculateNexc, &
    calculateTotalEnergy, predictorCorrector, printTimesDetailed
  
  implicit none

  private
  ! procedures
  public :: open_files_jpa, close_files_jpa, write_jpa, &
            open_file_etot, close_file_etot, write_total_energy, &
            open_file_nexc, close_file_nexc, write_nexc, &
            open_file_info, close_file_info, write_file_info, &
            write_file_info_header, write_file_info_fill_line_with_char, &
            open_file_timing, close_file_timing, write_timing, &
            file_pmat_exists, read_pmat, write_pmat, &
            file_pmat_mt_exists, read_pmat_mt, write_pmat_mt, &
            write_wavefunction

  !> Number of the unit to print timings
  integer                   :: file_time
  !> number of the unit to write the vector potential
  integer                   :: file_avec
  !> number of the unit to write the polarization field
  integer                   :: file_pvec
  !> number of the unit to write the current density
  integer                   :: file_jind
  !> number of the unit to write the number of excited electrons (per unit cell) 
  integer                   :: file_nexc
  !> number of the unit to write total energy
  integer                   :: file_etot
  !> number of the unit to write general information about the RT-TDDFT calculation
  integer                   :: file_info
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

  interface write_jpa
    module procedure :: write_jpa_single_line
    module procedure :: write_jpa_multiple_lines
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
    call getunit( unit )
    open( unit, file=trim(file_name), status='replace' )
  end subroutine

  !> Prints the current density \(\mathbf{J}\), or the polarization 
  !> \(\mathbf{P}\), or the vector potential \(\mathbf{A}\)
  subroutine write_jpa_multiple_lines( times, first, second, label )
    !> Array with the values of time \( t \)
    real(dp), intent(in) :: times(:)
    !> Array with the \( x, y, z \) components of \(\mathbf{J}\) , 
    !> \(\mathbf{P}\) or \(\mathbf{A}\) for each time \( t \)
    real(dp), intent(in) :: first(:, :)
    !> Same as before, but for the second array - usually \(\mathbf{A}\)
    real(dp), intent(in), optional :: second(:, :)
    !> String used to select the unit, where the data is printed out
    character(len=*), intent(in) :: label

    integer :: i, n
    logical :: twoArrays

    twoArrays = present( second )
    n = size( times )
    call assert( size(first, 2) == n, 'first array must have size = n along dim = 2')
    if ( twoArrays ) then 
      call assert( size(second, 2) == n, 'second array must have size = n along dim = 2')
    end if
    if( twoArrays ) then
      do i = 1, n
        call write_jpa_single_line( times(i), first(:,i), second(:,i), label )
      end do
    else 
      do i = 1, n
        call write_jpa_single_line( times(i), first(:,i), label=label )
      end do
    end if
  end subroutine 

  !> Prints the current density \(\mathbf{J}\), or the polarization 
  !> \(\mathbf{P}\), or the vector potential \(\mathbf{A}\)
  subroutine write_jpa_single_line( time, first, second, label )
    !> time \( t \)
    real(dp), intent(in) :: time
    !> Array with the \( x, y, z \) components of \(\mathbf{J}\) , 
    !> \(\mathbf{P}\) or \(\mathbf{A}\) for each time \( t \)
    real(dp), intent(in) :: first(:)
    !> Same as before, but for the second array - usually \(\mathbf{A}\)
    real(dp), intent(in), optional :: second(:)
    !> String used to select the unit, where the data is printed out
    character(len=*), intent(in) :: label

    integer :: unit

    call assert( trim(label)=='avec' .or. trim(label)=='pvec' .or. trim(label)=='jind', &
      'label must be one of the following strings: avec, pvec or jind' )
    call assert( size(first) == 3, 'first array must have size = 3')

    select case( trim(label) )
      case('avec') 
        unit = file_avec
        call assert( present( second ), '2 arrays must be passed for the case avec')
        call assert( size( second ) == 3, 'second array must have size = 3')
      case('pvec') 
        unit = file_pvec
      case('jind') 
        unit = file_jind
    end select

    if( present( second ) ) then
      write( unit, '(F9.3,6F20.12)' ) time, first(1), second(1), &
        & first(2), second(2), first(3), second(3)
    else
      write( unit, '(F9.3,3F20.12)' ) time, first(1), first(2), first(3)
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
    call getunit( file_time )
    open( file_time, file=trim(add_default_extension(filename_timing)), status='replace' )
  end subroutine

  subroutine close_file_timing
    close( file_time )
  end subroutine

  subroutine write_timing_initialization( timing_init )
    !> Time taken to initialize RT-TDDFT
    real(dp), intent(in) :: timing_init

    if( rank == 0 ) write( file_time, format_timing ) 'Initialization (sec):',timing_init
  end subroutine

  !> Subroutine to output the timings into `TIMING_RTTDDFT.OUT`
  subroutine write_timing_RTTDDFT_steps( itNumber, n, timing, screenshot_was_taken, molecular_dynamics )
    !> itNumber: The actual number of the counter that tells how many time steps 
    !> have already been executed
    integer, intent(in)             :: itNumber
    !> n:  Number of elements of the timing array
    integer, intent(in)             :: n
    !> timing: Array of timings. Each elements contains information
    !>   about how many seconds (timings) were spent in different parts of code
    type(Timing_RTTDDFT_and_MD), intent(in) :: timing(n)
    !> if `.True.`, a screenshot was taken at `itNumber`
    logical, intent(in)             :: screenshot_was_taken(n)
    !> Does timings about MD need to be printed?
    logical, intent(in), optional   :: molecular_dynamics
    

    integer                         :: ip, shift
    logical                         :: MD

    shift = itNumber - n
    MD = .False.
    if( present(molecular_dynamics) ) MD = molecular_dynamics

    if ( rank == 0 ) then
      do ip = 1, n
        associate( t_rttddft => timing(ip)%t_RTTDDFT )
        write(file_time,'(A30,I10)')'Time (sec) spent in iteration:',ip+shift
        write(file_time,format_timing) 'updatewvf:',t_rttddft%t_wvf
        write(file_time,format_timing) 'updatedens:',t_rttddft%t_dens
        if ( printTimesDetailed ) then
          write(file_time,format_timing) '-- rhovalk and genrhoir:',t_rttddft%t_dens_rho
          write(file_time,format_timing) '-- symrf:',t_rttddft%t_dens_symrf
          write(file_time,format_timing) '-- rfmtctof:',t_rttddft%t_dens_rfmtctof
          write(file_time,format_timing) '-- addrhocr:',t_rttddft%t_dens_addrhocr
          write(file_time,format_timing) '-- charge:',t_rttddft%t_dens_charge
          write(file_time,format_timing) '-- rhonorm:',t_rttddft%t_dens_rhonorm
        end if
        write(file_time,format_timing) 'updatepot:',t_rttddft%t_uppot
        if ( printTimesDetailed ) then
          write(file_time,format_timing) '-- poteff:',t_rttddft%t_poteff
          write(file_time,format_timing) '-- genveffig:',t_rttddft%t_genveffig
          write(file_time,format_timing) '-- genmeffig:',t_rttddft%t_genmeffig
        end if
        write(file_time,format_timing) 'UpdateCurrentDensity:',t_rttddft%t_curr
        write(file_time,format_timing) 'ObtainA:',t_rttddft%t_obtaina
        write(file_time,format_timing) 'updatehamiltonian:',t_rttddft%t_upham
        if ( printTimesDetailed ) then
          write(file_time,format_timing) '-- hmlint:',t_rttddft%t_hmlint
          write(file_time,format_timing) '-- other subs:',t_rttddft%t_ham
        end if
        if ( predictorCorrector )  &
          & write(file_time,format_timing) 'All cycles of predcorr:',t_rttddft%t_predcorr
        if ( calculateTotalEnergy .and. printTimesDetailed ) write(file_time,format_timing) 'Total Energy:',t_rttddft%t_toten
        if ( calculateNexc .and. printTimesDetailed ) write(file_time,format_timing)'nexc:',t_rttddft%t_nexc
        if ( screenshot_was_taken(ip) ) write(file_time,format_timing) 'Screenshots:',t_rttddft%t_screenshot
        end associate
        if( MD ) then
          associate( t_MD => timing(ip)%t_Ehrenfest )
          if( t_MD%MD_was_carried_out ) then
            write(file_time,format_timing) 'MD:', t_MD%t_MD_step
            if( printTimesDetailed ) then
              write(file_time,format_timing) '-- 1st part of forces:', t_MD%t_MD_1st
              write(file_time,format_timing) '-- 2nd part of forces:', t_MD%t_MD_2nd
              write(file_time,format_timing) '-- sum forces:', t_MD%t_MD_sumforces 
              write(file_time,format_timing) '-- move ions:', t_MD%t_MD_moveions
              write(file_time,format_timing) '-- update basis:', t_MD%t_MD_updateBasis
              write(file_time,format_timing) '-- update H, S:', t_MD%t_MD_hamoverl
              write(file_time,format_timing) '-- update pmat:', t_MD%t_MD_pmat
            end if
          end if
          end associate
        end if
        write(file_time,format_timing) 'time per iteration:',timing(ip)%t_iteration
      end do
    end if
  end subroutine

  logical function file_pmat_exists()
    inquire( file=trim(add_default_extension(filename_pmat)), exist=file_pmat_exists )
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

end module