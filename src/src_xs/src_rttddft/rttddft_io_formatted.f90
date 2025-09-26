module rttddft_io_formatted
  use asserts, only: assert
  use file_utils, only: add_default_extension, copy_text_file, delete_file, read_last_and_penultimate_lines_from_file
  use mod_misc, only: filext, githash, versionname
  use mod_rgrid, only: rgrid, gen_3d_rgrid
  use mod_xsf_format, only: add_xsf_extension, write_real_function_xsf
  use modinput, only: plot3d_type
  use modmpi, only: procs
  use precision, only: dp, i32, str_16, str_128, str_256
  use rttddft_CurrentDensity, only: Current_Density_Field
  use rttddft_electric_field, only: Electric_Field
  use rttddft_file_names
  use rttddft_Energy, only: Total_Energy
  use rttddft_Polarization, only: Polarization
  use rttddft_timings, only: Print_Timings, Timing_RTTDDFT_and_MD
  use rttddft_VectorField, only: Uniform_Vector_Field, x, y, z
  use rttddft_VectorPotential, only: Vector_Potential_Field
  use to_char_conversion, only: to_char

  implicit none

  private

  public :: close_file_etot, &
            close_file_info, &
            close_file_nexc, &
            close_file_timing, &
            close_files_vector_fields, &
            copy_files, &
            delete_jpa_files, &
            open_file_etot, &
            open_file_info, &
            open_file_nexc, &
            open_file_timing, &
            open_files_vector_fields, &
            read_vector_field, &
            write_density_to_file, &
            write_eigenvalues, &
            write_file_info_header, &
            write_file_info_fill_line_with_char, &
            write_file_info, &
            write_nexc, &
            write_occupations, &
            write_projection_coefficients, &
            write_timing, &
            write_total_energy, &
            write_vector_field

  !> Number of the unit to print timings
  integer(i32) :: file_time
  !> number of the unit to write the vector potential
  integer(i32) :: file_avec
  !> number of the unit to write the electric field
  integer(i32) :: file_evec
  !> number of the unit to write the polarization field
  integer(i32) :: file_pvec
  !> number of the unit to write the current density
  integer(i32) :: file_jind
  !> number of the unit to write the number of excited electrons (per unit cell) 
  integer(i32) :: file_nexc
  !> number of the unit to write total energy
  integer(i32) :: file_etot
  !> number of the unit to write general information about the RT-TDDFT calculation
  integer(i32) :: file_info
  !> Format of the timing outputs in RT-TDDFT
  character(len=*), parameter :: format_timing = '(A30,F12.6)'
  !> Format of the outputs: `CURRENT` and `POLARIZATION`
  character(len=*), parameter :: format_j_p = '(F12.4,3F20.12)'
  !> Format of the output: `VECTOR_POTENTIAL`
  character(len=*), parameter :: format_avec = '(F12.4,6F20.12)'

  interface write_timing
    module procedure :: write_timing_initialization
    module procedure :: write_timing_RTTDDFT_steps
  end interface

contains
  !> (private) Function to return the status of a file to open given the information if it is new or old
  pure function get_status_from_logical( new ) result(status)
    !> If `.true.`, a new file is created (overwriting an exisiting one, if this is the case)
    logical, intent(in) :: new
    character(len=:), allocatable :: status
    character(len=*), parameter :: status_new = "replace"
    character(len=*), parameter :: status_old = "old    "

    ! Cray 18 fails to concatenate trim to merge
    ! While status = trim( merge( status_new, status_old, new ) ) is valid
    ! we separate them
    status = merge( status_new, status_old, new )
    status = trim( status)
  end function

  !> (private) Function to return the position of a file to open given the information if it is new or old
  pure function get_position_from_logical( new ) result(position)
    !> If `.true.`, a new file is created (overwriting an exisiting one, if this is the case)
    logical, intent(in) :: new
    character(len=:), allocatable :: position
    character(len=*), parameter :: position_new = "rewind"
    character(len=*), parameter :: position_old = "append"

    position = merge( position_new, position_old, new )
  end function

  !> Read the last line, and if required the penultimate line too, of files with 
  !> the current \(\mathbf{J}\), or the polarization \(\mathbf{P}\), the electric field 
  !> \(\mathbf{E}\), or the vector potential \(\mathbf{A}\)
  subroutine read_vector_field( time, field_t, field_t_minus_dt, a_tot_t, a_tot_t_minus_dt )
    !> Time \(t\) contained in the last line
    real(dp), intent(out) :: time
    !> \(\mathbf{J}\) , \(\mathbf{P}\) or \(\mathbf{A}_{ind}\) at time \( t \)
    class(Uniform_Vector_Field), intent(out) :: field_t
    !> \(\mathbf{J}\) , \(\mathbf{P}\) or \(\mathbf{A}_{ind}\) at time \( t - \Delta t\)
    class(Uniform_Vector_Field), optional, intent(out) :: field_t_minus_dt
    !> Usually \(\mathbf{A}_{tot}\) at time \( t \)
    class(Vector_Potential_Field), optional, intent(out) :: a_tot_t
    !> \(\mathbf{J}\) , \(\mathbf{P}\) or \(\mathbf{A}\) at time \( t - \Delta t\)
    class(Vector_Potential_Field), optional, intent(out) :: a_tot_t_minus_dt

    character(len=str_256) :: last_line, penultimate_line, file_name
    if( present( field_t_minus_dt ) ) then 
      call assert( same_type_as( field_t, field_t_minus_dt ), '2nd argument must be of type(Vector_Field)')
    end if
    select type( field_t )
      type is( Vector_Potential_Field )
        file_name = filename_avec
        call assert( present( a_tot_t ), "a_tot_t must be passed" )
        if( present( field_t_minus_dt ) ) then 
          call assert( present( a_tot_t_minus_dt ), "a_tot_t_minus_dt must be passed")
        end if
      type is( Polarization )
        file_name = filename_pvec
      type is( Current_Density_Field )
        file_name = filename_jind
      type is( Electric_Field )
        file_name = filename_evec
      class default
        call assert( .false., 'unrecognized type passed to read_vector_field' )
    end select

    call read_last_and_penultimate_lines_from_file( add_default_extension( file_name), last_line, penultimate_line )
    
    if( present(field_t_minus_dt) ) then
      associate( v0 => field_t_minus_dt%components )
        if( present(a_tot_t_minus_dt) ) then
          associate( a0 => a_tot_t_minus_dt%components )
            read( penultimate_line, * ) time, v0(x), a0(x), v0(y), a0(y), v0(z), a0(z)
          end associate
        else
          read( penultimate_line, * ) time, v0
        end if
      end associate
    end if
    associate( v => field_t%components )
      if( present(a_tot_t) ) then
        associate( a => a_tot_t%components )
          read( last_line, * ) time, v(x), a(x), v(y), a(y), v(z), a(z)
        end associate
      else
        read( last_line, * ) time, v
      end if
    end associate
  end subroutine

  !> Copy files. Sources are files with name `fname` appended with [[add_default_extension]]
  !> and with an `extra_extension`. `fname` can be [[filename_avec]], [[filename_pvec]],
  !> [[filename_jind]], [[filename_evec]], [[filename_nexc]] (if `nexc` is `.true.`), and 
  !> [[filename_etot]] (if `etot` is `.true.`)
  subroutine copy_files( extra_extension, nexc, etot )
    !> Extension of source files
    character(len=*), intent(in) :: extra_extension
    !> If `.true.`, copy [[filename_nexc]]
    logical, intent(in) :: nexc
    !> If `.true.`, copy [[filename_etot]]
    logical, intent(in) :: etot
    
    call wrapper_copy_file( filename_avec, extra_extension )
    call wrapper_copy_file( filename_evec, extra_extension )
    call wrapper_copy_file( filename_pvec, extra_extension )
    call wrapper_copy_file( filename_jind, extra_extension )
    if( nexc ) call wrapper_copy_file( filename_nexc, extra_extension )
    if( etot ) call wrapper_copy_file( filename_etot, extra_extension )
  contains
    subroutine wrapper_copy_file( file_name, src_extra_extension )
      character(len=*), intent(in) :: file_name
      character(len=*), intent(in) :: src_extra_extension

      call copy_text_file( source_name=add_default_extension(file_name)//trim(src_extra_extension), &
        destination_name=add_default_extension(file_name) )
    end subroutine
  end subroutine

  !> Prints the current density \(\mathbf{J}\), or the polarization 
  !> \(\mathbf{P}\), or the vector potential \(\mathbf{A}\), or the electric 
  !> field \(\mathbf{E}\)
  subroutine write_vector_field( times, first, second )
    !> Array with the values of time \( t \)
    real(dp), intent(in) :: times(:)
    !> Array with the \( x, y, z \) components of \(\mathbf{J}\) , 
    !> \(\mathbf{P}\), \(\mathbf{E}\) or \(\mathbf{A}\) for each time \( t \)
    class(Uniform_Vector_Field), intent(in) :: first(:)
    !> Same as before, but for the second array - usually \(\mathbf{A}\)
    class(Uniform_Vector_Field), optional :: second(:)

    integer(i32) :: i, n, unit
    logical :: twoArrays

    twoArrays = present( second )
    n = size( times )

    call assert( size( first ) == n, 'first array must have size = n')
    if( twoArrays ) call assert( size( second ) == n, 'second array must have size = n')
    
    select type( first )
      type is( Vector_Potential_Field )
        unit = file_avec
        call assert( twoArrays, '2nd argument must be passed for the case of Vector_Field')
        call assert( same_type_as( first, second ), '2nd argument must be of type(Vector_Field)')
      type is( Polarization )
        unit = file_pvec
      type is( Current_Density_Field )
        unit = file_jind
      type is( Electric_Field )
        unit = file_evec
      class default
        call assert( .false., 'unrecognized type passed to write_vector_field')
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

  !> Open files for writing jind, pvec, avec, and e_field
  subroutine open_files_vector_fields( new )
    !> If `.true.`, open new files (rewriting old ones).
    !> If `.false.`, append to existing files
    logical, intent(in) :: new

    open( newunit=file_jind, file=add_default_extension(filename_jind), status=get_status_from_logical( new ), position=get_position_from_logical( new ) )
    open( newunit=file_pvec, file=add_default_extension(filename_pvec), status=get_status_from_logical( new ), position=get_position_from_logical( new ) )
    open( newunit=file_avec, file=add_default_extension(filename_avec), status=get_status_from_logical( new ), position=get_position_from_logical( new ) )
    open( newunit=file_evec, file=add_default_extension(filename_evec), status=get_status_from_logical( new ), position=get_position_from_logical( new ) )
  end subroutine

  subroutine close_files_vector_fields()
    close( file_jind )
    close( file_pvec )
    close( file_avec )
    close( file_evec )
  end subroutine

  subroutine delete_jpa_files()
    integer(i32) :: i_error
    call delete_file( add_default_extension( filename_jind ), i_error )
    call delete_file( add_default_extension( filename_pvec ), i_error )
    call delete_file( add_default_extension( filename_avec ), i_error )
    call delete_file( add_default_extension( filename_evec ), i_error )
  end subroutine

  subroutine open_file_etot( new )
    !> If `.true.`, open new files (rewriting old ones).
    !> If `.false.`, append to existing files
    logical, intent(in) :: new
    open( newunit=file_etot, file=add_default_extension(filename_etot), &
      status=get_status_from_logical( new ), position=get_position_from_logical( new ) )
  end subroutine

  subroutine close_file_etot()
    close( file_etot )
  end subroutine

  !> Print out the total energy
  subroutine write_total_energy( print_header, time_array, e_tot_array )
    !> Tells if a header must be printed (useful when the file is opened for 
    !> the 1st time)
    logical, intent(in) :: print_header
    !> Array with the values of time \( t \)
    real(dp), intent(in) :: time_array(:)
    !> Array with the energies (total energy, XC, Madelung, etc.)
    type(Total_Energy), intent(in) :: e_tot_array(:)

    integer(i32) :: i

    if ( print_header ) then
      write(file_etot,'(A9,8A20)') 'Time','ETOT','Madelung','Eigenvalues-Core',&
        & 'Eigenvalues-Valence','Exchange','Correlation','XC-potential',&
        & 'Coulomb pot. energy'
    end if
    associate( n => size(time_array) )
      call assert( size(e_tot_array) == n, 'e_tot_array must contain n elements')
      do i = 1, n
        write(file_etot,'(F9.3,8F20.10)') time_array(i), &
          & e_tot_array(i)%total_energy(), e_tot_array(i)%madelung, &
          & e_tot_array(i)%eigenvalues_core, e_tot_array(i)%hamiltonian,&
          & e_tot_array(i)%exchange, e_tot_array(i)%correlation, &
          & e_tot_array(i)%integral_vxc_times_density, e_tot_array(i)%Coulomb
      end do
    end associate
  end subroutine

  subroutine open_file_nexc( new )
    !> If `.true.`, open new files (rewriting old ones).
    !> If `.false.`, append to existing files
    logical, intent(in) :: new
    open( newunit=file_nexc, file=add_default_extension(filename_nexc), status=get_status_from_logical( new ), position=get_position_from_logical( new ) )
  end subroutine

  subroutine close_file_nexc()
    close( file_nexc )
  end subroutine

  !> Prints the number of excitations
  subroutine write_nexc( print_header, time_array, n_exc_array, n_gs_array )
    !> printHeader: If we need to print a header (useful when we open the file for the 1st time)
    logical, intent(in)   :: print_header
    !> Array with the values of time \( t \)
    real(dp), intent(in)  :: time_array(:)
    !> Number of electrons which were excited
    real(dp), intent(in)  :: n_exc_array(:)
    !> Number of electrons on the groundstate
    real(dp), intent(in)  :: n_gs_array(:)

    integer(i32) :: i, n
    character(len=*), parameter :: format_time = 'F9.3'
    character(len=*), parameter :: format_n = 'F20.10'
    character(len=*), parameter :: format_header = '(A9,3A20)'
    character(len=*), parameter :: format_line = '(' // format_time // ',3' // format_n // ')'

    n = size( time_array )
    call assert( size( n_exc_array ) == n, 'n_exc_array must have n elements')
    call assert( size( n_gs_array ) == n, 'n_gs_array must have n elements')

    if ( print_header ) write( file_nexc, format_header ) 'Time','N.Elec.GS', 'N.XS', 'Sum'
    do i = 1, n
      write( file_nexc, format_line ) time_array(i), n_gs_array(i), n_exc_array(i), n_gs_array(i)+n_exc_array(i)
    end do
  end subroutine

  subroutine open_file_info()
    open( newunit=file_info, file=add_default_extension( filename_info ), status="replace", action="write" )
  end subroutine

  subroutine close_file_info()
    close( file_info )
  end subroutine

  subroutine write_file_info( string, string_format )
    !> string to be printed out
    character(len=*), intent(in) :: string
    !> format for printing out `string`
    character(len=*), intent(in), optional :: string_format

    character(len=:), allocatable :: format_local
    if( present(string_format) ) then
      format_local = string_format
    else
      format_local = '(A)'
    end if
    write( file_info, format_local ) string
  end subroutine

  subroutine write_file_info_fill_line_with_char( ch )
    character, intent(in) :: ch
    call printline( file_info, ch )
  end subroutine

  subroutine write_file_info_header()
    character(len=str_128) :: string
    character(len=str_16) :: dat, tim

    call write_file_info( 'Real-time TDDFT calculation started' )
    call write_file_info( 'EXCITING '//trim( versionname )//' started' ) 
    if ( len( trim( githash ) ) > 0 ) call write_file_info( 'version hash id: '//trim( githash ) )
#ifdef MPI
    write( string, '(A,I6,A)') 'MPI version using ', procs, ' processor(s)'
    call write_file_info( string )
#endif
    call date_and_time( date=dat, time=tim )
    write( string, '("Date (DD-MM-YYYY) : ", A2, "-", A2, "-", A4)') &
    &  dat (7:8), dat (5:6), dat (1:4)
    call write_file_info( string )
    write( string, '("Time (hh:mm:ss)   : ", A2, ":", A2, ":", A2)') &
    &  tim (1:2), tim (3:4), tim (5:6)
    call write_file_info( string )
    call write_file_info( 'All units are atomic (Hartree, Bohr, etc.)' )
  end subroutine

  subroutine open_file_timing( new )
    !> If `.true.`, open new files (rewriting old ones).
    !> If `.false.`, append to existing files
    logical, intent(in) :: new
    open( newunit=file_time, file=add_default_extension(filename_timing), &
      status=get_status_from_logical( new ), position=get_position_from_logical( new ) )
  end subroutine

  subroutine close_file_timing()
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
    
    integer(i32) :: ip, shift, n
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
        call write_nonzero_timing( '-- basis:', t_rttddft%dens%basis )
        call write_nonzero_timing( 'updatepot:', t_rttddft%pot%total )
        call write_nonzero_timing( '-- poteff:', t_rttddft%pot%poteff )
        call write_nonzero_timing( '-- genveffig:', t_rttddft%pot%genveffig )
        call write_nonzero_timing( '-- genmeffig:', t_rttddft%pot%genmeffig )
        call write_nonzero_timing( 'UpdateCurrentDensity:', t_rttddft%current_density )
        call write_nonzero_timing( 'ObtainA:', t_rttddft%vector_potential )
        call write_nonzero_timing( 'Berry-phase related:', t_rttddft%td_berry )
        call write_nonzero_timing( 'updatehamiltonian:', t_rttddft%ham%total )
        call write_nonzero_timing( '-- hmlint:', t_rttddft%ham%hmlint )
        call write_nonzero_timing( 'All cycles of predcorr:', t_rttddft%pred_corr )
        call write_nonzero_timing( 'Total Energy:', t_rttddft%energy )
        call write_nonzero_timing( 'nexc:', t_rttddft%n_exc )
        call write_nonzero_timing( 'Screenshots:', t_rttddft%screenshot )
        call write_nonzero_timing( 'Print:', t_rttddft%t_print )
      end associate
      if( MD ) then
        associate( t_MD => timing(ip)%t_Ehrenfest )
          call write_nonzero_timing( 'MD:', t_MD%t_MD_step )
          call write_nonzero_timing( '-- 1st part of forces:', t_MD%t_MD_1st )
          call write_nonzero_timing( '-- 2nd part of forces:', t_MD%t_MD_2nd )
          call write_nonzero_timing( '-- sum forces:', t_MD%t_MD_sumforces )
          call write_nonzero_timing( '-- move ions:', t_MD%t_MD_moveions )
          call write_nonzero_timing( '-- update basis:', t_MD%t_MD_updateBasis )
          call write_nonzero_timing( '-- update H:', t_MD%ham )
          call write_nonzero_timing( '-- update S:', t_MD%overlap%total )
          call write_nonzero_timing( '-- update pmat:', t_MD%pmat )
        end associate
      end if
      write( file_time, format_timing ) 'time per iteration:', timing(ip)%t_iteration
    end do
  end subroutine

  !> Write timing only if it is nonzero (> tol)
  subroutine write_nonzero_timing( description, timing )
    !> Action name
    character(len = *), intent(in) :: description
    !> Action duration
    real(dp), intent(in) :: timing

    real(dp), parameter :: tol = 1.0e-6_dp

    if ( timing > tol ) write( file_time, format_timing ) description, timing

  end subroutine

  !> (Private) Get `filename_projection_coefficients` including the iteration number
  function get_filename_projection_coefficient( it ) result(name)
    !> Iteration number
    integer(i32), intent(in) :: it
    !> File name (to return)
    character(len=:), allocatable :: name

    name = add_default_extension( filename_projection_coefficients // to_char( it ) )
  end function

  !> Output projection coefficients to a text file
  subroutine write_projection_coefficients( it, print_absolute_value, print_format, proj_coeff )
    !> Iteration number
    integer(i32), intent(in) :: it
    !> If `.true.`, print `abs**2` of each projection coefficient instead of the complex number
    logical, intent(in) :: print_absolute_value
    !> Fortran format of a real number
    character(len=*), intent(in) :: print_format
    !> Projection coefficients
    complex(dp), contiguous, intent(in) :: proj_coeff(:, :, :)

    character(len=:), allocatable :: format_lines
    character(len=*), parameter :: format_header_line = '(A5,I10)'
    integer(i32) :: unit, ist, ik, n_states_gnd, last_kpt

    last_kpt = ubound( proj_coeff, 3 )
    n_states_gnd = size( proj_coeff, 1 )
    format_lines = '(' // trim( adjustl( to_char( merge( n_states_gnd, 2*n_states_gnd, print_absolute_value ) ) ) ) // trim( print_format ) // ')'

    open( newunit=unit, file = get_filename_projection_coefficient( it ), action = 'write' )
    do ik = 1, size( proj_coeff, 3 )
      write( unit, format_header_line) 'ik: ', ik
      do ist = 1, size( proj_coeff, 2 )
        if( print_absolute_value ) then
          write( unit, format_lines ) abs( proj_coeff(:, ist, ik) )**2
        else
          write( unit, format_lines ) proj_coeff(:, ist, ik)
        end if
      end do
    end do
    close( unit )
  end subroutine

  !> (Private) Get `filename_eigenvalues` including the iteration number
  function get_filename_eigenvalues( it ) result(name)
    !> Iteration number
    integer(i32), intent(in) :: it
    !> File name (to return)
    character(len=:), allocatable :: name

    name = add_default_extension( filename_eigenvalues // to_char( it ) )
  end function

  !> Output eigenvalues to a text file
  subroutine write_eigenvalues( it, eigenvalues, dimensions )
    !> Iteration number
    integer(i32), intent(in) :: it
    !> KS eigenvalues
    real(dp), contiguous, intent(in) :: eigenvalues(:, :)
    !> Size of eigenvalues (along 1st dim.) to be printed out
    integer(i32), contiguous, intent(in) :: dimensions(:)

    call write_array_along_states_and_kpoints( get_filename_eigenvalues( it ), eigenvalues, dimensions )
  end subroutine

  !> (Private) Get `filename_occupations` including the iteration number
  function get_filename_occupations( it ) result(name)
    !> Iteration number
    integer(i32), intent(in) :: it
    !> File name (to return)
    character(len=:), allocatable :: name

    name = add_default_extension( filename_occupations // to_char( it ) )
  end function

  !> Output eigenvalues to text and/or binary file(s)
  subroutine write_occupations( it, occupations, write_txt, write_binary, print_format )
    !> Iteration number
    integer(i32), intent(in) :: it
    !> Occupation factors
    real(dp), contiguous, intent(in) :: occupations(:, :)
    !> If `.true.`, write an output with text format
    logical, intent(in) :: write_txt
    !> If `.true.`, write an output with binary format
    logical, intent(in) :: write_binary
    !> Fortran format of a real number
    character(len=*), intent(in) :: print_format

    integer(i32) :: ik
    character(len=:), allocatable :: backup

    if( write_binary ) then
      backup = filext
      filext = "_" // to_char(it) // filext
      do ik = 1, size( occupations, 2 )
        call putoccsv(ik, occupations(:, ik))
      end do
      filext = backup
    end if
    if( write_txt ) call write_array_along_states_and_kpoints( get_filename_occupations(it), &
      occupations, real_number_format=print_format )
  end subroutine

  !> (Private) Auxiliary function to write eigenvalues and occupations
  subroutine write_array_along_states_and_kpoints( file_name, array, dimensions, real_number_format )
    !> File name
    character(len=*), intent(in) :: file_name
    !> Array to be printed out
    real(dp), contiguous, intent(in) :: array(:, :)
    !> List with sizes of `array` (along 1st dim.) to print out
    integer(i32), contiguous, optional, intent(in) :: dimensions(:)
    !> Format of a real number. When absent, `real_number_format_default` is used
    character(len=*), optional, intent(in) :: real_number_format
    
    integer(i32) :: unit, ik, ist, m, n
    integer(i32), allocatable :: dims(:)
    character(len=*), parameter :: integer_format = 'I7'
    character(len=*), parameter :: real_number_format_default = 'F20.12'
    character(len=*), parameter :: format_header_line = '(A5,' // integer_format // ')'
    character(len=:), allocatable :: format_lines
    
    m = size( array, 1 )
    n = size( array, 2 )
    if( present(dimensions) ) then
      call assert( size(dimensions) == n, "dimensions must have size n")
      call assert( all( dimensions <= m ), "each element in dimensions must be <= m" )
      dims = dimensions
    else 
      allocate( dims(n), source=m )
    end if

    if( present(real_number_format) ) then
      format_lines = '(' // integer_format // ',' // trim(real_number_format) // ')'
    else
      format_lines = '(' // integer_format // ',' // trim(real_number_format_default) // ')'
    end if
    
    open( newunit = unit, file = file_name, action = 'write' )
    do ik = 1, n
      write( unit, format_header_line ) 'ik = ', ik
      do ist = 1, dims(ik)
        write( unit, format_lines ) ist, array(ist, ik)
      end do
      write( unit, * ) ''
    end do
    close(unit)
  end subroutine

  !> Return a label, indicating `density` or `delta-density`
  pure function get_density_label( delta_rho ) result(label)
    !> If `.true.`, select `delta-density`
    logical, intent(in) :: delta_rho
    !> File name (to return)
    character(len=:), allocatable :: label

    if( delta_rho ) then
      label = filename_density_changes
    else 
      label = filename_density
    end if
  end function

  !> (Private) Return the output name where the density will be stored
  function get_filename_density( it, delta_rho ) result(name)
    !> Iteration number
    integer(i32), intent(in) :: it
    !> If `.true.`, select `delta-density`
    logical, intent(in) :: delta_rho
    !> File name (to return)
    character(len=:), allocatable :: name

    name = add_xsf_extension( get_density_label( delta_rho ) // "_" // to_char( it )  )
  end function

  !> Write the electron density (or changes in electron density) to an output file
  subroutine write_density_to_file( it, rho_MT, rho_interstitial, delta_rho, plot3d, my_rank_writes )
    !> Iteration number
    integer, intent(in) :: it
    !> Real-valued function in MT (lmmaxvr, nrmtmax, natmtot): it can be \(n\) or \(\Delta n\)
    real(dp), contiguous, intent(in) :: rho_MT(:, :, :)
    !> Real-valued function in IR region (ngrtot): it can be \(n\) or \(\Delta n\)
    real(dp), contiguous, intent(in) :: rho_interstitial(:)
    !> If `.true.`, `rho_MT` and `rho_insterstitial` refer to \(\Delta n\), instead of \(n\)
    logical, intent(in) :: delta_rho
    !> Grid data fot 3D density plots
    type(plot3d_type), intent(in), pointer :: plot3d
    !> If `.true.`, my MPI rank is supposed to write outputs
    logical, intent(in) :: my_rank_writes

    integer(i32) :: l_max, lm_max
    real(dp), allocatable :: rho_rgrid(:)
    type(rgrid) :: grid_3D

    grid_3D = gen_3d_rgrid( plot3d, 0 )
    allocate( rho_rgrid(grid_3D%npt) )
    lm_max = size( rho_MT, 1 )
    l_max = int( sqrt( real(lm_max, dp) ), i32 ) - 1
    ! Attention, this function has an MPI collective call
    call transform_real_function_to_rgrid( grid_3D, l_max, rho_MT, rho_interstitial, rho_rgrid )
    if( my_rank_writes ) call write_real_function_xsf( grid_3D, rho_rgrid, get_density_label( delta_rho ) , get_filename_density( it, delta_rho ))
  end subroutine

  !> Convert real-valued coordinate-space-defined 3d function from the 
  !> IR-MT representation to the coordinate representation
  subroutine transform_real_function_to_rgrid( grid, lmax, function_mt, function_ir, function_rgrid )
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

    call rfarray( lmax, size( function_mt, 1 ), function_mt, function_ir, grid%npt, grid%vpl, function_rgrid )
  end subroutine 

end module