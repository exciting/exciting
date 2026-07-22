module gw_info
#include "asserts.fpp"
  use gw_io, only: open_file, file_format_text
  use mod_mpi_gw, only: indexes_parallelization
  use precision, only: dp, i32, str_64, str_128

  implicit none

  private

  ! Strings informing progress status 
  character(len=*), public, parameter :: progress_coulomb = 'coulomb'
  character(len=*), public, parameter :: progress_sgi = 'sgi'
  character(len=*), public, parameter :: progress_ipw = 'ipw'
  character(len=*), public, parameter :: progress_calculate_epsilon = 'calc_eps'
  character(len=*), public, parameter :: progress_polarizability_tetrahedron = 'qdepwtet'
  character(len=*), public, parameter :: progress_write_epsilon = 'write_eps'
  character(len=*), public, parameter :: progress_epsilon_loop_kpoints = 'eps_loop_kpt'
  character(len=*), public, parameter :: progress_epsilon_loop_blocks = 'eps_loop_blocks'

  ! Strings to name symbols
  character(len=*), public, parameter :: q_points_numbering = 'dummy index to list q-points'
  character(len=*), public, parameter :: q_points_numbering_abbr = '(xq)'
  character(len=*), public, parameter :: q_points_indexes = 'q-points indexes (as used in the code)'
  character(len=*), public, parameter :: q_points_indexes_abbr = 'q'
  character(len=*), public, parameter :: k_points_numbering = 'dummy index to list k-points'
  character(len=*), public, parameter :: k_points_numbering_abbr = '(xk)'
  character(len=*), public, parameter :: k_points_indexes = 'k-points indexes (as used in the code)'
  character(len=*), public, parameter :: k_points_indexes_abbr = 'k'
  character(len=*), public, parameter :: empty_bands_indexes = 'unoccupied states'
  character(len=*), public, parameter :: empty_bands_abbr = 'm'
  character(len=*), public, parameter :: bands_indexes = 'states (may include core states, depending on the choice given in input.xml)'
  character(len=*), public, parameter :: bands_abbr = 'm'

  !> Unit used for the general GW output file (`GW_INFO.OUT`)
  integer(i32), public, protected :: fgw
  !> Default name of the general GW output file
  character(len=*), parameter :: filename_gwinfo = 'GW_INFO.OUT'

  public :: open_gwinfo, write_to_gwinfo, write_to_gwinfo_boxmessage, &
    write_to_gwinfo_parallelization_info, write_to_gwinfo_progress, &
    write_to_gwinfo_progress_bar, write_to_gwinfo_table_with_index_map

contains

!> Open the `GW_INFO.OUT`
subroutine open_gwinfo( )
  
  call open_file( filename_gwinfo , 'write', file_format_text, fgw )

end subroutine


!> Write a string into `GW_INFO.OUT`
subroutine write_to_gwinfo( string )
  character(len=*), intent(in) :: string

  write( fgw, * ) string

end subroutine


!> Write information about the progress of subtasks in a `taskGroup` calculation
subroutine write_to_gwinfo_progress( string )
  character(len=*), intent(in) :: string
  
  select case( trim(string) )
    case( progress_coulomb )
      call write_to_gwinfo('*** Reading Coulomb matrix from file')
    case( progress_sgi )
      call write_to_gwinfo('*** Reading matrix with orthogonal IPW from file')
    case( progress_ipw )
      call write_to_gwinfo('*** Obtaining overlap between IPW and PW')
    case( progress_calculate_epsilon )
      call write_to_gwinfo('*** Calculating dielectric matrix')
    case( progress_polarizability_tetrahedron )
      call write_to_gwinfo('****** Evaluating polarizability factors with the tetrahedron method')
    case( progress_epsilon_loop_kpoints )
      call write_to_gwinfo('****** Loop over k-points')
    case( progress_epsilon_loop_blocks )
      call write_to_gwinfo('********* Loop over block-multiplications')
    case( progress_write_epsilon )
      call write_to_gwinfo('*** Writing dielectric matrix into file')
    case default
      CALL_ASSERT(.false., 'Unknown option in write_to_gwinfo_progress')
  end select

end subroutine


!> Write information about the progress of subtasks in a `taskGroup` calculation
subroutine write_to_gwinfo_progress_bar( time_spent, fraction, identation_level )
  !> Time in sec. spent in the current task
  real(dp), intent(in) :: time_spent
  !> Fraction (number between 0 and 1) of current task that has been completed
  real(dp), intent(in) :: fraction
  !> Define the identation level of the progress bar
  integer(i32), intent(in) :: identation_level
  
  character(len=str_128) :: string
  character(len=*), parameter :: format_str = '(1A,1F6.2,1A,1F15.1,1A,1F15.1,1A)'
  character(len=:), allocatable :: identation_string
  integer(i32), parameter :: identation_multiplier = 3
  character(len=1), parameter :: identation_char = '*'

  identation_string = repeat( identation_char, identation_multiplier*identation_level )
  CALL_ASSERT( fraction >= 0 .and. fraction <= 1, message='fraction must be a real number between 0 and 1')
  CALL_ASSERT( identation_level >= 0, 'identation_level must be positive' )
  write( string, format_str ) identation_string // ' completing ', 100*fraction, '%. Time: ', time_spent, &
    ' (sec). Estimated total time:', time_spent/fraction, ' (sec)'
  call write_to_gwinfo( trim(string) )

end subroutine


!> Write a table containing an index map to `GW_INFO.OUT`
subroutine write_to_gwinfo_table_with_index_map( index_map, symbols, names )
  !> index map to be written to the output file
  integer(i32), contiguous, intent(in) :: index_map(:)
  !> symbols (abbreviations) that define physical quantities represented by the index map
  character(len=*), contiguous, intent(in) :: symbols(:)
  !> full names of the physical quantities represented by the sets of indexes
  character(len=*), contiguous, intent(in) :: names(:)
  
  integer(i32) :: i, n
  character(len=str_64) :: string
  character(len=*), parameter :: spacing = '8'
  character(len=*), parameter :: format_int = '(2I' // spacing // ')'
  character(len=*), parameter :: format_str = '(2A' // spacing // ')'

  CALL_ASSERT( size( symbols ) == 2, 'Array symbols must have size 2' )
  CALL_ASSERT( size( names ) == 2, 'Array names must have size 2' )

  n = size( index_map )
  call write_to_gwinfo( '*** Table with index map' )
  call write_to_gwinfo( '*** 1st quantity is: ' // trim( names(1) ) // ', with symbol: ' // trim( symbols(1) ) )
  call write_to_gwinfo( '*** 2nd quantity is: ' // trim( names(2) ) // ', with symbol: ' // trim( symbols(2) ) )
  call write_to_gwinfo( '' )

  call write_to_gwinfo( 'Map' )
  write( string, format_str ) trim( symbols(1) ), trim( symbols(2) )
  call write_to_gwinfo( string )

  do i = 1, n
    write( string, format_int ) i, index_map(i)
    call write_to_gwinfo( string )
  end do

  call write_to_gwinfo( '' )
  
end subroutine


!> Write information about parallelization to `GW_INFO.OUT`
subroutine write_to_gwinfo_parallelization_info( ranks, sets, symbols, names )
  !> Array of ranks
  integer(i32), contiguous, intent(in) :: ranks(:)
  !> Array of indexes representing how m physical quantities (this is size along the 1st. dim.)
  !> are distributed over n ranks (n is the size along the 2nd dim.)
  type(indexes_parallelization), contiguous, intent(in) :: sets(:, :)
  !> symbols (abbreviations) that define physical quantities represented by the sets of indexes
  character(len=*), contiguous, intent(in) :: symbols(:)
  !> full names of the physical quantities represented by the sets of indexes
  character(len=*), contiguous, intent(in) :: names(:)

  integer(i32) :: i, j, n, m
  character(len=*), parameter :: spacing = '8'
  character(len=*), parameter :: format_integer = '(1I' // spacing // ')'
  character(len=*), parameter :: format_string = '(1A' // spacing // ')'

  m = size( symbols, 1 )
  n = size( ranks )
  CALL_ASSERT( size(sets, 1) == m, 'sets must have m components along 1st dim'  )
  CALL_ASSERT( size(sets, 2) == n, 'sets must have n components along 2nd dim' ) 
  CALL_ASSERT( size(names, 1) == m, 'names and symbols must have same size'  )
  
  call write_to_gwinfo('')
  call write_to_gwinfo('*** Parallelization')
  call write_to_gwinfo('*** Overview of how indexes are distributed among MPI processes')
  do j = 1, m
    write( fgw, * ) '*** '// trim( adjustl( symbols(j) ) ) //'i, '// &
      trim( adjustl( symbols(j) ) )//'f : first, last ' // trim( adjustl( names(j) ) )
  end do
  write( fgw, format_string, advance='no' ) 'global'
  do j = 1, m
    write( fgw, format_string, advance='no' ) trim( adjustl( symbols(j) ) ) //'i' 
    write( fgw, format_string, advance='no' ) trim( adjustl( symbols(j) ) ) //'f' 
  end do
  call write_to_gwinfo('')

  write( fgw, format_string, advance='no' ) '---'
  do j = 1, m
    write( fgw, format_integer, advance='no' ) minval( [ (sets(j, i)%my_first, i = 1, n) ])
    write( fgw, format_integer, advance='no' ) maxval( [ (sets(j, i)%my_last, i = 1, n) ])
  end do
  write( fgw, '(A,/,A)' ) '',''

  write( fgw, format_string, advance='no' ) 'rank'
  do j = 1, m
    write( fgw, format_string, advance='no' ) trim( adjustl( symbols(j) ) ) //'i' 
    write( fgw, format_string, advance='no' ) trim( adjustl( symbols(j) ) ) //'f' 
  end do
  call write_to_gwinfo('')
  
  do i = 1, n
    write( fgw, format_integer, advance='no' ) ranks(i)
    do j = 1, m
      write( fgw, format_integer, advance='no' ) sets(j, i)%my_first
      write( fgw, format_integer, advance='no' ) sets(j, i)%my_last
    end do
    write( fgw, * ) ''
  end do
  write( fgw, '(A,/,A)' ) '',''

end subroutine


!> Write a string surrounded by a box of characters into `GW_INFO.OUT`
subroutine write_to_gwinfo_boxmessage( char, string )
  character, intent(in) :: char
  character(len=*), intent(in) :: string

  call BoxMSG( fgw, char, string )

end subroutine

end module