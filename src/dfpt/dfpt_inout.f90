!> This module contains i/o functionalities for DFPT calculations.
module dfpt_inout
  use dfpt_variables

  use precision, only: dp
  use FoX_wxml

  implicit none
  private

  ! GLOBAL VARIABLES
  !> human readable info output unit
  integer, public :: info_out_unit = 0
  !> xml info output file
  type(xmlf_t), public :: info_out_xml
  !> upper case prefix for file names
  character(64), public :: prefix_upper = 'DFPT'
  !> lower case prefix for file names
  character(64), public :: prefix_lower = 'dfpt'
  !> info file is opened
  logical, public :: info_is_open = .false.

  public :: dfpt_io_set_prefix
  public :: dfpt_io_read_zfun, dfpt_io_write_zfun
  public :: dfpt_io_info_init, dfpt_io_info_finit, dfpt_io_info_scf_init, dfpt_io_info_scf_finit, dfpt_io_info_scf, dfpt_io_info_string

  contains

    !> This subroutine sets the standard prefix for output files.
    subroutine dfpt_io_set_prefix( upper, lower )
      !> upper/lower case prefix
      character(*), intent(in) :: upper, lower
      if( trim( upper ) /= '' ) then
        write( prefix_upper, '(a)' ) trim( upper )
      else
        prefix_upper = ''
      end if
      if( trim( lower ) /= '' ) then
        write( prefix_lower, '(a)' ) trim( lower )
      else
        prefix_lower = ''
      end if
    end subroutine dfpt_io_set_prefix

    !> This subroutine writes a complex unit cell function to a binary file.
    !> See also [[dfpt_io_read_zfun(subroutine)]].
    subroutine dfpt_io_write_zfun( zfmt, zfir, fname, success, &
        file_extension, directory )
      use m_getunit
      !> complex muffin-tin function given as a complex spherical harmonics expansion
      complex(dp), intent(in) :: zfmt(:,:,:)
      !> complex interstitial function on a real space grid
      complex(dp), intent(in) :: zfir(:)
      !> file name
      character(*), intent(in) :: fname
      !> `.true.` if reading was successful
      logical, intent(out) :: success
      !> file extension
      character(*), optional, intent(in) :: file_extension
      !> path to directory (default: current directory)
      character(*), optional, intent(in) :: directory

      character(64) :: fxt
      character(256) :: dirname
      integer :: stat, un, lmmax, nrmtmax, natmtot, ngrtot

      fxt = ''
      if( present( file_extension ) ) write( fxt, '("_",a)' ) trim( file_extension )
      dirname = '.'
      if( present( directory ) ) write( dirname, '(a)' ) trim( directory )
      dirname = trim( dirname )//'/'

      lmmax   = size( zfmt, dim=1 )  ! maximum number of (l,m) pairs
      nrmtmax = size( zfmt, dim=2 )  ! maximum number or radial points
      natmtot = size( zfmt, dim=3 )  ! total number of atoms
      ngrtot  = size( zfir, dim=1 )  ! total number of real-space grid points

      call getunit( un )
      open( un, file=trim( dirname )//trim( fname )//trim( fxt )//'.OUT', action='write', form='unformatted', iostat=stat )
      success = (stat == 0)
      if( .not. success ) return
      write( un, iostat=stat ) lmmax, nrmtmax, natmtot, ngrtot
      success = success .and. (stat == 0)
      write( un, iostat=stat ) zfmt, zfir
      success = success .and. (stat == 0)
      close( un )
    end subroutine dfpt_io_write_zfun

    !> This subroutine reads a complex unit cell function from a binary file.
    !> See also [[dfpt_io_write_zfun(subroutine)]].
    subroutine dfpt_io_read_zfun( zfmt, zfir, fname, success, &
        file_extension, directory, error_message )
      use m_getunit
      !> complex muffin-tin function given as a complex spherical harmonics expansion
      complex(dp), intent(out) :: zfmt(:,:,:)
      !> complex interstitial function on a real space grid
      complex(dp), intent(out) :: zfir(:)
      !> file name
      character(*), intent(in) :: fname
      !> `.true.` if reading was successful
      logical, intent(out) :: success
      !> file extension
      character(*), optional, intent(in) :: file_extension
      !> path to directory (default: current directory)
      character(*), optional, intent(in) :: directory
      !> error message in case of unsuccessful reading
      character(*), optional, intent(out) :: error_message

      character(64) :: fxt
      character(256) :: dirname
      character(512) :: errmsg
      integer :: stat, un, lmmax, nrmtmax, natmtot, ngrtot

      fxt = ''
      if( present( file_extension ) ) write( fxt, '("_",a)' ) trim( file_extension )
      dirname = '.'
      if( present( directory ) ) write( dirname, '(a)' ) trim( directory )
      dirname = trim( dirname )//'/'

      call getunit( un )
      open( un, file=trim( dirname )//trim( fname )//trim( fxt )//'.OUT', action='read', form='unformatted', iostat=stat )
      success = (stat == 0)
      if( .not. success ) return
      read( un, iostat=stat ) lmmax, nrmtmax, natmtot, ngrtot
      success = success .and. (stat == 0)

      errmsg = '(dfpt_io_read_zfun)'//new_line( 'a' )//'file: '//trim( dirname )//trim( fname )//trim( fxt )//'.OUT'
      if( lmmax /= size( zfmt, dim=1 ) ) then
        write( errmsg, '(a,"`lmmax` of function in file (",i4",) and in argument function (",i4,") differ.")' ) &
          trim( errmsg )//new_line( 'a' ), lmmax, size( zfmt, dim=1 )
        success = .false.
      end if
      if( nrmtmax /= size( zfmt, dim=2 ) ) then
        write( errmsg, '(a,"`nrmtmax` of function in file (",i4",) and in argument function (",i4,") differ.")' ) &
          trim( errmsg )//new_line( 'a' ), nrmtmax, size( zfmt, dim=2 )
        success = .false.
      end if
      if( natmtot /= size( zfmt, dim=3 ) ) then
        write( errmsg, '(a,"`natmtot` of function in file (",i4",) and in argument function (",i4,") differ.")' ) &
          trim( errmsg )//new_line( 'a' ), natmtot, size( zfmt, dim=3 )
        success = .false.
      end if
      if( ngrtot /= size( zfir, dim=1 ) ) then
        write( errmsg, '(a,"`ngrtot` of function in file (",i4",) and in argument function (",i4,") differ.")' ) &
          trim( errmsg )//new_line( 'a' ), ngrtot, size( zfir, dim=1 )
        success = .false.
      end if

      if( .not. success ) then
        close( un )
        if( present( error_message ) ) error_message = trim( errmsg )
        return
      end if
      read( un, iostat=stat ) zfmt, zfir
      success = success .and. (stat == 0)
      if( .not. success .and. present( error_message ) ) &
        error_message = '(dfpt_io_read_zfun): Uknown error in reading function from file.'
      close( un )
    end subroutine dfpt_io_read_zfun

    !> This subroutine opens and initializes the general info output file (human readable and xml).
    subroutine dfpt_io_info_init( &
        file_extension, directory )
      use m_getunit
      use mod_misc, only: versionname, githash
      use modinput
      !> file extension
      character(*), optional, intent(in) :: file_extension
      !> path to directory (default: current directory)
      character(*), optional, intent(in) :: directory

      character(10) :: date, time
      character(64) :: fxt
      character(256) :: dirname
      character(512) :: string
      integer :: stat

      fxt = ''
      if( present( file_extension ) ) write( fxt, '("_",a)' ) trim( file_extension )
      dirname = '.'
      if( present( directory ) ) write( dirname, '(a)' ) trim( directory )
      dirname = trim( dirname )//'/'

      ! create files
      ! human readable
      call getunit( info_out_unit )
      open( info_out_unit, file=trim( dirname )//trim( prefix_upper )//'_INFO'//trim( fxt )//'.OUT', action='write', form='formatted', iostat=stat )
      info_is_open = (stat == 0)
      ! xml
      call xml_OpenFile( trim( dirname )//trim( prefix_lower )//'_info'//trim( fxt )//'.xml', info_out_xml, replace=.true., pretty_print=.true. )
      call xml_AddXMLPI( info_out_xml, 'xml-stylesheet', 'href="'//trim( input%xsltpath )//'/info.xsl" type="text/xsl"' )

      ! write version and time
      call date_and_time( date=date, time=time )
      ! human readable
      call printline( info_out_unit, '=' )
      write( string, '("EXCITING ", a, " DFPT ", a," calculation started")' ) trim( versionname ), trim( prefix_upper )
      call printtext( info_out_unit, '=', trim( string ) )
      if( len( trim( githash ) ) > 0 ) then
        write( string, '("version hash id: ",a)' ) githash
        call printtext( info_out_unit, '=', trim( string ) )
      end if
      call printtext( info_out_unit, '=', '' )
      write( string, '("Date (DD-MM-YYYY) : ",A2,"-",A2,"-",A4)' ) date(7:8), date(5:6), date(1:4)
      call printtext( info_out_unit, '=', trim( string ) )
      write( string, '("Time (hh:mm:ss)   : ",A2,":",A2,":",A2)' ) time(1:2), time(3:4), time(5:6)
      call printtext( info_out_unit, '=', trim( string ) )
      call printline( info_out_unit, '=' )
      ! xml
      call xml_NewElement( info_out_xml, 'info' )
      call xml_AddAttribute( info_out_xml, 'versionname', trim( versionname ) )
      call xml_AddAttribute( info_out_xml, 'versionhash', trim( githash ) )
      call xml_AddAttribute( info_out_xml, 'title', trim( input%title ) )
      write( string, '(A2,"-",A2,"-",A4)' ) date(7:8), date(5:6), date(1:4)
      call xml_AddAttribute( info_out_xml, 'date', trim( string ) )
      write( string, '(A2,":",A2,":",A2)' ) time(1:2), time(3:4), time(5:6)
      call xml_AddAttribute( info_out_xml, 'time', trim( string ) )
    end subroutine dfpt_io_info_init

    !> This subroutine finalizes and closes the general info output file (human readable and xml).
    subroutine dfpt_io_info_finit
      use mod_misc, only: versionname

      integer :: stat
      character(512) :: string

      if( .not. info_is_open ) return

      ! human readable
      call printline( info_out_unit, '=' )
      write( string, '("EXCITING ", a, " DFPT ", a," calculation stopped")' ) trim( versionname ), trim( prefix_upper )
      call printtext( info_out_unit, '=', trim( string ) )
      call printline( info_out_unit, '=' )
      close( info_out_unit, iostat=stat )
      info_is_open = (stat /= 0)
      ! xml
      call xml_EndElement( info_out_xml, 'info' )
      call xml_close( info_out_xml )
    end subroutine dfpt_io_info_finit

    !> This subroutine initializes the self-consistency loop output in the general 
    !> info outout file (human readable and xml).
    subroutine dfpt_io_info_scf_init( fromfile )
      !> density response initialized from file
      logical, intent(in) :: fromfile

      call printbox( info_out_unit, '*', 'scf loop for Sternheimer equation started' )
      call xml_NewElement( info_out_xml, 'scf' )
      if( fromfile ) then
        write( info_out_unit, '("Initial density response read from file.")' )
      end if
      write( info_out_unit, * )
    end subroutine dfpt_io_info_scf_init

    !> This subroutine finalizes the self-consistency loop output in the general 
    !> info outout file (human readable and xml).
    subroutine dfpt_io_info_scf_finit
      call printbox( info_out_unit, '*', 'scf loop for Sternheimer equation stopped' )
      call xml_EndElement( info_out_xml, 'scf' )
      write( info_out_unit, * )
    end subroutine dfpt_io_info_scf_finit

    !> This subroutine writes an iteration of the self-consistency cycle to the general
    !> info outout file (human readable and xml).
    subroutine dfpt_io_info_scf( iiter, delta, time, &
        defermi )
      use modinput
      !> iteration number
      integer, intent(in) :: iiter
      !> maximum change in potential/density response
      real(dp), intent(in) :: delta
      !> time spent for iteration
      real(dp), intent(in) :: time
      !> Fermi energy response
      real(dp), optional, intent(in) :: defermi

      character(512) :: string
      
      write( info_out_unit, '("Sternheimer SCF loop iteration ",i3)' ) iiter
      call xml_NewElement( info_out_xml, 'iteration' )
      write( string, '(i8)' ) iiter
      call xml_AddAttribute( info_out_xml, 'counter', trim( adjustl( string ) ) )
      write( string, '(g25.16)' ) delta
      if( input%groundstate%mixerswitch == 1 ) then
        write( info_out_unit, '("maximum change in potential response:",T50,g16.6)' ) delta
        call xml_AddAttribute( info_out_xml, 'potential_change', trim( adjustl( string ) ) )
      else
        write( info_out_unit, '("maximum change in density response:",T50,g16.6)' ) delta
        call xml_AddAttribute( info_out_xml, 'density_change', trim( adjustl( string ) ) )
      end if
      if( present( defermi ) ) then
        write( string, '(g25.16)' ) defermi
        write( info_out_unit, '("Fermi energy response:",T50,g16.6)' ) defermi
        call xml_AddAttribute( info_out_xml, 'Fermi_energy_response', trim( adjustl( string ) ) )
      end if
      write( info_out_unit, '("Time spent for iteration (seconds):",T50,f16.2)' ) time
      write( string, '(g25.16)' ) time
      call xml_AddAttribute( info_out_xml, 'duration', trim( adjustl( string ) ) )
      write( info_out_unit, * )
      call xml_EndElement( info_out_xml, 'iteration' )
    end subroutine dfpt_io_info_scf

    !> Write string to human readable output file.
    subroutine dfpt_io_info_string( string )
      !> string to write
      character(*), intent(in) :: string

      write( info_out_unit, '(a)' ) trim( adjustl( string ) )//new_line( 'a' )
    end subroutine dfpt_io_info_string

end module dfpt_inout
