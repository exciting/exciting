!> Collection if I/O routines for electron-phonon calculations.
module eph_inout
  use precision, only: dp
  use modmpi, only: terminate_if_false
  use asserts, only: assert

  implicit none
  private

  public :: eph_io_write_energies, eph_io_write_ephmat
  
contains

  !> Write energies at given reciprocal space points to file.
  !> 
  !> Points can be given as a list of vectors or as a path object.
  subroutine eph_io_write_energies( bvec, e, fname, format, plist, path )
    use bz_path, only: bz_path_type
    use xjson, only: to_json
    !> reciprocal lattice vectors
    real(dp), intent(in) :: bvec(3, 3)
    !> energies
    real(dp), intent(in) :: e(:,:)
    !> file name (without suffix!)
    character(len=*), intent(in) :: fname
    !> output format;   
    !> currently supported: `text` (plain text), `json` (JSON dictionary)
    character(len=*), intent(in) :: format
    !> list of vectors representing reciprocal space points
    real(dp), optional, intent(in) :: plist(:,:)
    !> BZ path
    type(bz_path_type), optional, intent(in) :: path
  
    integer :: un, nst, np, ist, ip, stat
    character(len=8) :: sfx

    nst = size( e, dim=1 )
    np = size( e, dim=2 )

    ! check input
    call assert( .not. (present(plist) .and. present(path)), &
      'Not both `plist` and `path` should be present.' )

    sfx = ''
    select case (trim( adjustl( format ) ))
      case ('text')
        sfx = '.dat'
      case ('json')
        sfx = '.json'
      case default
        call terminate_if_false( .false., '(eph_io_write_energies) &
          Unsupported format `'//trim( adjustl( format ) )//'`.' )
    end select

    if (present(plist)) then
      call terminate_if_false( size( plist, dim=2 ) == np, '(eph_io_write_energies) &
        Energies `e` and point list `plist` must contain same number of points.' )
      ! point list text output
      if (sfx == '.dat') then
        open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
        call terminate_if_false( stat == 0, '(eph_io_write_energies) &
          Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
        do ip = 1, np
          write( un, '("#",100g26.16e3)' ) plist(:, ip)
          do ist = 1, nst
            write( un, '(g26.16e3)' ) e(ist, ip)
          end do
          write( un, * )
        end do
        close( un )
      ! point list JSON output
      else if (sfx == '.json') then
        open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
        call terminate_if_false( stat == 0, '(eph_io_write_energies) &
          Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
        write( un, '("{",a,": ",a,", ")', advance='no' ) '"bvec"', to_json( bvec )
        write( un, '(a,": ",a,", ")', advance='no' ) '"points"', to_json( plist )
        write( un, '(a,": ",a,"}")' ) '"energies"', to_json( e )
        close( un )
      end if
    else if (present(path)) then
      call terminate_if_false( path%num_points == np, '(eph_io_write_energies) &
        Energies `e` and path `path` must contain same number of points.' )
      ! path text output
      if (sfx == '.dat') then
        open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
        call terminate_if_false( stat == 0, '(eph_io_write_energies) &
          Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
        do ist = 1, nst
          do ip = 1, np
            write( un, '(2g26.16e3)' ) path%points(ip)%distance, e(ist, ip)
          end do
          write( un, * )
        end do
        close( un )
      ! path JSON output
      else if (sfx == '.json') then
        open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat)
        call terminate_if_false( stat == 0, '(eph_io_write_energies) &
          Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
        write( un, '("{",a,": ",a,", ")', advance='no' ) '"path"', path%to_json()
        write( un, '(a,": ",a,"}")' ) '"bands"', to_json( e )
        close( un )
      end if
    end if
  end subroutine eph_io_write_energies

  !> Write electron-phonon matrix elements at given reciprocal space points to file.
  subroutine eph_io_write_ephmat( bvec, vkl, vkql, vql, elengyk, elengykq, phengyq, g, gavg, fname, format )
    use xjson, only: to_json
    !> reciprocal lattice vectors
    real(dp), intent(in) :: bvec(3, 3)
    !> list of \({\bf k}\)-vectors
    real(dp), intent(in) :: vkl(:,:)
    !> list of \({\bf k}+{\bf q}\)-vectors
    real(dp), intent(in) :: vkql(:,:,:)
    !> list of \({\bf q}\)-vectors
    real(dp), intent(in) :: vql(:,:)
    !> electron energies \(\epsilon_{n{\bf k}}\)
    real(dp), intent(in) :: elengyk(:,:)
    !> electron energies \(\epsilon_{m{\bf k}+{\bf q}}\)
    real(dp), intent(in) :: elengykq(:,:,:)
    !> phonon frequencies \(\omega_{\nu{\bf q}}\)
    real(dp), intent(in) :: phengyq(:,:)
    !> matrix elements \(g_{mn,\nu}({\bf k},{\bf q})\)
    complex(dp), intent(in) :: g(:,:,:,:,:)
    !> absolute value of \(g\) averaged over degenerate states
    real(dp), intent(in) :: gavg(:,:,:,:,:)
    !> file name (without suffix!)
    character(len=*), intent(in) :: fname
    !> output format;   
    !> currently supported: `text` (plain text), `json` (JSON dictionary)
    character(len=*), intent(in) :: format

    integer :: un, nk, nq, nstk, nstkq, nmode, ik, iq, ist, jst, imode, stat
    character(len=8) :: sfx
    
    nk = size( vkl, dim=2 )
    nq = size( vql, dim=2 )
    nstk = size( elengyk, dim=1 )
    nstkq = size( elengykq, dim=1 )
    nmode = size( phengyq, dim=1 )

    ! check input
    call assert( size( elengyk, dim=2 ) == nk, '2nd dimension of `elengyk` must equal 2nd dimension of `vkl`.' )
    call assert( size( elengykq, dim=2 ) == nk, '2nd dimension of `elengykq` must equal 2nd dimension of `vkl`.' )
    call assert( size( elengykq, dim=3 ) == nq, '3rd dimension of `elengykq` must equal 2nd dimension of `vql`.' )
    call assert( size( phengyq, dim=2 ) == nq, '2nd dimension of `phengyq` must equal 2nd dimension of `vql`.' )
    call assert( size( g, dim=1 ) == nstkq, '1st dimension of `g` must equal 1st dimension of `elengykq`.' )
    call assert( size( g, dim=2 ) == nstk, '2nd dimension of `g` must equal 1st dimension of `elengyk`.' )
    call assert( size( g, dim=3 ) == nmode, '3rd dimension of `g` must equal 1st dimension of `phengyq`.' )
    call assert( size( g, dim=4 ) == nk, '4th dimension of `g` must equal 2nd dimension of `vkl`.' )
    call assert( size( g, dim=5 ) == nq, '5th dimension of `g` must equal 2nd dimension of `vql`.' )
    call assert( all( shape(g) == shape(gavg) ), '`g` and `gavg` must have same shape.' )

    sfx = ''
    select case (trim( adjustl( format ) ))
      case ('text')
        sfx = '.dat'
      case ('json')
        sfx = '.json'
      case default
        call terminate_if_false( .false., '(eph_io_write_ephmat) &
          Unsupported format `'//trim( adjustl( format ) )//'`.' )
    end select

    open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
    call terminate_if_false( stat == 0, '(eph_io_write_ephmat) &
      Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )

    ! text output
    if (sfx == '.dat') then
      write( un, '("#",a8,2a9,6a26)' ) '<m,k+q|', '|n,k>', 'nu', 'e_n(k)', 'e_m(k+q)', 'w_nu(q)', 'Re{g_mn,nu(k,q)}', 'Im{g_mn,nu(k,q)}', '|g_mn,nu(k,q)|'
      do iq = 1, nq
        do ik = 1, nk
          write( un, '("# vkl  = ",3g26.16e3)' ) vkl(:, ik)
          write( un, '("# vql  = ",3g26.16e3)' ) vql(:, iq)
          write( un, '("# vkql = ",3g26.16e3)' ) vkql(:, ik, iq)
          do imode = 1, nmode
            do jst = 1, nstk
              do ist = 1, nstkq
                write( un, '(3i9,6g26.16e3)' ) ist, jst, imode, &
                  elengykq(ist, ik, iq), elengyk(jst, ik), phengyq(imode, iq), &
                  g(ist, jst, imode, ik, iq), gavg(ist, jst, imode, ik, iq)
              end do
            end do
          end do
          write( un, * )
        end do
      end do
    ! JSON output
    else if (sfx == '.json') then
      write( un, '("{",a,": ",a,", ")', advance='no' ) '"bvec"', to_json( bvec )
      write( un, '(a,": ",a,", ")', advance='no' ) '"vkl"', to_json( vkl )
      write( un, '(a,": ",a,", ")', advance='no' ) '"vql"', to_json( vql )
      write( un, '(a,": ",a,", ")', advance='no' ) '"vkql"', to_json( vkql )
      write( un, '(a,": ",a,", ")', advance='no' ) '"elengyk"', to_json( elengyk )
      write( un, '(a,": ",a,", ")', advance='no' ) '"elengykq"', to_json( elengykq )
      write( un, '(a,": ",a,", ")', advance='no' ) '"phengyq"', to_json( phengyq )
      write( un, '(a,": ",a,", ")', advance='no' ) '"g"', to_json( g )
      write( un, '(a,": ",a,"}")' ) '"|g|"', to_json( gavg )
    end if

    close( un )
  end subroutine eph_io_write_ephmat

end module eph_inout
