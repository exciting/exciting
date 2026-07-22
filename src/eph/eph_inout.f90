!> Collection if I/O routines for electron-phonon calculations.
module eph_inout
  use precision, only: dp
  use modmpi, only: terminate_if_false
#include "asserts.fpp"
  implicit none
  private

  !> allowed integration methods for Fan-Migdal self-energy
  character(len=32), parameter :: eph_else_allowed_methods(2) = [character(32) :: 'smearing', 'kramers-kronig']

  public :: eph_io_write_energies, &
            eph_io_write_ephmat, &
            eph_io_write_el_self_energy, eph_io_read_el_self_energy, &
            eph_io_write_else_sfun, &
            eph_io_write_quasi_particle_energies, &
            eph_io_write_coupling_strength, eph_io_write_integrated_coupling_strength
  
contains

  !================================================================================ 
  ! ENERGIES ON RECIPROCAL SPACE POINTS
  !
  !> Write energies at given reciprocal space points to file.
  !> 
  !> Points can be given as a list of vectors or as a path object.
  subroutine eph_io_write_energies( bvec, i1, e, fname, format, plist, path, xlabel, xvalue )
    use bz_path, only: bz_path_type
    use xjson, only: to_json
    !> reciprocal lattice vectors
    real(dp), intent(in) :: bvec(3, 3)
    !> index of first energy
    integer, intent(in) :: i1
    !> energies
    real(dp), intent(in) :: e(i1:,:)
    !> file name (without suffix!)
    character(len=*), intent(in) :: fname
    !> output format;   
    !> currently supported: `text` (plain text), `json` (JSON dictionary)
    character(len=*), intent(in) :: format
    !> list of vectors representing reciprocal space points
    real(dp), optional, intent(in) :: plist(:,:)
    !> BZ path
    type(bz_path_type), optional, intent(in) :: path
    !> label for optional extra column
    character(len=*), optional, intent(in) :: xlabel
    !> optional extra column
    real(dp), optional, intent(in) :: xvalue(i1:,:)
  
    integer :: un, nst, np, ist, ip, stat
    character(len=8) :: sfx
    character(len=:), allocatable :: xlab

    xlab = ''
    if (present(xlabel)) xlab = trim(xlabel)

    nst = size( e, dim=1 )
    np = size( e, dim=2 )

    ! check input
    CALL_ASSERT( present(plist) .or. present(path), 'Either `plist` or `path` must be present.' )
    CALL_ASSERT( .not. (present(plist) .and. present(path)), 'Not both `plist` and `path` should be present.' )
    if (present(xvalue)) then
      CALL_ASSERT( present(xlabel), '`xlabel` must be present, if `xvalue` is present.' )
    end if

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
        write( un, '("#",a6,2a26)' ) 'band', 'energy (Hartree)', trim(xlab)
        do ip = 1, np
          write( un, '("#",100g26.16e3)' ) plist(:, ip)
          do ist = i1, i1+nst-1
            if (present(xvalue)) then
              write( un, '(i6,2g26.16e3)' ) ist, e(ist, ip), xvalue(ist, ip)
            else
              write( un, '(i6,g26.16e3)' ) ist, e(ist, ip)
            end if
          end do
          write( un, * )
        end do
        close( un )
      ! point list JSON output
      else if (sfx == '.json') then
        open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
        call terminate_if_false( stat == 0, '(eph_io_write_energies) &
          Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
        write( un, '("{",a,": ",a)', advance='no' ) '"bvec"', to_json( bvec )
        write( un, '(", ",a,": ",a)', advance='no' ) '"points"', to_json( plist )
        write( un, '(", ",a,": ",a)', advance='no' ) '"energies"', to_json( e )
        if (present(xvalue)) write( un, '(", ",a,": ",a)', advance='no' ) '"'//trim(xlab)//'"', to_json( xvalue )
        write( un, '("}")' )
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
        write( un, '("#",3a26)' ) 'distance on path (1/bohr)', 'energy (Hartree)', trim(xlab)
        do ist = i1, i1+nst-1
          write( un, '("# band ",i6)' ) ist
          do ip = 1, np
            if (present(xvalue)) then
              write( un, '(3g26.16e3)' ) path%points(ip)%distance, e(ist, ip), xvalue(ist, ip)
            else
              write( un, '(2g26.16e3)' ) path%points(ip)%distance, e(ist, ip)
            end if
          end do
          write( un, * )
        end do
        close( un )
      ! path JSON output
      else if (sfx == '.json') then
        open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat)
        call terminate_if_false( stat == 0, '(eph_io_write_energies) &
          Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
        write( un, '("{",a,": ",a)', advance='no' ) '"path"', path%to_json()
        write( un, '(", ",a,": ",a)', advance='no' ) '"bands"', to_json( e )
        if (present(xvalue)) write( un, '(", ",a,": ",a)', advance='no' ) '"'//trim(xlab)//'"', to_json( xvalue )
        write( un, '("}")' )
        close( un )
      end if
    end if
  end subroutine eph_io_write_energies
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! EPH MATRIX ON RECIPROCAL SPACE POINTS
  !
  !> Write electron-phonon matrix elements at given reciprocal space points to file.
  subroutine eph_io_write_ephmat( bvec, vkl, vkql, vql, i1, elengyk, elengykq, m1, phengyq, g, gavg, fname, format )
    use xjson, only: to_json
    !> reciprocal lattice vectors
    real(dp), intent(in) :: bvec(3, 3)
    !> list of \({\bf k}\)-vectors
    real(dp), intent(in) :: vkl(:,:)
    !> list of \({\bf k}+{\bf q}\)-vectors
    real(dp), intent(in) :: vkql(:,:,:)
    !> list of \({\bf q}\)-vectors
    real(dp), intent(in) :: vql(:,:)
    !> index of first band
    integer, intent(in) :: i1
    !> electron energies \(\epsilon_{n{\bf k}}\)
    real(dp), intent(in) :: elengyk(i1:,:)
    !> electron energies \(\epsilon_{m{\bf k}+{\bf q}}\)
    real(dp), intent(in) :: elengykq(i1:,:,:)
    !> index of first mode
    integer, intent(in) :: m1
    !> phonon frequencies \(\omega_{\nu{\bf q}}\)
    real(dp), intent(in) :: phengyq(m1:,:)
    !> matrix elements \(g_{mn,\nu}({\bf k},{\bf q})\)
    complex(dp), intent(in) :: g(i1:,i1:,m1:,:,:)
    !> absolute value of \(g\) averaged over degenerate states
    real(dp), intent(in) :: gavg(i1:,i1:,m1:,:,:)
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
    CALL_ASSERT( size( elengyk, dim=2 ) == nk, '2nd dimension of `elengyk` must equal 2nd dimension of `vkl`.' )
    CALL_ASSERT( size( elengykq, dim=2 ) == nk, '2nd dimension of `elengykq` must equal 2nd dimension of `vkl`.' )
    CALL_ASSERT( size( elengykq, dim=3 ) == nq, '3rd dimension of `elengykq` must equal 2nd dimension of `vql`.' )
    CALL_ASSERT( size( phengyq, dim=2 ) == nq, '2nd dimension of `phengyq` must equal 2nd dimension of `vql`.' )
    CALL_ASSERT( size( g, dim=1 ) == nstkq, '1st dimension of `g` must equal 1st dimension of `elengykq`.' )
    CALL_ASSERT( size( g, dim=2 ) == nstk, '2nd dimension of `g` must equal 1st dimension of `elengyk`.' )
    CALL_ASSERT( size( g, dim=3 ) == nmode, '3rd dimension of `g` must equal 1st dimension of `phengyq`.' )
    CALL_ASSERT( size( g, dim=4 ) == nk, '4th dimension of `g` must equal 2nd dimension of `vkl`.' )
    CALL_ASSERT( size( g, dim=5 ) == nq, '5th dimension of `g` must equal 2nd dimension of `vql`.' )
    CALL_ASSERT( all( shape(g) == shape(gavg) ), '`g` and `gavg` must have same shape.' )

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
          do imode = m1, m1+nmode-1
            do jst = i1, i1+nstk-1
              do ist = i1, i1+nstkq-1
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
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! ELECTRON SELF-ENERGY
  !
  !> Write electron self-energy \(\Sigma_n({\bf k},\omega,T)\) to binary file.
  subroutine eph_io_write_el_self_energy( file, ik, fst, freqs, temps, el_energy, selfen_fm, selfen_dw, selfen_fm_hilo, selfen_dw_hilo, integration, swidth )
    use block_data_file, only: block_data_file_type
    !> binary file
    type(block_data_file_type), intent(inout) :: file
    !> index of \({\bf k}\)-point
    integer, intent(in) :: ik
    !> first state \(n\) for which self-energy is given
    integer, intent(in) :: fst
    !> frequencies \(\omega\)
    real(dp), intent(in) :: freqs(:)
    !> temperatures \(T\)
    real(dp), intent(in) :: temps(:)
    !> electron energies \(\epsilon_{n{\bf k}}\)
    real(dp), intent(in) :: el_energy(:)
    !> Fan-Migdal self-energy \(\Sigma^{\rm FM}_n({\bf k},\omega,T)\)
    complex(dp), intent(in) :: selfen_fm(:,:,:)
    !> Debye-Waller self-energy \(\Sigma^{\rm DW}_n({\bf k},T)\)
    real(dp), intent(in) :: selfen_dw(:,:)
    !> high and low energy Fan-Migdal self-energy \(\Sigma^{\rm FM,hilo}_n({\bf k},T)\)
    real(dp), intent(in) :: selfen_fm_hilo(:,:)
    !> high and low energy Debye-Waller self-energy \(\Sigma^{\rm DW,hilo}_n({\bf k},T)\)
    real(dp), intent(in) :: selfen_dw_hilo(:,:)
    !> integration method used to compute self-energy
    character(len=*), intent(in) :: integration
    !> imaginary broadening \(\eta\) used to compute self-energy
    real(dp), intent(in) :: swidth
  
    integer :: nfreq, ntemp, nst, itemp, icode

    integer, allocatable :: shp(:)
    complex(dp), allocatable :: record(:,:,:)

    ! set matrix sizes
    nfreq = size( freqs )
    nst = size( el_energy )
    ntemp = size( temps )

    ! read shape of data block in file
    shp = file%get_block_shape()

    ! check input
    CALL_ASSERT( size(shp) == 3, 'File with data blocks on rank 3 expected.' )
    CALL_ASSERT( shp(1) == nfreq+2, '1st dimension of data block in file must equal number of frequencies + 2.' )
    CALL_ASSERT( shp(2) == nst+1, '2nd dimension of data block in file must equal number of bands + 1.' )
    CALL_ASSERT( shp(3) == ntemp, '3rd dimension of data block in file must equal number of temperatures.' )
    CALL_ASSERT( ik > 0, 'k-point index `ik` must be positive.' )
    CALL_ASSERT( size( selfen_fm, dim=1 ) == nfreq, '1st dimension of `selfen_fm` must equal number of frequencies.' )
    CALL_ASSERT( size( selfen_fm, dim=2 ) == nst, '2nd dimension of `selfen_fm` must equal number of bands.' )
    CALL_ASSERT( size( selfen_fm, dim=3 ) == ntemp, '3rd dimension of `selfen_fm` must equal number of temperatures.' )
    CALL_ASSERT( size( selfen_dw, dim=1 ) == nst, '1st dimension of `selfen_dw` must equal number of bands.' )
    CALL_ASSERT( size( selfen_dw, dim=2 ) == ntemp, '2nd dimension of `selfen_dw` must equal number of temperatures.' )
    CALL_ASSERT( size( selfen_fm_hilo, dim=1 ) == nst, '1st dimension of `selfen_fm_hilo` must equal number of bands.' )
    CALL_ASSERT( size( selfen_fm_hilo, dim=2 ) == ntemp, '2nd dimension of `selfen_fm_hilo` must equal number of temperatures.' )
    CALL_ASSERT( size( selfen_dw_hilo, dim=1 ) == nst, '1st dimension of `selfen_dw_hilo` must equal number of bands.' )
    CALL_ASSERT( size( selfen_dw_hilo, dim=2 ) == ntemp, '2nd dimension of `selfen_dw_hilo` must equal number of temperatures.' )

    allocate( record(shp(1), shp(2), shp(3)) )

    ! get integer code of used integration method
    do icode = size( eph_else_allowed_methods ), 1, -1
      if (trim( adjustl( integration ) ) == trim( adjustl( eph_else_allowed_methods(icode) ) )) exit
    end do

    ! layout
    ! FM          | freqs
    ! DW, FM HILO | temp, fst
    ! epsilon     | method, swidth
    do itemp = 1, ntemp
      record(1:nfreq, 1:nst, itemp) = selfen_fm(:, :, itemp)
      record(nfreq+1, 1:nst, itemp) = cmplx( selfen_dw(:, itemp), el_energy, dp )
      record(nfreq+2, 1:nst, itemp) = cmplx( selfen_fm_hilo(:, itemp), selfen_dw_hilo(:, itemp), dp )
      record(1:nfreq, nst+1, itemp) = cmplx( freqs, 0, dp )
      record(nfreq+1, nst+1, itemp) = cmplx( temps(itemp), fst, dp )
      record(nfreq+2, nst+1, itemp) = cmplx( icode, swidth, dp )
    end do

    call file%write( ik, record )

    deallocate( record )
  end subroutine eph_io_write_el_self_energy

  !> Read electron self-energy \(\Sigma_n({\bf k},\omega,T)\) from binary file.
  subroutine eph_io_read_el_self_energy( file, ik, fst, freqs, temps, el_energy, selfen_fm, selfen_dw, selfen_fm_hilo, selfen_dw_hilo, integration, swidth )
    use block_data_file, only: block_data_file_type
    !> binary file
    type(block_data_file_type), intent(inout) :: file
    !> index of \({\bf k}\)-point
    integer, intent(in) :: ik
    !> first state \(n\) for which self-energy is given
    integer, intent(out) :: fst
    !> frequencies \(\omega\)
    real(dp), allocatable, intent(out) :: freqs(:)
    !> temperatures \(T\)
    real(dp), allocatable, intent(out) :: temps(:)
    !> electron energies \(\epsilon_{n{\bf k}}\)
    real(dp), allocatable, intent(out) :: el_energy(:)
    !> Fan-Migdal self-energy \(\Sigma^{\rm FM}_n({\bf k},\omega,T)\)
    complex(dp), allocatable, intent(out) :: selfen_fm(:,:,:)
    !> Debye-Waller self-energy \(\Sigma^{\rm DW}_n({\bf k},T)\)
    real(dp), allocatable, intent(out) :: selfen_dw(:,:)
    !> high and low energy Fan-Migdal self-energy \(\Sigma^{\rm FM,hilo}_n({\bf k},T)\)
    real(dp), allocatable, intent(out) :: selfen_fm_hilo(:,:)
    !> high and low energy Debye-Waller self-energy \(\Sigma^{\rm Dw,hilo}_n({\bf k},T)\)
    real(dp), allocatable, intent(out) :: selfen_dw_hilo(:,:)
    !> integration method used to compute self-energy
    character(len=:), allocatable, intent(out) :: integration
    !> imaginary broadening \(\eta\) used to compute self-energy
    real(dp), intent(out) :: swidth
  
    integer :: nfreq, ntemp, nst, lst, itemp, icode

    integer, allocatable :: shp(:)
    complex(dp), allocatable :: record(:,:,:)

    ! layout
    ! FM          | freqs
    ! DW, epsilon | temp, fst
    ! HILO        | method, swidth

    ! read shape of data block in file
    shp = file%get_block_shape()

    ! check input
    CALL_ASSERT( size(shp) == 3, 'File with data blocks on rank 3 expected.' )
    CALL_ASSERT( ik > 0, 'k-point index `ik` must be positive.' )

    ! set matrix sizes
    nfreq = shp(1) - 2
    nst = shp(2) - 1
    ntemp = shp(3)

    allocate( record(shp(1), shp(2), shp(3)) )
    call file%read( ik, record )

    fst          = nint(record(nfreq+1, nst+1, 1)%im)
    lst          = fst + nst - 1
    icode        = nint(record(nfreq+2, nst+1, 1)%re)
    swidth       = record(nfreq+2, nst+1, 1)%im

    allocate( freqs(nfreq), temps(ntemp), el_energy(fst:lst) ) 
    allocate( selfen_fm(nfreq, fst:lst, ntemp), selfen_dw(fst:lst, ntemp), selfen_fm_hilo(fst:lst, ntemp), selfen_dw_hilo(fst:lst, ntemp) )

    el_energy    = record(nfreq+1, 1:nst, 1)%im
    freqs        = record(1:nfreq, nst+1, 1)%re
    do itemp = 1, ntemp
      selfen_fm(:, :, itemp) = record(1:nfreq, 1:nst, itemp)
      selfen_dw(:, itemp)    = record(nfreq+1, 1:nst, itemp)%re
      selfen_fm_hilo(:, itemp)  = record(nfreq+2, 1:nst, itemp)%re
      selfen_dw_hilo(:, itemp)  = record(nfreq+2, 1:nst, itemp)%im
      temps(itemp) = record(nfreq+1, nst+1, itemp)%re
    end do
    if (icode > 0 .and. icode <= size(eph_else_allowed_methods)) then
      integration = trim( adjustl( eph_else_allowed_methods(icode) ) )
    else
      integration = 'unknown'
    end if

    deallocate( record )
  end subroutine eph_io_read_el_self_energy
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! ELECTRON SELF-ENERGY & SPECTRAL FUNCTION
  !> Write electron self-energy \(\Sigma_{nn}({\bf k},\omega)\) and the corresponding
  !> spectral function \(A_n({\bf k},\omega)\) to file.
  subroutine eph_io_write_else_sfun( freqs, selfen, sfun, format, fname, el_energy, selfen_dw )
    use xjson, only: to_json
    !> frequencies \(\omega\) for each band \(n\)
    real(dp), intent(in) :: freqs(:,:)
    !> total self-energy \(\Sigma_{nn}({\bf k},\omega)\) for each band \(n\)
    complex(dp), intent(in) :: selfen(:,:)
    !> spectral function \(A_n({\bf k},\omega)\) for each band \(n\)
    real(dp), intent(in) :: sfun(:,:)
    !> output format;   
    !> currently supported: `text` (plain text), `json` (JSON dictionary)
    character(len=*), intent(in) :: format
    !> file name
    character(len=*), intent(in) :: fname
    !> electron energies \(\epsilon_{n{\bf k}}\) for each band \(n\)
    real(dp), optional, intent(in) :: el_energy(:)
    !> Debye-Waller self-energy contribution \(\Sigma^{\rm DW}_{nn}({\bf k})\) for each band \(n\)
    real(dp), optional, intent(in) :: selfen_dw(:)

    integer :: un, stat, nfreq, nst, ifreq, ist
    character(len=8) :: sfx

    nfreq = size( freqs, dim=1 )
    nst = size( selfen, dim=2 )

    sfx = ''
    select case (trim( adjustl( format ) ))
      case ('text')
        sfx = '.dat'
      case ('json')
        sfx = '.json'
      case default
        call terminate_if_false( .false., '(eph_io_write_else_sfun) &
          Unsupported format `'//trim( adjustl( format ) )//'`.' )
    end select

    open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
    call terminate_if_false( stat == 0, '(eph_io_write_else_sfun) &
      Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )

    ! text output
    if (sfx == '.dat') then
      write( un, '("#",a25,3a26)' ) 'frequency', 'Re(self-energy)', 'Im(self-energy)', 'spectral function'
      do ist = 1, nst
        if (present(el_energy)) write( un, '("# electron energy         :",g26.16e3)' ) el_energy(ist)
        if (present(selfen_dw)) write( un, '("# Debye-Waller self-energy:",g26.16e3)' ) selfen_dw(ist)
        do ifreq = 1, nfreq
          write( un, '(4g26.16e3)' ) freqs(ifreq, min(ist, size(freqs, dim=2))), selfen(ifreq, ist), sfun(ifreq, ist)
        end do
        write( un, * )
        write( un, * )
      end do
    ! JSON output
    else if (sfx == '.json') then
      write( un, '("{",a,": ",a,", ")', advance='no' ) '"frequencies"', to_json( freqs )
      write( un, '(a,": ",a,", ")', advance='no' ) '"self-energy"', to_json( selfen )
      write( un, '(a,": ",a)', advance='no' ) '"spectral function"', to_json( sfun )
      if (present(el_energy)) write( un, '(", ",a,": ",a)', advance='no' ) '"electron energy"', to_json( el_energy )
      if (present(selfen_dw)) write( un, '(", ",a,": ",a)', advance='no' ) '"DW self-energy"', to_json( selfen_dw )
      write( un, '("}")' )
    end if

    close( un )
  end subroutine eph_io_write_else_sfun
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! QUASI PARTICLE ENERGIES
  !
  !> Write renormalized quasi-particle energies and different self-energy contributions to file. 
  !> Output is similar to EVALQP.DAT from GW module.
  subroutine eph_io_write_quasi_particle_energies( kset, i1, eSP, eQP, sFM, sDW, sFM_hilo, sDW_hilo, Z, fname, format )
    use mod_kpointset, only: k_set
    use xjson, only: to_json
    !> k-point set
    type(k_set), intent(in) :: kset
    !> index of first energy
    integer, intent(in) :: i1
    !> single particle energies \(\epsilon_{n{\bf k}}\)
    real(dp), intent(in) :: eSP(i1:,:)
    !> quasi particle energies \(\epsilon^{\rm QP}_{n{\bf k}}\)
    complex(dp), intent(in) :: eQP(i1:,:)
    !> Fan-Migdal self-energy \(\Sigma^{\rm FM}_{nn}({\bf k})\)
    complex(dp), intent(in) :: sFM(i1:,:)
    !> Debye-Waller self-energy \(\Sigma^{\rm FM}_{nn}({\bf k})\)
    real(dp), intent(in) :: sDW(i1:,:)
    !> high and low energy Fan-Migdal self-energy \(\Sigma^{\rm FM, hilo}_{nn}({\bf k})\)
    real(dp), intent(in) :: sFM_hilo(i1:,:)
    !> high and low energy Debye-Waller self-energy \(\Sigma^{\rm DW, hilo}_{nn}({\bf k})\)
    real(dp), intent(in) :: sDW_hilo(i1:,:)
    !> quasi particle strength \(Z_{n{\bf k}\)
    complex(dp), intent(in) :: Z(i1:,:)
    !> file name (without suffix!)
    character(len=*), intent(in) :: fname
    !> output format;   
    !> currently supported: `text` (plain text), `json` (JSON dictionary)
    character(len=*), intent(in) :: format
    
    integer :: un, stat, nst, ik, ist
    character(len=8) :: sfx

    nst = size( eSP, dim=1 )

    sfx = ''
    select case (trim( adjustl( format ) ))
      case ('text')
        sfx = '.dat'
      case ('json')
        sfx = '.json'
      case default
        call terminate_if_false( .false., '(eph_io_write_quasi_particle_energies) &
          Unsupported format `'//trim( adjustl( format ) )//'`.' )
    end select

    ! text output
    if (sfx == '.dat') then
      open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
      call terminate_if_false( stat == 0, '(eph_io_write_quasi_particle_energies) &
        Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
      do ik = 1, kset%nkpt
        write( un, '("k-point #",i6,":",4f12.6)' ) ik, kset%vkl(:, ik), kset%wkpt(ik)
        write( un, '(a6,11(a16,x))' ) 'state', 'E_KS[Ha]', 'Re(E_EPH)[Ha]', 'Im(E_EPH)[Ha]', 'Re(S_FM)[Ha]', 'Im(S_FM)[Ha]', 'S_DW[Ha]', 'S_FM_hilo[Ha]', 'S_DW_hilo[Ha]', 'DE_EPH[Ha]', 'Re(Znk)', 'Im(Znk)'
        do ist = i1, i1+nst-1
          write( un, '(i4,2x,11(f16.8,x))' ) ist, eSP(ist, ik), eQP(ist, ik), sFM(ist, ik), sDW(ist, ik), sFM_hilo(ist, ik), sDW_hilo(ist, ik), eQP(ist, ik)%re-eSP(ist, ik), Z(ist, ik)
        end do
        write( un, * )
      end do
      close( un )
    ! JSON output
    else if (sfx == '.json') then
      open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
      call terminate_if_false( stat == 0, '(eph_io_write_quasi_particle_energies) &
        Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
      write( un, '("{",a,": ",a)', advance='no' ) '"bvec"', to_json( kset%bvec )
      write( un, '(", ",a,": ",a)', advance='no' ) '"vkl"', to_json( kset%vkl(:, 1:kset%nkpt) )
      write( un, '(", ",a,": ",a)', advance='no' ) '"SP energies"', to_json( eSP )
      write( un, '(", ",a,": ",a)', advance='no' ) '"QP energies"', to_json( eQP )
      write( un, '(", ",a,": ",a)', advance='no' ) '"Sigma FM"', to_json( sFM )
      write( un, '(", ",a,": ",a)', advance='no' ) '"Sigma DW"', to_json( sDW )
      write( un, '(", ",a,": ",a)', advance='no' ) '"Sigma FM (hilo)"', to_json( sFM_hilo )
      write( un, '(", ",a,": ",a)', advance='no' ) '"Sigma DW (hilo)"', to_json( sDW_hilo )
      write( un, '(", ",a,": ",a)', advance='no' ) '"QP strength"', to_json( Z )
      write( un, '("}")' )
      close( un )
    end if
  end subroutine eph_io_write_quasi_particle_energies
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! EPH COUPLING STRENGTH
  !
  !> Write EPH coupling strength at given reciprocal space points to file.
  !> 
  !> Points can be given as a list of vectors or as a path object.
  subroutine eph_io_write_coupling_strength( bvec, i1, e, lambda, fname, format, plist, path )
    use bz_path, only: bz_path_type
    use xjson, only: to_json
    !> reciprocal lattice vectors
    real(dp), intent(in) :: bvec(3, 3)
    !> index of first energy
    integer, intent(in) :: i1
    !> energies
    real(dp), intent(in) :: e(i1:,:)
    !> coupling strength
    real(dp), intent(in) :: lambda(i1:,:,:)
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
    CALL_ASSERT( .not. (present(plist) .and. present(path)), 'Not both `plist` and `path` should be present.' )
    CALL_ASSERT( size( lambda, dim=1 ) == nst, 'First dimension of arrays `e` and `lambda` must match.' )
    CALL_ASSERT( size( lambda, dim=2 ) == np, 'Second dimension of arrays `e` and `lambda` must match.' )

    sfx = ''
    select case (trim( adjustl( format ) ))
      case ('text')
        sfx = '.dat'
      case ('json')
        sfx = '.json'
      case default
        call terminate_if_false( .false., '(eph_io_write_coupling_strength) &
          Unsupported format `'//trim( adjustl( format ) )//'`.' )
    end select

    if (present(plist)) then
      call terminate_if_false( size( plist, dim=2 ) == np, '(eph_io_write_coupling_strength) &
        Energies `e` and point list `plist` must contain same number of points.' )
      ! point list text output
      if (sfx == '.dat') then
        open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
        call terminate_if_false( stat == 0, '(eph_io_write_coupling_strength) &
          Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
        if (size(lambda, dim=3) == 1) then
          write( un, '("#",a6,2a26)' ) 'band', 'energy (Hartree)', 'coupling (quasielastic)'
        else
          write( un, '("#",a6,3a26)' ) 'band', 'energy (Hartree)', 'coupling (absorption)', 'coupling (emission)'
        end if
        do ip = 1, np
          write( un, '("#",100g26.16e3)' ) plist(:, ip)
          do ist = i1, i1+nst-1
            write( un, '(i6,100g26.16e3)' ) ist, e(ist, ip), lambda(ist, ip, :)
          end do
          write( un, * )
        end do
        close( un )
      ! point list JSON output
      else if (sfx == '.json') then
        open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
        call terminate_if_false( stat == 0, '(eph_io_write_coupling_strength) &
          Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
        write( un, '("{")' )
        write( un, '(a,": ",a,", ")' ) '"bvec"', to_json( bvec )
        write( un, '(a,": ",a,", ")' ) '"points"', to_json( plist )
        write( un, '(a,": ",a,", ")' ) '"energies"', to_json( e )
        if (size(lambda, dim=3) == 1) then
          write( un, '(a,": ",a)' ) '"coupling (quasielastic)"', to_json( lambda(:, :, 1) )
        else
          write( un, '(a,": ",a,", ")' ) '"coupling (absorption)"', to_json( lambda(:, :, 1) )
          write( un, '(a,": ",a)' ) '"coupling (emission)"', to_json( lambda(:, :, 2) )
        end if
        write( un, '("}")' )
        close( un )
      end if
    else if (present(path)) then
      call terminate_if_false( path%num_points == np, '(eph_io_write_coupling_strength) &
        Energies `e` and path `path` must contain same number of points.' )
      ! path text output
      if (sfx == '.dat') then
        open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
        call terminate_if_false( stat == 0, '(eph_io_write_coupling_strength) &
          Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
        if (size(lambda, dim=3) == 1) then
          write( un, '("#",3a26)' ) 'distance on path (1/bohr)', 'energy (Hartree)', 'coupling (quasielastic)'
        else
          write( un, '("#",4a26)' ) 'distance on path (1/bohr)', 'energy (Hartree)', 'coupling (absorption)', 'coupling (emission)'
        end if
        do ist = i1, i1+nst-1
          write( un, '("# band ",i6)' ) ist
          do ip = 1, np
            write( un, '(100g26.16e3)' ) path%points(ip)%distance, e(ist, ip), lambda(ist, ip, :)
          end do
          write( un, * )
        end do
        close( un )
      ! path JSON output
      else if (sfx == '.json') then
        open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat)
        call terminate_if_false( stat == 0, '(eph_io_write_coupling_strength) &
          Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
        write( un, '("{")' )
        write( un, '(a,": ",a,", ")' ) '"path"', path%to_json()
        write( un, '(a,": ",a,", ")' ) '"bands"', to_json( e )
        if (size(lambda, dim=3) == 1) then
          write( un, '(a,": ",a)' ) '"coupling strength (quasielastic)"', to_json( lambda(:, :, 1) )
        else
          write( un, '(a,": ",a,", ")' ) '"coupling strength (absorption)"', to_json( lambda(:, :, 1) )
          write( un, '(a,": ",a)' ) '"coupling strength (emission)"', to_json( lambda(:, :, 2) )
        end if
        write( un, '("}")' )
        close( un )
      end if
    end if
  end subroutine eph_io_write_coupling_strength

  !> Write integrated EPH coupling strength at given set of frequencies to file.
  subroutine eph_io_write_integrated_coupling_strength( freqs, dos, lambda, fname, format, cumulative_lambda )
    use xjson, only: to_json
    !> frequencies
    real(dp), intent(in) :: freqs(:)
    !> density of states
    real(dp), intent(in) :: dos(:)
    !> integrated coupling strength
    real(dp), intent(in) :: lambda(:,:)
    !> file name (without suffix!)
    character(len=*), intent(in) :: fname
    !> output format;   
    !> currently supported: `text` (plain text), `json` (JSON dictionary)
    character(len=*), intent(in) :: format
    !> cumulative coupling strength (phonon coupling only)
    real(dp), optional, intent(in) :: cumulative_lambda(:,:)
  
    integer :: un, nfreq, ifreq, stat
    character(len=8) :: sfx

    nfreq = size( freqs )

    ! check input
    CALL_ASSERT( size( lambda, dim=1 ) == nfreq, 'First dimension of arrays `freqs` and `lambda` must match.' )
    CALL_ASSERT( size( dos ) == nfreq, 'Size of arrays `freqs` and `dos` must match.' )
    if (present(cumulative_lambda)) then
      CALL_ASSERT( size( cumulative_lambda, dim=1 ) == nfreq, 'First dimension of arrays `freqs` and `cumulative_lambda` must match.' )
    end if

    sfx = ''
    select case (trim( adjustl( format ) ))
      case ('text')
        sfx = '.dat'
      case ('json')
        sfx = '.json'
      case default
        call terminate_if_false( .false., '(eph_io_write_integrated_coupling_strength) &
          Unsupported format `'//trim( adjustl( format ) )//'`.' )
    end select

    ! text output
    if (sfx == '.dat') then
      open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
      call terminate_if_false( stat == 0, '(eph_io_write_integrated_coupling_strength) &
        Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
      if (size(lambda, dim=2) == 1) then
        write( un, '("#",3a26)', advance='no' ) 'frequency (Hartree)', 'DOS (states/cell/Hartree)', 'coupling (quasielastic)'
      else
        write( un, '("#",4a26)', advance='no' ) 'frequency (Hartree)', 'DOS (states/cell/Hartree)', 'coupling (absorption)', 'coupling (emission)'
      end if
      if (present(cumulative_lambda)) then
        if (size(cumulative_lambda, dim=2) == 1) then
          write( un, '(a26)' ) 'cumulative (quasielastic)'
        else
          write( un, '(2a26)' ) 'cumulative (absorption)', 'cumulative (emission)'
        end if
      else
        write( un, * )
      end if
      do ifreq = 1, nfreq
        if (present(cumulative_lambda)) then
          write( un, '(100g26.16e3)' ) freqs(ifreq), dos(ifreq), lambda(ifreq, :), cumulative_lambda(ifreq, :)
        else
          write( un, '(100g26.16e3)' ) freqs(ifreq), dos(ifreq), lambda(ifreq, :)
        end if
      end do
      close( un )
    ! JSON output
    else if (sfx == '.json') then
      open( newunit=un, file=trim( adjustl( fname ) )//trim( sfx ), action='write', form='formatted', iostat=stat )
      call terminate_if_false( stat == 0, '(eph_io_write_integrated_coupling_strength) &
        Failed to open file `'//trim( adjustl( fname ) )//trim( sfx )//'`.' )
      write( un, '("{")' )
      write( un, '(a,": ",a,", ")' ) '"frequencies"', to_json( freqs )
      write( un, '(a,": ",a,", ")' ) '"dos"', to_json( dos )
      if (size(lambda, dim=2) == 1) then
        write( un, '(a,": ",a)' ) '"coupling (quasielastic)"', to_json( lambda(:, 1) )
      else
        write( un, '(a,": ",a,", ")' ) '"coupling (absorption)"', to_json( lambda(:, 1) )
        write( un, '(a,": ",a)', advance='no' ) '"coupling (emission)"', to_json( lambda(:, 2) )
      end if
      if (present(cumulative_lambda)) then
        write( un, '(", ")' )
        if (size(cumulative_lambda, dim=2) == 1) then
          write( un, '(a,": ",a)' ) '"cumulative coupling (quasielastic)"', to_json( cumulative_lambda(:, 1) )
        else
          write( un, '(a,": ",a,", ")' ) '"cumulative coupling (absorption)"', to_json( cumulative_lambda(:, 1) )
          write( un, '(a,": ",a)' ) '"cumulative coupling (emission)"', to_json( cumulative_lambda(:, 2) )
        end if
      else
        write( un, * )
      end if
      write( un, '("}")' )
      close( un )
    end if
  end subroutine eph_io_write_integrated_coupling_strength
  !-------------------------------------------------------------------------------- 

end module eph_inout
