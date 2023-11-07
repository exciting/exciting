!> I/O utility procedures for phonon related tasks.
module phonons_io_util
  use precision, only: dp
  use asserts, only: assert

  implicit none
  private

  public :: ph_io_q_string, ph_io_qsap_string
  public :: ph_io_write_dynmat_col, ph_io_read_dynmat_col, ph_io_read_dynmat_grid
  public :: ph_io_write_borncharge, ph_io_read_borncharge
  public :: ph_io_write_dielten, ph_io_read_dielten

contains

  !> Get `Q####_####_####` string for a given wave vector.
  function ph_io_q_string( ivq, ngridq ) result( string )
    !> integer grid coordinates of point
    integer, intent(in) :: ivq(3)
    !> grid division
    integer, intent(in) :: ngridq(3)
    !> string
    character(:), allocatable :: string

    integer :: i, j, m(3), n(3)
    character(15) :: tmp

    integer, external :: gcd

    m = 0; n = 0
    do i = 1, 3
      if( ivq(i) == 0 ) cycle
      j = gcd( abs( ivq(i) ), ngridq(i) )
      m(i) = abs( ivq(i) / j )
      n(i) = abs( ngridq(i) / j )
    end do
    write( tmp, '("Q",2i2.2,"_",2i2.2,"_",2i2.2)' ) m(1), n(1), m(2), n(2), m(3), n(3)
    string = tmp
  end function ph_io_q_string

  !> `Q####_####_####_S##_A###_P#` string for a given wave vector, species, atom and polarization direction.
  function ph_io_qsap_string( ivq, ngridq, is, ia, ip ) result( string )
    !> integer grid coordinates of point
    integer, intent(in) :: ivq(3)
    !> grid division
    integer, intent(in) :: ngridq(3)
    !> species index
    integer, intent(in) :: is
    !> atom index
    integer, intent(in) :: ia
    !> polarization direction
    integer, intent(in) :: ip
    !> string
    character(:), allocatable :: string

    character(11) :: tmp

    string = ph_io_q_string( ivq, ngridq )
    write( tmp, '("S",i2.2,"_A",i3.3,"_P",i1.1)' ) is, ia, ip
    string = string // '_' // trim( tmp )
  end function ph_io_qsap_string

  !> write dynamical matrix column to file
  subroutine ph_io_write_dynmat_col( dyn, fxt, success, &
      directory )
    use mod_atoms, only: natmtot, nspecies, natoms, idxas
    !> dynamical matrix column
    complex(dp), intent(in) :: dyn(3, natmtot)
    !> file extension
    character(*), intent(in) :: fxt
    !> `.true.` if writing was successful
    logical, intent(out) :: success
    !> path to directory (default: current directory)
    character(*), optional, intent(in) :: directory

    integer :: stat, un, is, ia, ias, ip
    real(dp) :: a, b
    character(256) :: dirname

    dirname = '.'
    if( present( directory ) ) write( dirname, '(a)' ) trim( directory )
    dirname = trim( dirname )//'/'

    open( newunit=un, file=trim( dirname )//'DYN_'//trim( fxt )//'.OUT', action='write', form='formatted', iostat=stat )
    success = (stat == 0)
    if( .not. success ) return
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        do ip = 1, 3
          a = dble( dyn(ip, ias) )
          b = aimag( dyn(ip, ias) )
          if( abs(a) < 1e-12_dp ) a = 0.0_dp
          if( abs(b) < 1e-12_dp ) b = 0.0_dp
          write( un, '(2g18.10," : is = ",i4,", ia = ",i4,", ip = ",i4)', iostat=stat ) a, b, is, ia, ip
          success = success .and. (stat == 0)
        end do
      end do
    end do
    close( un )
  end subroutine ph_io_write_dynmat_col

  !> read dynamical matrix column from file
  subroutine ph_io_read_dynmat_col( dyn, fxt, success, &
      directory )
    use mod_atoms, only: natmtot, nspecies, natoms, idxas
    !> dynamical matrix col
    complex(dp), intent(out) :: dyn(3, natmtot)
    !> file extension
    character(*), intent(in) :: fxt
    !> `.true.` if reading was successful
    logical, intent(out) :: success
    !> path to directory (default: current directory)
    character(*), optional, intent(in) :: directory

    integer :: stat, un, is, ia, ias, ip, js, ja, jp
    real(dp) :: a, b
    character(256) :: dirname

    dirname = '.'
    if( present( directory ) ) write( dirname, '(a)' ) trim( directory )
    dirname = trim( dirname )//'/'

    open( newunit=un, file=trim( dirname )//'DYN_'//trim( fxt )//'.OUT', action='read', status='old', form='formatted', iostat=stat )
    success = (stat == 0)
    if( .not. success ) return
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        do ip = 1, 3
          js = 0; ja = 0; jp = 0;
          read( un, '(2g18.10,tr1,3(tr7,i4))', iostat=stat ) a, b, js, ja, jp
          success = success .and. (stat == 0)
          if( js /= is .or. ja /= ia .or. jp /= ip ) then
            write( *, * )
            write( *, '("Error (ph_util_read_dynmat): Incompatible file content in file")' )
            write( *, '(a)' ) trim(dirname)//'DYN_'//trim(fxt)//'.OUT'
            write( *, '("expected:       is = ",i4,", ia = ",i4,", ip = ",i4)', iostat=stat ) is, ia, ip
            write( *, '("read from file: is = ",i4,", ia = ",i4,", ip = ",i4)', iostat=stat ) js, ja, jp
            success = .false.
            return
          end if
          dyn(ip, ias) = cmplx( a, b, dp )
        end do
      end do
    end do
    close( un )
  end subroutine ph_io_read_dynmat_col

  !> Read dynamical matrices for a given set of \({\bf q}\)-vectors on a regular grid.
  subroutine ph_io_read_dynmat_grid( ivq, ngridq, nspecies, natoms, dynq, success, directory )
    !> location of \({bf q}\) vectors on integer grid
    integer, intent(in) :: ivq(:, :)
    !> integer divisions of \({\bf q}\)-grid
    integer, intent(in) :: ngridq(3)
    !> number of atomic species
    integer, intent(in) :: nspecies
    !> number of atoms per species
    integer, intent(in) :: natoms(:)
    !> dynamical matrices at \({\bf q}\)
    complex(dp), allocatable, intent(out) :: dynq(:, :, :)
    !> `.true.` if reading was successful
    logical, intent(out) :: success
    !> path to directory (default: current directory)
    character(*), optional, intent(in) :: directory

    integer :: nqpt, natmtot, iq, is, ia, ip, i
    character(256) :: dirname
    character(:), allocatable :: fxt
    complex(dp), allocatable :: dyn_col(:, :)

    dirname = '.'
    if( present( directory ) ) write( dirname, '(a)' ) trim( directory )
    dirname = trim( dirname )//'/'

    nqpt = size( ivq, dim=2 )
    natmtot = sum( natoms(1:nspecies) )

    allocate( dynq(3*natmtot, 3*natmtot, nqpt) )
    allocate( dyn_col(3, natmtot) )
    success = .true.

    do iq = 1, nqpt
      i = 0
      do is = 1, nspecies
        do ia = 1, natoms(is)
          do ip = 1, 3
            i = i + 1
            ! get file extension
            fxt = ph_io_qsap_string( ivq(:, iq), ngridq, is, ia, ip )
            ! read dynamical matrix column from file
            call ph_io_read_dynmat_col( dyn_col, fxt, success, directory=dirname )
            if( .not. success ) return
            dynq(:, i, iq) = reshape( dyn_col, [3*natmtot] )
          end do
        end do
      end do
    end do

    deallocate( dyn_col )
  end subroutine ph_io_read_dynmat_grid

  !> Write Born effective charge tensors \({\bf Z}^\ast_\kappa\) to file.
  subroutine ph_io_write_borncharge( borncharge, fname, success, &
      sumrule_correction )
    use mod_atoms, only: natmtot, nspecies, natoms, idxas, spsymb, atposc
    use mod_lattice, only: ainv
    !> Born effective charge \({\bf Z}^\ast_\kappa\) 
    real(dp), intent(in) :: borncharge(3, 3, natmtot)
    !> file name
    character(*), intent(in) :: fname
    !> `.true.` if writing was successful
    logical, intent(out) :: success
    !> acoustic sum rule correction
    real(dp), optional, intent(in) :: sumrule_correction(3, 3)

    integer :: ierr, un, is, ia, ias, ip
    real(dp) :: vl(3)
    logical :: dosum

    real(dp), allocatable :: borncharge_sum(:,:,:)
    
    dosum = present( sumrule_correction )

    ! try to open file
    open( newunit=un, file=trim( adjustl( fname ) ), action='write', form='formatted', iostat=ierr )
    success = (ierr == 0)
    if( .not. success ) return

    ! write file
    write( un, '("# Born effective charge tensors for all atoms.")' )
    write( un, '("# Rows correspond to E-field direction.")' )
    write( un, '("# Columns correspond to atom displacement direction.")' )
    if( dosum ) &
      write( un, '("# Acoustic sum rule has been imposed.")' )
    write( un, '("#")' )

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        call r3mv( ainv, atposc(:, ia, is), vl )
        write( un, '("# species ",i2," atom ",i3," (",a,i3") : ",3f13.6)' ) is, ia, trim( spsymb(is) ), ia, vl
        do ip = 1, 3
          write( un, '(3f20.10)', iostat=ierr ) borncharge(ip, :, ias)
          success = (ierr == 0)
          if( .not. success ) return
        end do
      end do
    end do

    if( dosum ) then
      write( un, '("# Acoustic sum rule correction (add to each tensor above to get original value)")' )
      do ip = 1, 3
        write( un, '(3f20.10)' ) sumrule_correction(ip, :)
      end do
    end if

    ! close file
    close( un, iostat=ierr )
    success = (ierr == 0)
    if( .not. success ) return
  end subroutine ph_io_write_borncharge

  !> Read Born effective charge tensors \({\bf Z}^\ast_\kappa\) from file.
  subroutine ph_io_read_borncharge( borncharge, fname, success, &
      sumrule_correction )
    use os_utils, only: path_exists
    use mod_atoms, only: natmtot
    !> Born effective charge \({\bf Z}^\ast_\kappa\) 
    real(dp), allocatable, intent(out) :: borncharge(:,:,:)
    !> file name
    character(*), intent(in) :: fname
    !> `.true.` if reading was successful
    logical, intent(out) :: success
    !> acoustic sum rule correction
    real(dp), optional, intent(out) :: sumrule_correction(3, 3)

    integer :: ierr, un, ias, ip
    logical :: dosum
    character(1024) :: line
    
    dosum = present( sumrule_correction )

    success = (path_exists( trim( adjustl( fname ) ), ierr ) .or. ierr /= 0)
    if( .not. success ) return

    if( allocated( borncharge ) ) deallocate( borncharge )
    allocate( borncharge(3, 3, natmtot), source=0.0_dp )
    
    ! try to open file
    open( newunit=un, file=trim( adjustl( fname ) ), action='read', status='old', form='formatted', iostat=ierr )
    success = (ierr == 0)
    if( .not. success ) return
    
    ! read file
    ip = 0; ias = 0
    do while( ierr == 0 )
      read( un, '(a)', iostat=ierr ) line
      line = trim( adjustl( line ) )
      if( line(1:1) == '#' .or. line == '' ) cycle
      if( mod( ip, 3 ) == 0 ) ias = ias + 1
      ip = mod( ip, 3 ) + 1
      read( line, *, iostat=ierr ) borncharge(ip, :, ias)
      if( ias == natmtot .and. ip == 3 ) exit
    end do
    success = (ias == natmtot .and. ip == 3 .and. ierr == 0)
    if( dosum .and. success ) then
      do while( ierr == 0 )
        read( un, '(a)', iostat=ierr ) line
        line = trim( adjustl( line ) )
        if( line(1:1) == '#' .or. line == '' ) cycle
        if( mod( ip, 3 ) == 0 ) ias = ias + 1
        ip = mod( ip, 3 ) + 1
        read( line, *, iostat=ierr ) sumrule_correction(ip, :)
        if( ias == natmtot+1 .and. ip == 3 ) exit
      end do
      success = (ias == natmtot+1 .and. ip == 3 .and. ierr == 0)
    end if

    if( .not. success ) return

    ! close file
    close( un, iostat=ierr )
    success = (ierr == 0)
    if( .not. success ) return
  end subroutine ph_io_read_borncharge

  !> Write high-frequency dielectric tensors \({\bf \epsilon}^\infty\) to file.
  subroutine ph_io_write_dielten( dielten, fname, success )
    !> dielectric tensor \({\bf \epsilon}^\infty\) 
    real(dp), intent(in) :: dielten(3, 3)
    !> file name
    character(*), intent(in) :: fname
    !> `.true.` if writing was successful
    logical, intent(out) :: success

    integer :: ierr, un, ip

    ! try to open file
    open( newunit=un, file=trim( adjustl( fname ) ), action='write', form='formatted', iostat=ierr )
    success = (ierr == 0)
    if( .not. success ) return

    ! write file
    write( un, '("# High frequency dielectric tensor (clamped nuclei).")' )
    write( un, '("#")' )

    do ip = 1, 3
      write( un, '(3f20.10)', iostat=ierr ) dielten(ip, :)
      success = (ierr == 0)
      if( .not. success ) return
    end do

    ! close file
    close( un, iostat=ierr )
    success = (ierr == 0)
    if( .not. success ) return
  end subroutine ph_io_write_dielten

  !> Read high-frequency dielectric tensors \({\bf \epsilon}^\infty\) from file.
  subroutine ph_io_read_dielten( dielten, fname, success )
    use os_utils, only: path_exists
    !> dielectric tensor \({\bf \epsilon}^\infty\) 
    real(dp), intent(out) :: dielten(3, 3)
    !> file name
    character(*), intent(in) :: fname
    !> `.true.` if reading was successful
    logical, intent(out) :: success

    integer :: ierr, un, ip
    character(1024) :: line

    success = (path_exists( trim( adjustl( fname ) ), ierr ) .or. ierr /= 0)
    if( .not. success ) return

    ! try to open file
    open( newunit=un, file=trim( adjustl( fname ) ), action='read', status='old', form='formatted', iostat=ierr )
    success = (ierr == 0)
    if( .not. success ) return
    
    ! read file
    ip = 0
    do while( ierr == 0 )
      read( un, '(a)', iostat=ierr ) line
      line = trim( adjustl( line ) )
      if( line(1:1) == '#' .or. line == '' ) cycle
      ip = ip + 1
      read( line, *, iostat=ierr ) dielten(ip, :)
      if( ip == 3 ) exit
    end do
    success = (ip == 3 .and. ierr == 0)
    if( .not. success ) return

    ! close file
    close( un, iostat=ierr )
    success = (ierr == 0)
    if( .not. success ) return
  end subroutine ph_io_read_dielten

end module phonons_io_util
