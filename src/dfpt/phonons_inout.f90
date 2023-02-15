!> This module contains i/o functionalities for DFPT phonons calculations.
!> See also module [[dfpt_inout(module)]].
module phonons_inout
  use, intrinsic :: iso_fortran_env, only: stdout => output_unit

  use dfpt_variables
  use dfpt_inout
  use phonons_variables

  use precision, only: dp
  use modmpi
  use FoX_wxml

  implicit none
  private

  character(:), allocatable, public :: ph_io_qi_dirname

  public :: ph_io_irrep_fxt, ph_io_canon_fxt
  public :: ph_io_genqidir
  public :: ph_io_find_done_parts
  public :: ph_io_info_init, ph_io_info_finit, &
            ph_io_info_symmetries_and_parallelization, ph_io_info_irrep_patterns
  public :: ph_io_write_dforce_const, ph_io_read_dforce_const
  public :: ph_io_print_calc_info

  contains

    !> get the file extension for a given \({\bf q}\) point, irrep, and irrep member
    subroutine ph_io_irrep_fxt( iq, iirrep, dirrep, fxt )
      !> index of \({\bf q}\) point
      integer, intent(in) :: iq
      !> index of the irrep
      integer, intent(in) :: iirrep
      !> irrep member
      integer, intent(in) :: dirrep
      !> file extension
      character(:), allocatable, intent(out) :: fxt

      integer :: i, j, ivq(3), m(3), n(3)
      character(64) :: tmp

      integer, external :: gcd

      ivq = ph_qset%ivk(:, iq)
      m = 0; n = 0
      do i = 1, 3
        if( ivq(i) == 0 ) cycle
        j = gcd( abs( ivq(i) ), ph_qset%ngridk(i) )
        m(i) = abs( ivq(i) / j )
        n(i) = abs( ph_qset%ngridk(i) / j )
      end do
      write( tmp, '("Q",2i2.2,"_",2i2.2,"_",2i2.2,"_I",i3.3)' ) m(1), n(1), m(2), n(2), m(3), n(3), iirrep
      if( dirrep > 0 ) write( tmp, '(a,"_",i1.1)' ) trim( tmp ), dirrep
      fxt = trim( tmp )
    end subroutine ph_io_irrep_fxt

    !> get file extension for a given \({\bf q}\) point, atom and Cartesian direction
    subroutine ph_io_canon_fxt( iq, is, ia, ip, fxt)
      !> index of \({\bf q}\) point
      integer, intent(in) :: iq
      !> species index
      integer, intent(in) :: is
      !> atom index
      integer, intent(in) :: ia
      !> Cartesian direction
      integer, intent(in) :: ip
      !> file extension
      character(:), allocatable, intent(out) :: fxt

      integer :: i, j, ivq(3), m(3), n(3)
      character(64) :: tmp

      integer, external :: gcd

      ivq = ph_qset%ivk(:, iq)
      m = 0; n = 0
      do i = 1, 3
        if( ivq(i) == 0 ) cycle
        j = gcd( abs( ivq(i) ), ph_qset%ngridk(i) )
        m(i) = abs( ivq(i) / j )
        n(i) = abs( ph_qset%ngridk(i) / j )
      end do
      write( tmp, '("Q",2i2.2,"_",2i2.2,"_",2i2.2,"_S",i2.2,"_A",i3.3,"_P",i1.1)' ) m(1), n(1), m(2), n(2), m(3), n(3), is, ia, ip
      fxt = trim( tmp )
    end subroutine ph_io_canon_fxt

    !> generate directory for independent \(({\bf q},I)\) part
    subroutine ph_io_genqidir(iq, iirrep, qidirname, &
        comm )
      use os_utils, only: make_directory
      !> index of \({\bf q}\) point
      integer, intent(in) :: iq
      !> index of the irrep
      integer, intent(in) :: iirrep
      !> name of directory
      character(:), allocatable, intent(out) :: qidirname
      !> MPI communicator (default: `mpiglobal`)
      type(mpiinfo), optional, intent(in) :: comm

      integer :: i
      type(mpiinfo) :: mpi

      mpi = mpiglobal
      if( present( comm ) ) mpi = comm

      call ph_io_irrep_fxt( iq, iirrep, 0, qidirname )
      qidirname = './'//qidirname//'/'
      i = make_directory( qidirname, comm=mpi ) 

      if( i /= 0 ) then
        if( mpi%rank == 0 ) then
          write( *, * )
          write( *, '("Error (ph_io_genqidir): Error in creating (q,I) subdirectory.")' )
          write( *, '(" Directory:  ",a)' ) trim( qidirname )
        end if
        call terminate
      end if
    end subroutine ph_io_genqidir

    !> Find which parts have already been done, i.e., for which parts
    !> the dynamical matrix has been computed.
    subroutine ph_io_find_done_parts( qset, basis, qi_done )
      use phonons_symmetry, only: irrep_basis
      use mod_kpointset, only: k_set
      use mod_atoms, only: natmtot
      use os_utils, only: path_exists
      !> set of \({\bf q}\) points
      type(k_set), intent(in) :: qset
      !> irrep basis per \({\bf q}\) point
      type(irrep_basis), intent(in) :: basis(qset%nkpt)
      !> `.true.` if \(({\bf q},I)\) pair has been precomputed
      logical, intent(out) :: qi_done(qset%nkpt, 3*natmtot)

      integer :: iq, iirrep, dirrep, ierr
      logical :: done
      character(:), allocatable :: qidirname, fxt

      qi_done = .false.

      do iq = 1, qset%nkpt
        do iirrep = 1, basis(iq)%nirrep
          call ph_io_irrep_fxt( iq, iirrep, 0, qidirname )
          done = .true.
          do dirrep = 1, basis(iq)%irreps(iirrep)%dim
            call ph_io_irrep_fxt( iq, iirrep, dirrep, fxt )
            done = done .and. path_exists( trim( qidirname )//'/DYN_'//trim( fxt )//'.OUT', ierr )
          end do
          qi_done(iq, iirrep) = done
        end do
      end do
    end subroutine ph_io_find_done_parts

    !> This subroutine opens and initializes the general info output file (human readable and xml)
    !> for a given\({\bf q}\) point and irrep.
    subroutine ph_io_info_init( iq, iirrep, &
        directory )
      !> index of \({\bf q}\) point
      integer, intent(in) :: iq
      !> index of the irrep
      integer, intent(in) :: iirrep
      !> path to directory (default: current directory)
      character(*), optional, intent(in) :: directory

      character(256) :: dirname, string
      character(:), allocatable :: fxt

      dirname = '.'
      if( present( directory ) ) write( dirname, '(a)' ) trim( directory )
      dirname = trim(dirname)//'/'

      ! get file extension
      call ph_io_irrep_fxt( iq, iirrep, 0, fxt )

      call dfpt_io_info_init( file_extension=fxt, directory=dirname )

      write( info_out_unit, * )
      write( string, '(i6,3f10.6,i6)' ) iq, ph_qset%vkl(:, iq), iirrep
      write( info_out_unit, '("Running calculation for q-point ",a," (",a,") ",a,"and irreducible representation ",a,".")') &
        trim( adjustl( string(1:6) ) ), trim( adjustl( string(7:36) ) ), new_line( 'a' ), trim( adjustl( string(37:) ) )
      write( string, '(i6)' ) iq
      call xml_AddAttribute( info_out_xml, 'q_idx', trim( adjustl( string ) ) )
      write( string, '(3g25.16)' ) ph_qset%vkl(:, iq)
      call xml_AddAttribute( info_out_xml, 'vql', trim( adjustl( string ) ) )
      write( string, '(i6)' ) iirrep
      call xml_AddAttribute( info_out_xml, 'irrep_idx', trim( adjustl( string ) ) )
    end subroutine ph_io_info_init

    !> This subroutine finalizes and closes the general info output file (human readable and xml).
    subroutine ph_io_info_finit
      call dfpt_io_info_finit
    end subroutine ph_io_info_finit

    !> This subroutine writes symmetry information and parallelization settings
    !> to the general info output file (human readable and xml).
    subroutine ph_io_info_symmetries_and_parallelization( qset, basis, parts, ipart, parts_all )
      use phonons_symmetry
      use phonons_parallelization, only: ph_part
      use mod_kpointset, only: k_set
      !> set of \({\bf q}\) points
      type(k_set), intent(in) :: qset
      !> irrep basis per \({\bf q}\) point
      type(irrep_basis), intent(in) :: basis(qset%nkpt)
      !> array of (remaining) independent parts
      type(ph_part), intent(in) :: parts(:)
      !> index of current part
      integer, intent(in) :: ipart
      !> array of all independent parts
      type(ph_part), intent(in) :: parts_all(:)

      character(256) :: string
      character(:), allocatable :: long
      integer :: iq, ip, ii, irank, totload, remload, maxload, nidle

      ! write symmetries
      call printbox( info_out_unit, '*', 'symmetries and parallelization settings' )
      write( info_out_unit, * )
      if( ph_canonical ) then
        write( info_out_unit, '("Canonical displacement patterns are used.")' )
        write( info_out_unit, '("k-points are not reduced by symmetries.")' )
        write( info_out_unit, '("Response quantities are not symmetrized.")' )
        string = 'canonical'
      else
        write( info_out_unit, '("Irreducible representations are used.")' )
        if( ph_kset%isreduced ) then
          write( info_out_unit, '("k-points are reduced by symmetries.")' )
        else
          write( info_out_unit, '("k-points are not reduced by symmetries.")' )
        end if
        write( info_out_unit, '("Response quantities are symmetrized.")' )
        string = 'irreps'
      end if
      write( info_out_unit, * )
      write( info_out_unit, '("q-points and symmetries:")' )
      write( info_out_unit, * )
      do ip = 1, size( parts_all )
        if( parts(ipart)%iq == parts_all(ip)%iq .and. parts(ipart)%iirrep == parts_all(ip)%iirrep ) exit
      end do
      call write_qi_parts( info_out_unit, qset, basis, parts_all, ip )
      call xml_NewElement( info_out_xml, 'symmetries' )
      call xml_AddAttribute( info_out_xml, 'basis', trim( string ) )
      write( string, '(i8)' ) qset%nkpt
      call xml_AddAttribute( info_out_xml, 'num_qpoints', trim( adjustl( string ) ) )
      iq = 0
      do ip = 1, size( parts_all )
        if( iq /= parts_all(ip)%iq ) then
          iq = parts_all(ip)%iq
          call xml_NewElement( info_out_xml, 'qpoint' )
          write( string, '(3g25.16)' ) qset%vkl(:, iq)
          call xml_AddAttribute( info_out_xml, 'vql', trim( adjustl( string ) ) )
          write( string, '(i8)' ) basis(iq)%nsym
          call xml_AddAttribute( info_out_xml, 'num_sym_q', trim( adjustl( string ) ) )
          write( string, '(i8)' ) parts_all(ip)%nkpt
          call xml_AddAttribute( info_out_xml, 'num_kpoints_q', trim( adjustl( string ) ) )
          write( string, '(i8)' ) basis(iq)%nirrep
          call xml_AddAttribute( info_out_xml, 'num_irreps', trim( adjustl( string ) ) )
          long = ''
          do ii = 1, basis(iq)%nirrep
            write( string, '(1x,i1)' ) basis(iq)%irreps(ii)%dim
            long = long // trim( adjustl( string ) )
          end do
          call xml_AddAttribute( info_out_xml, 'irrep_dim', long )
          call xml_EndElement( info_out_xml, 'qpoint' )
        end if
      end do
      call xml_EndElement( info_out_xml, 'symmetries' )
      write( info_out_unit, * )

      ! write parallelization settings
      totload = sum( [(parts_all(ip)%get_load(), ip=1, size( parts_all ))] )
      remload = sum( [(parts(ip)%get_load(), ip=1, size( parts ))] )
      call xml_NewElement( info_out_xml, 'parallelization' )
      write( string, '(i12)' ) size( parts )
      write( info_out_unit, '("current task / total (remaining) number of tasks:",T55,i5,"/",i8," (",a,")")' ) &
        ipart, size( parts_all ), trim( adjustl( string ) )
      write( string, '(i8)' ) ipart
      call xml_AddAttribute( info_out_xml, 'this_task', trim( adjustl( string ) ) )
      write( string, '(i8)' ) size( parts_all )
      call xml_AddAttribute( info_out_xml, 'num_total_tasks', trim( adjustl( string ) ) )
      write( string, '(i8)' ) size( parts )
      call xml_AddAttribute( info_out_xml, 'num_remaining_tasks', trim( adjustl( string ) ) )
      write( string, '(i12)' ) remload
      write( info_out_unit, '("load of current task / total (remaining) load:",T55,i5,"/",i8," (",a,")")' ) &
        parts(ipart)%get_load(), totload, trim( adjustl( string ) )
      write( string, '(i8)' ) parts(ipart)%get_load()
      call xml_AddAttribute( info_out_xml, 'this_load', trim( adjustl( string ) ) )
      write( string, '(i8)' ) totload
      call xml_AddAttribute( info_out_xml, 'total_load', trim( adjustl( string ) ) )
      write( string, '(i8)' ) remload
      call xml_AddAttribute( info_out_xml, 'remaining_load', trim( adjustl( string ) ) )
      write( info_out_unit, '("procs for current task / total number of procs:",T55,i5,"/",i8)' ) &
        parts(ipart)%mpi%procs, ph_numprocs
      write( string, '(i8)' ) parts(ipart)%mpi%procs
      call xml_AddAttribute( info_out_xml, 'num_this_task_procs', trim( adjustl( string ) ) )
      write( string, '(i8)' ) mpiglobal%procs
      call xml_AddAttribute( info_out_xml, 'num_total_procs', trim( adjustl( string ) ) )
      maxload = maxval( [(sum( [(parts(ip)%load_of_rank( irank ), ip=1, size( parts ))] ), &
                         irank=0, ph_numprocs-1)] )
      nidle = count( [(sum( [(parts(ip)%load_of_rank( irank ), ip=1, size( parts ))] ) == 0, &
                      irank=0, ph_numprocs-1)] )
      write( info_out_unit, '("average / maximum load per proc:",T55,i5,"/",i8)' ) &
        ceiling( dble( remload ) / ph_numprocs ), maxload
      write( info_out_unit, '("number of idle procs:",T55,i5)' ) nidle
      call xml_EndElement( info_out_xml, 'parallelization' )
    end subroutine ph_io_info_symmetries_and_parallelization

    !> This subroutine writes the displacement patterns of the given irrep
    !> to the general info output file (human readable and xml).
    subroutine ph_io_info_irrep_patterns( irr )
      use phonons_symmetry
      use mod_atoms, only: nspecies, natoms, idxas, spsymb
      !> irreducible representation
      type(irrep), intent(in) :: irr

      complex(dp), parameter :: eps = cmplx( 1e-64_dp, 1e-64_dp, dp )
      
      character(256) :: string
      integer :: id, is, ia, ias

      ! write displacement patterns
      call printbox( info_out_unit, '*', 'irreducible representations / atomic displacement patterns' )
      write( info_out_unit, * )
      call xml_NewElement( info_out_xml, 'irrep' )
      write( string, '(i8)' ) irr%dim
      call xml_AddAttribute( info_out_xml, 'dimension', trim( adjustl( string ) ) )
      do id = 1, irr%dim
        write( info_out_unit, '("irrep member ",i1)') id
        call xml_NewElement( info_out_xml, 'irrep_member')
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            write( info_out_unit, '(1x,a2,3x,i3,3x,sp,2f14.10,"i")' ) trim( adjustl( spsymb(is) ) ), ia, irr%pat(1, ias, id) + eps
            write( info_out_unit, '(12x,sp,2f14.10,"i")' ) irr%pat(2, ias, id) + eps
            write( info_out_unit, '(12x,sp,2f14.10,"i")' ) irr%pat(3, ias, id) + eps
            call xml_NewElement( info_out_xml, 'atom' )
            call xml_AddAttribute( info_out_xml, 'species', trim( adjustl( spsymb(is) ) ) )
            write( string, '(i8)' ) ia
            call xml_AddAttribute( info_out_xml, 'atom_index', trim( adjustl( string ) ) )
            write( string, '(sp,2g25.16)' ) irr%pat(1, ias, id) + eps
            call xml_AddAttribute( info_out_xml, 'displacement_pattern_x', trim( adjustl( string ) ) )
            write( string, '(sp,2g25.16)' ) irr%pat(2, ias, id) + eps
            call xml_AddAttribute( info_out_xml, 'displacement_pattern_y', trim( adjustl( string ) ) )
            write( string, '(sp,2g25.16)' ) irr%pat(3, ias, id) + eps
            call xml_AddAttribute( info_out_xml, 'displacement_pattern_z', trim( adjustl( string ) ) )
            call xml_EndElement( info_out_xml, 'atom' )
          end do
        end do
        write( info_out_unit, * )
        call xml_EndElement( info_out_xml, 'irrep_member' )
      end do
      call xml_EndElement( info_out_xml, 'irrep' )
    end subroutine ph_io_info_irrep_patterns

    !> This subroutine prints information and parallelization
    !> possibilities for a DFPT phonon calculation.
    subroutine ph_io_print_calc_info
      use m_getunit
      integer :: ip, iq, un, stat, totload, remload, load
      character(64) :: str

      if( mpiglobal%rank /= 0 ) return

      call getunit( un )
      open( un, file=trim( prefix_upper )//'_RUN_INFO.OUT', action='write', form='formatted', iostat=stat )

      call printbox( un, '=', 'DFPT PHONON CALCULATION INFORMATION' )
      write( un, * )
      if( ph_canonical ) then
        write( un, '("Canonical displacement patterns will be used.")' )
        write( un, '("k-points will not be reduced by symmetries.")' )
        write( un, '("Response quantities will not be symmetrized.")' )
      else
        write( un, '("Irreducible representations will be used.")' )
        if( dfpt_kset%isreduced ) then
          write( un, '("k-points will be reduced by symmetries.")' )
        else
          write( un, '("k-points will not be reduced by symmetries.")' )
        end if
        write( un, '("Response quantities will be symmetrized.")' )
      end if

      call printbox( un, '*', 'q-points and displacements' )
      write( un, * )
      write( un, '("number of q-points:    ",i12)' ) ph_qset%nkpt
      write( un, '("number of phonon modes:",i12)' ) size( ph_irrep_basis(1)%irreps(1)%pat(:, :, 1) )
      write( un, * )
      call write_qi_parts( un, ph_qset, ph_irrep_basis(1:), ph_parts_all, 0 )

      call printbox( un, '*', 'independent parts and computational load' )
      write( un, * )
      totload = sum( [(ph_parts_all(ip)%get_load(), ip=1, size( ph_parts_all ))] )
      remload = sum( [(ph_parts(ip)%get_load(), ip=1, size( ph_parts ))] )
      write( str, '(i12)' ) size( ph_parts )
      write( un, '("total (remaining) number of independent parts: ",i6," (",a,")")' ) size( ph_parts_all ), trim( adjustl( str ) )
      write( str, '(i12)' ) remload
      write( un, '("total (remaining) computational load:    ",i12," (",a,")")' ) totload, trim( adjustl( str ) )
      write( un, * )

      write( un, '("remaining parts")' )
      write( un, '("  ip   iq irrep dim load procs")' )
                  !"iiii iiii iiiii iii iiii iiiii"
      write( un, '("------------------------------")' )
      do ip = 1, size( ph_parts )
        write( un, '(i4,x,i4,x,i5,x,i3,x,i4,x,i5)' ) &
          ip, ph_parts(ip)%iq, ph_parts(ip)%iirrep, ph_parts(ip)%dirrep, &
          ph_parts(ip)%get_load(), ph_parts(ip)%num_ranks
      end do
      write( un, * )

      write( un, '("done parts")' )
      write( un, '("  iq irrep dim load")' )
                  !"iiii iiiii iii iiii"
      write( un, '("-------------------")' )
      do ip = 1, size( ph_parts_all )
        if( .not. ph_parts_done(ph_parts_all(ip)%iq, ph_parts_all(ip)%iirrep) ) cycle
        write( un, '(i4,x,i5,x,i3,x,i4,x,a4,x,i5)' ) &
          ph_parts_all(ip)%iq, ph_parts_all(ip)%iirrep, ph_parts_all(ip)%dirrep, &
          ph_parts_all(ip)%get_load()
      end do
      write( un, * )

      call printbox( un, '*', 'parallelization settings' )
      write( un, * )
      write( un, '("number of processes:             ",i12)' ) ph_numprocs
      write( un, '("optimal average load per process:",f12.1)' ) dble( remload ) / ph_numprocs
      write( un, '("minimum load per process:        ",i12)' ) &
        minval( [(sum( pack( [(ph_parts(ip)%load_of_rank( iq ), ip=1, size( ph_parts ))], &
                             [(ph_parts(ip)%is_my_rank( iq ), ip=1, size( ph_parts ))] ) ), &
                 iq=0, ph_numprocs-1)] )
      write( un, '("maximum load per process:        ",i12)' ) &
        maxval( [(sum( pack( [(ph_parts(ip)%load_of_rank( iq ), ip=1, size( ph_parts ))], &
                             [(ph_parts(ip)%is_my_rank( iq ), ip=1, size( ph_parts ))] ) ), &
                 iq=0, ph_numprocs-1)] )
      write( un, '("number of idle processes:        ",i12)' ) &
        count( [(all( [(ph_parts(ip)%load_of_rank( iq ), ip=1, size( ph_parts ))] <= 0 ), &
                iq=0, ph_numprocs-1)] )

      write( un, * )
      write( un, '("proc load (imbalance)")' )
                  !"iiii iiii (sfff.f)"
      write( un, '("----------------------")' )
      do iq = 0, ph_numprocs-1
        load = sum( [(ph_parts(ip)%load_of_rank( iq ), ip=1, size( ph_parts ))] )
        if( load == 0 ) then
          write( un, '(i4,x,"idle")' ) iq
        else
          write( str, '(sp,f7.1)' ) dble( load ) - dble( remload ) / ph_numprocs
          write( un, '(i4,x,i4,x,"(",a,")")' ) iq, load, trim( adjustl( str ) )
        end if
      end do

      write( un, * )
      write( un, '("execution schedule")' )
      call write_schedule( un, ph_schedule )

      close( un, iostat=stat )
    end subroutine ph_io_print_calc_info

    !> Write constant part of force response (independent of displacement
    !> patterns and density / potential response) to file.
    subroutine ph_io_write_dforce_const( dforce, success )
      use m_getunit
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, spsymb
      !> constant part of fore response
      complex(dp), intent(in) :: dforce(3, natmtot, 3)
      !> `.true.` on success
      logical, intent(out) :: success

      integer :: un, stat, is, ia, ias, ip
      complex(dp) :: row(3)

      call getunit(un)

      open( un, file=trim( prefix_upper )//'_DFORCE_CONST.OUT', action='write', form='formatted', iostat=stat )
      success = (stat == 0)
      if( .not. success ) return

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          write( un, '("#",x,a,x,i3)' ) trim( adjustl( spsymb(is) ) ), ia
          do ip = 1, 3
            row = dforce(ip, ias, :)
            where( abs( dble( row ) ) < 1e-12 ) row = cmplx( 0.0_dp, aimag( row ), dp )
            where( abs( aimag( row ) ) < 1e-12 ) row = cmplx( dble( row ), 0.0_dp, dp )
            write( un, '(3(2g21.14,3x))', iostat=stat ) row
            success = success .and. (stat == 0)
          end do
        end do
      end do

      close( un )
    end subroutine ph_io_write_dforce_const

    !> Read constant part of force response (independent of displacement
    !> patterns and density / potential response) from file.
    subroutine ph_io_read_dforce_const( dforce, success )
      use m_getunit
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      !> constant part of fore response
      complex(dp), intent(out) :: dforce(3, natmtot, 3)
      !> `.true.` on success
      logical, intent(out) :: success

      integer :: un, stat, is, ia, ias, ip
      character(256) :: line
      real(dp) :: row(6)

      call getunit(un)

      open( un, file=trim( prefix_upper )//'_DFORCE_CONST.OUT', action='read', form='formatted', iostat=stat )
      success = (stat == 0)
      if( .not. success ) return

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          read( un, * ) line
          do ip = 1, 3
            read( un, *, iostat=stat ) row
            dforce(ip, ias, :) = cmplx( row(1::2), row(2::2), dp )
            success = success .and. (stat == 0)
          end do
        end do
      end do

      close( un )
    end subroutine ph_io_read_dforce_const

    !> Print q-vectors and symmetries and irrep information for all parts to unit.
    subroutine write_qi_parts( un, qset, basis, parts, ipart )
      use phonons_symmetry, only: irrep_basis
      use phonons_parallelization, only: ph_part
      use mod_kpointset, only: k_set
      !> output unit
      integer, intent(in) :: un
      !> set of \({\bf q}\) points
      type(k_set), intent(in) :: qset
      !> irrep basis per \({\bf q}\) point
      type(irrep_basis), intent(in) :: basis(qset%nkpt)
      !> array of (remaining) independent parts
      type(ph_part), intent(in) :: parts(:)
      !> index of current part
      integer, intent(in) :: ipart

      integer :: iq, ip, ii

      write( un, '("  iq  q-vector (lattice coordinates)  N_sym  N_k  N_irrep (size)")' )
                  !"iiii (ff.ffffff,ff.ffffff,ff.ffffff)  iiiii iiii  iiiiiii (i,i,...)"
      write( un, '("--------------------------------------------------------------------------------")' )
      do iq = 1, qset%nkpt
        do ip = 1, size( parts )
          if( parts(ip)%iq == iq ) exit
        end do
        if( ip > size( parts ) ) cycle
        write( un, '(i4," (",2(f9.6,","),f9.6,")  ",i5,x,i4,2x,i7,x)', advance='no' ) &
          iq, qset%vkl(:, iq), basis(iq)%nsym, parts(ip)%nkpt, basis(iq)%nirrep
        write( un, '("(")', advance='no' )
        do ii = 1, basis(iq)%nirrep
          write( un, '(i1)', advance='no' ) basis(iq)%irreps(ii)%dim
          if( ipart > 0 .and. ipart <= size( parts ) ) then
            if( parts(ipart)%iq == iq .and. parts(ipart)%iirrep == ii ) write( un, '("*")', advance='no' )
          end if
          if( ii < basis(iq)%nirrep ) write( un, '(",")', advance='no' )
        end do
        write( un, '(")")' )
      end do
    end subroutine write_qi_parts

    subroutine write_schedule( un, schedule )
      integer, intent(in) :: un
      integer, intent(in) :: schedule(:,:)

      integer, parameter :: block_width = 2
      integer, parameter :: load_step = 5
      character(*), parameter :: char_list = '#:%@x0+o'

      integer :: nproc, maxload, npart, iproc, ipart, fst, lst, width, i, j, k
      character(1) :: this_char, prev_char
      character(16) :: fmt

      integer, allocatable :: blocks(:,:)
      character(1), allocatable :: char_schedule(:,:)
      character(:), allocatable :: line, best_chars

      nproc = size( schedule, dim=1 )
      maxload = size( schedule, dim=2 )
      npart = maxval( schedule )

      allocate( char_schedule(nproc, block_width*maxload), source=' ' )
      allocate( blocks(4, npart) )
      blocks(1, :) =  huge( 1 )
      blocks(2, :) = -huge( 1 )
      blocks(3, :) =  huge( 1 )
      blocks(4, :) = -huge( 1 )

      do iproc = 1, nproc
        fst = 1; lst = 1
        prev_char = ' '
        do while( lst <= maxload )
          ! next block
          ipart = schedule(iproc, fst)
          do while( schedule(iproc, lst) == ipart )
            lst = lst + 1
            if( lst > maxload ) exit
          end do
          lst = lst - 1
          if( ipart > 0 .and. ipart <= npart ) then
            blocks(1, ipart) = min( iproc, blocks(1, ipart) )
            blocks(2, ipart) = max( iproc, blocks(2, ipart) )
            blocks(3, ipart) = min( fst, blocks(3, ipart) )
            blocks(4, ipart) = max( lst, blocks(4, ipart) )
          end if
          ! character range
          i = (fst - 1) * block_width + 1
          j = lst * block_width
          ! find character
          if( ipart == 0 ) then
            this_char = ' '
          else
            if( iproc == 1 ) then
              best_chars = find_best_chars( char_list, prev_char )
            else
              line = ''
              do k = max(1, i-1), min(block_width*maxload, j+1)
                line = line//char_schedule(iproc-1, k)
              end do
              if( all( schedule(iproc-1, fst:lst) == ipart ) ) then
                best_chars = char_schedule(iproc-1, i)
              else
                best_chars = find_best_chars( char_list, line )
                best_chars = find_best_chars( best_chars, prev_char )
              end if
            end if
            this_char = best_chars(1:1)
            prev_char = this_char
          end if
          ! assign to char array
          char_schedule(iproc, i:j) = this_char

          fst = lst + 1
          lst = fst
        end do
      end do

      ! assign names to blocks
      do ipart = 1, npart
        if( blocks(2, ipart) <= 0 .or. blocks(1, ipart) > nproc ) cycle
        write( fmt, '(i12)' ) ipart
        line = trim( adjustl( fmt ) )
        this_char = char_schedule(blocks(1, ipart), blocks(4, ipart)*block_width)
        i = (blocks(1, ipart) + blocks(2, ipart)) / 2
        j = (blocks(3, ipart) - 1) * block_width
        width = (blocks(4, ipart) - blocks(3, ipart) + 1) * block_width
        if( len( line ) < width ) line = 'P'//line
        do while( len( line ) < width )
          if( mod( len( line ), 2 ) == 0 ) then
            line = this_char//line
          else
            line = line//this_char
          end if
        end do
        do k = 1, width
          char_schedule(i, j+k) = line(k:k)
        end do
      end do

      ! write to unit
      write( fmt, '(i12)' ) load_step*block_width - 1
      fmt = "(i" // trim( adjustl( fmt ) ) // ",""|"")"
      write( un, '(" proc|")', advance='no' )
      do i = load_step, maxload, load_step
        write( un, trim(fmt), advance='no' ) i
      end do
      write( un, * )
      do iproc = 1, nproc
        write( un, '(i5,"|")', advance='no' ) iproc-1
        do i = 1, block_width*maxload
          write( un, '(a)', advance='no' ) char_schedule(iproc, i)
        end do
        write( un, * )
      end do

      contains

        function find_best_chars( allowed, forbidden ) result( chars )
          character(*), intent(in) :: allowed, forbidden
          character(:), allocatable :: chars

          integer :: i, j
          integer, allocatable :: cnt(:)

          allocate( cnt(len(allowed)), source=0 )

          do j = 1, len(forbidden)
            do i = 1, len(allowed)
              if( allowed(i:i) == forbidden(j:j) ) then
                cnt(i) = cnt(i) + 1
                exit
              end if
            end do
          end do

          chars = ''
          j = minval( cnt )
          do i = 1, len(allowed)
            if( cnt(i) == j ) chars = chars//allowed(i:i)
          end do
          
          deallocate( cnt )
        end function find_best_chars
    end subroutine write_schedule

end module phonons_inout
