!> This module contains variable types describing independent parts
!> of a DFPT phonon calculation and procedures to distribute them
!> among multiple MPI processes.
!> 
!> The set of \({\bf q}\) points and the respective irreducible representations
!> (irreps) \(I\) mainly determine the parallelization settings of the phonon
!> calculation. See the module [[phonons_symmetry(module)]] for more information about
!> irreps. Each pair \(({\bf q},I)\) forms an independent part of the full calculation.
!> The total computational load \(\texttt{load}({\bf q},I)\) of each part is proportional
!> to the product of the number of irrep members and the number of \({\bf k}\) points, i.e.,
!> \[ \texttt{load}({\bf q},I) = d_I \times |\text{IBZ}({\bf q})| \;. \]
!> See the subroutine [[ph_par_distribute(subroutine)]] for more details about the
!> distribution of the total load \(\texttt{total_load} = \sum_{{\bf q},I} \texttt{load}({\bf q},I)\)
!> among all available processe.
module phonons_parallelization
  use precision, only: dp
  use asserts, only: assert
  use modmpi

  implicit none
  private

  ! DERIVED DATA TYPES
  !> variable type describing an independent part of the phonon calculation
  !> and its parallelization settings
  type ph_part
    !> index of \({\bf q}\) point
    integer :: iq
    !> index of irrep \(I\)
    integer :: iirrep
    !> number of \({\bf k}\) points in \(\text{IBZ}({\bf q})\)
    integer :: nkpt
    !> dimension \(d_I\) of irrep
    integer :: dirrep
    !> global range of ranks that do this part
    integer :: first_rank, last_rank
    !> number of ranks that to this part
    integer :: num_ranks
    !> first and last irrep member per rank
    integer, allocatable :: first_d(:), last_d(:)
    !> first and last \({\bf k}\) point per rank
    integer, allocatable :: first_k(:), last_k(:)
    !> MPI communicator of all ranks working on this part
    type(mpiinfo) :: mpi
    !> MPI communicator of all ranks working on the first \({\bf k}\) point
    type(mpiinfo) :: mpik
    !> MPI communicators of all ranks working on the same irrep member
    type( mpiinfo), allocatable :: mpid(:)

    contains

      !> simple setter functions
      procedure, public :: set_first_rank, set_last_rank
      !> get computational load of part
      procedure, public :: get_load
      !> check if global rank is working on this part
      procedure, public :: is_my_rank
      !> first and last \({\bf k}\) point index the given rank is working on
      procedure, public :: first_k_of_rank, last_k_of_rank
      !> first and last irrep member the given rank is working on
      procedure, public :: first_d_of_rank, last_d_of_rank
      !> offset for a given \({\bf k}\) point and irrep member for file operations
      procedure, public :: get_dk_offset
      !> computational load for a given rank
      procedure, public :: load_of_rank 
      !> distribute the total load of that part among contributing processes
      procedure, public :: distribute_load
  end type ph_part
  ! constructor
  interface ph_part
    module procedure :: ph_part_gen
  end interface

  public :: ph_part
  public :: ph_par_distribute, ph_par_parts_from_schedule
  public :: ph_par_zscatter, ph_par_zgather

  contains

    !> This subroutine finds the independent parts of the full phonon calculation
    !> and distributes them among the given number of processes.
    !>
    !> The full phonon calculation is made of independent parts that can be computed
    !> in parallel without interdependencies. Each pair \(({\bf q},I)\) of a \({\bf q}\)-point
    !> and an irrep \(I\) forms such an independent part.
    !> This subroutine distributes the independent parts among the given number of processes
    !> `numprocs`. Some (precomputed) parts can be excluded using the map `qi_done`.
    !> If `qi_done(iq, iirrep) = .true.`, this part will not be distributed.
    !>
    !> Each part has a computational load which is determined by the product of the
    !> number of \({\bf k}\)-points in the \({\rm IBZ}({\bf q})\) and the number of
    !> members \(\mu\) in the irrep. In principle, parallelization can be done over
    !> \(({\bf q},I)\) and \(({\bf k},\mu)\) pairs simultaneously. Parallelization over
    !> part internal load \(({\bf k},\mu)\) is prefered over parallelization over 
    !> parts \(({\bf q},I)\).
    !>
    !> The resulting distribution pattern will be returned in `schedule`.
    !> See [[gen_schedule(subroutine)]] for details.
    subroutine ph_par_distribute( qset, basis, numprocs, minppp, maxppp, qi_done, parts, schedule, &
        all_parts )
      use phonons_symmetry, only: irrep_basis
      use mod_kpointset, only: k_set, generate_k_vectors, delete_k_vectors
      use mod_symmetry, only: maxsymcrys, nsymcrys, lsplsymc
      use mod_atoms, only: natmtot
      use modinput
      !> set of \({\bf q}\) points
      type(k_set), intent(in) :: qset
      !> irrep basis for each \({\bf q}\) point
      type(irrep_basis), intent(in) :: basis(qset%nkpt)
      !> number of processes over which the calculation shoud be distributed
      integer, intent(in) :: numprocs
      !> minimum number of processes per part
      integer, intent(in) :: minppp
      !> maximum number of processes per part
      integer, intent(in) :: maxppp
      !> flag if \(({\bf q},I)\) pair is already done
      logical, intent(in) :: qi_done(qset%nkpt, 3*natmtot)
      !> list of remaining independent parts of the calculation
      type(ph_part), allocatable, intent(out) :: parts(:)
      !> schedule telling which process is working on which part in which order
      integer, allocatable, intent(out) :: schedule(:,:)
      !> list of all independent parts of the calculation
      type(ph_part), allocatable, optional, intent(out) :: all_parts(:)

      integer :: i, j, k, iq, ip, np, np_all, iirrep, dirrep, iproc
      integer :: nsymcrys_, lsplsymc_(maxsymcrys)
      type(k_set) :: kset
      
      integer, allocatable :: pidx(:), piq(:), pirrep(:), pnkpt(:), pdirrep(:), pload(:)

      call assert( minppp <= numprocs, &
        'Minimal number of processes per part `minppp` must not be greater than number of processes `numprocs`.' )

      ! allocate temporary arrays for parts information
      allocate( pidx(   3*natmtot*qset%nkpt), source=0 )
      allocate( piq(    3*natmtot*qset%nkpt) )
      allocate( pirrep( 3*natmtot*qset%nkpt) )
      allocate( pnkpt(  3*natmtot*qset%nkpt) )
      allocate( pdirrep(3*natmtot*qset%nkpt) )
      allocate( pload(0:3*natmtot*qset%nkpt), source=0 )

      np = 0
      np_all = 0

      ! store global variables
      nsymcrys_ = nsymcrys
      lsplsymc_ = lsplsymc

      ! start loop over q-points
      do iq = 1, qset%nkpt
        ! find IBZ(q)
        nsymcrys = basis(iq)%nsym
        do i = 1, basis(iq)%nsym
          lsplsymc(i) = lsplsymc_(basis(iq)%isym(i))
        end do
        call generate_k_vectors( kset, qset%bvec, &
               input%groundstate%ngridk, &
               input%groundstate%vkloff, &
               input%groundstate%reducek, &
               uselibzint=.false. )

        ! start loop over irreps
        do iirrep = 1, basis(iq)%nirrep
          dirrep = basis(iq)%irreps(iirrep)%dim
          np_all = np_all + 1
          piq(np_all) = iq
          pirrep(np_all) = iirrep
          pnkpt(np_all) = kset%nkpt
          pdirrep(np_all) = dirrep
          pload(np_all) = pnkpt(np_all) * pdirrep(np_all)
          ! add to list of remaining parts if not yet done
          if( .not. qi_done(iq, iirrep) ) then
            np = np + 1
            pidx(np) = np_all
          end if
        end do
         ! end loop over irreps
        call delete_k_vectors( kset )
      end do
      ! end loop over q-points

      ! assign all parts if needed
      if( present( all_parts ) ) then
        if( allocated( all_parts ) ) deallocate( all_parts )
        allocate( all_parts(np_all) )
        do ip = 1, np_all
          all_parts(ip) = ph_part( piq(ip), pirrep(ip), pdirrep(ip), pnkpt(ip) )
        end do
      end if

      ! return if there is nothing to do
      if( np == 0 ) then
        deallocate( pidx, piq, pirrep, pnkpt, pdirrep, pload )
        allocate( parts(0) )
        allocate( schedule(1, 0) )
        return
      end if

      ! bring remaining parts to the front
      piq(1:np) = piq(pidx(1:np))
      pirrep(1:np) = pirrep(pidx(1:np))
      pnkpt(1:np) = pnkpt(pidx(1:np))
      pdirrep(1:np) = pdirrep(pidx(1:np))
      pload(1:np) = pload(pidx(1:np))

      ! generate schedule
      schedule = gen_schedule( numprocs, np, pload(1:np), minppp, maxppp )

      ! create and distribute parts
      if( allocated( parts ) ) deallocate( parts )
      allocate( parts(np) )
      do ip = 1, np
        parts(ip) = ph_part( piq(ip), pirrep(ip), pdirrep(ip), pnkpt(ip) )
        do iproc = 1, numprocs
          if( any( schedule(iproc, :) == ip ) ) then
            call parts(ip)%set_first_rank( iproc-1 )
            exit
          end if
        end do
        do iproc = numprocs, 1, -1
          if( any( schedule(iproc, :) == ip ) ) then
            call parts(ip)%set_last_rank( iproc-1 )
            exit
          end if
        end do
        call parts(ip)%distribute_load()
      end do

      ! create local MPI communicators
#ifdef MPI
      do ip = 1, np
        ! communicator with all procs working on this part
        k = 1
        call barrier
        i = MPI_UNDEFINED
        if( parts(ip)%is_my_rank( mpiglobal%rank ) ) i = k
        call MPI_comm_split( mpiglobal%comm, i, mpiglobal%rank, j, mpiglobal%ierr )
        call set_mpi( parts(ip)%mpi, j )
        ! communicator with all procs working on the first k-point in this part
        k = 1
        call barrier
        i = MPI_UNDEFINED
        if( (k >= parts(ip)%first_k_of_rank( mpiglobal%rank )) .and. (k <= parts(ip)%last_k_of_rank( mpiglobal%rank )) ) i = k
        call MPI_comm_split( mpiglobal%comm, i, mpiglobal%rank, j, mpiglobal%ierr )
        call set_mpi( parts(ip)%mpik, j )
        ! communicators with all procs working on the same irrep members in this part
        do k = 1, pdirrep(ip)
          call barrier
          i = MPI_UNDEFINED
          if( (k >= parts(ip)%first_d_of_rank( mpiglobal%rank )) .and. (k <= parts(ip)%last_d_of_rank( mpiglobal%rank )) ) i = k
          call MPI_comm_split( mpiglobal%comm, i, mpiglobal%rank, j, mpiglobal%ierr )
          call set_mpi( parts(ip)%mpid(k), j )
        end do
      end do
#endif

      ! restore global variables
      nsymcrys = nsymcrys_
      lsplsymc = lsplsymc_

      ! free memory
      deallocate( pidx, piq, pirrep, pnkpt, pdirrep, pload )
    end subroutine ph_par_distribute

    !> Get parts the process calling process works on from the
    !> given row of the schedule (see [[gen_schedule(function)]]).
    function ph_par_parts_from_schedule( schedule ) result( parts )
      integer, intent(in) :: schedule(:)
      integer, allocatable :: parts(:)

      integer :: i, ip, np
      integer, allocatable :: res(:)

      allocate( res(size( schedule )), source=0 )
      np = 0; ip = 0
      do i = 1, size( schedule )
        if( schedule(i) /= ip ) then
          ip = schedule(i)
          if( ip == 0 ) cycle
          np = np + 1
          res(np) = ip
        end if
      end do

      if( allocated( parts ) ) deallocate( parts )
      allocate( parts, source=res(1:np) )

      deallocate( res )
    end function ph_par_parts_from_schedule

    !> simple wrapper for scattering complex array using MPI
    subroutine ph_par_zscatter( data, root, mpi, rlen, &
        roff )
      !> complex array to scatter over processes
      complex(dp), intent(inout) :: data(*)
      !> rank of root process
      integer, intent(in) :: root
      !> mpi opbject / communicator
      type(mpiinfo) :: mpi
      !> number of elements per process
      integer, intent(in) :: rlen(0:)
      !> element offset per process
      integer, optional, intent(in) :: roff(0:)

#ifdef MPI
      if( (mpi%rank < 0) .or. (mpi%procs <= 0) ) return

      ! variable record lengths
      if( size( rlen ) > 1 ) then
        call terminate_if_false( present( roff ), '(ph_par_zscatter) &
          Offsets must be provided for variable record lengths.' )
        call terminate_if_false( (size( rlen ) >= mpi%procs) .and. (size( roff ) >= mpi%procs), '(ph_par_zscatter) &
          Less records that processes.' )
        if( mpi%rank == root ) then
          call MPI_Scatterv( data, rlen, roff, MPI_DOUBLE_COMPLEX, MPI_IN_PLACE, rlen(mpi%rank), MPI_DOUBLE_COMPLEX, &
                 root, mpi%comm, mpi%ierr )
        else
          call MPI_Scatterv( data, rlen, roff, MPI_DOUBLE_COMPLEX, data, rlen(mpi%rank), MPI_DOUBLE_COMPLEX, &
                 root, mpi%comm, mpi%ierr )
        end if
      else
        if( mpi%rank == root ) then
          call MPI_Scatter( data, rlen(0), MPI_DOUBLE_COMPLEX, MPI_IN_PLACE, rlen(0), MPI_DOUBLE_COMPLEX, &
                 root, mpi%comm, mpi%ierr )
        else
          call MPI_Scatter( data, rlen(0), MPI_DOUBLE_COMPLEX, data, rlen(0), MPI_DOUBLE_COMPLEX, &
                 root, mpi%comm, mpi%ierr )
        end if
      end if
#endif
    end subroutine ph_par_zscatter

    !> simple wrapper for gathering complex array using MPI
    subroutine ph_par_zgather( data, root, mpi, rlen, &
        roff )
      !> complex array to gtaher from processes
      complex(dp), intent(inout) :: data(*)
      !> rank of root process
      integer, intent(in) :: root
      !> mpi object / communicator
      type(mpiinfo) :: mpi
      !> number of elements per process
      integer, intent(in) :: rlen(0:)
      !> element offset per process
      integer, optional, intent(in) :: roff(0:)

#ifdef MPI
      if( (mpi%rank < 0) .or. (mpi%procs <= 0) ) return

      ! variable record lengths
      if( size( rlen ) > 1 ) then
        call terminate_if_false( present( roff ), '(ph_par_zgather) &
          Offsets must be provided for variable record lengths.' )
        call terminate_if_false( (size( rlen ) >= mpi%procs) .and. (size( roff ) >= mpi%procs), '(ph_par_zgather) &
          Less records that processes.' )
        if( mpi%rank == root ) then
          call MPI_Gatherv( MPI_IN_PLACE, rlen(mpi%rank), MPI_DOUBLE_COMPLEX, data, rlen, roff, MPI_DOUBLE_COMPLEX, &
                 root, mpi%comm, mpi%ierr )
        else
          call MPI_Gatherv( data, rlen(mpi%rank), MPI_DOUBLE_COMPLEX, data, rlen, roff, MPI_DOUBLE_COMPLEX, &
                 root, mpi%comm, mpi%ierr )
        end if
      else
        if( mpi%rank == root ) then
          call MPI_Gather( MPI_IN_PLACE, rlen(0), MPI_DOUBLE_COMPLEX, data, rlen(0), MPI_DOUBLE_COMPLEX, &
                 root, mpi%comm, mpi%ierr )
        else
          call MPI_Gather( data, rlen(0), MPI_DOUBLE_COMPLEX, data, rlen(0), MPI_DOUBLE_COMPLEX, &
                 root, mpi%comm, mpi%ierr )
        end if
      end if
#endif
    end subroutine ph_par_zgather

    !> `ph_part` constructor
    function ph_part_gen( iq, iirrep, dirrep, nkpt ) &
        result( part )
      !> index of \({\bf q}\) point
      integer, intent(in) :: iq
      !> index of irrep
      integer, intent(in) :: iirrep
      !> dimension of the irrep
      integer, intent(in) :: dirrep
      !> number of \({\bf k}\) points in \(\text{IBZ}({\bf q})\)
      integer, intent(in) :: nkpt
      type(ph_part) :: part

      part%iq     = iq
      part%iirrep = iirrep
      part%dirrep = dirrep
      part%nkpt   = nkpt
      part%mpi    = mpiglobal
      part%mpik   = mpiglobal
      allocate( part%mpid(part%dirrep), source=mpiglobal )
      call part%set_first_rank(0)
      call part%set_last_rank(0)
    end function

    !> get load of part
    function get_load( this ) result( load )
      class(ph_part) :: this
      integer :: load
      load = this%nkpt * this%dirrep
    end function

    !> set first process that does this part
    subroutine set_first_rank( this, rank )
      class(ph_part), intent(inout) :: this
      !> global rank of process
      integer, intent(in) :: rank
      
      this%first_rank = max( 0, rank )
      this%num_ranks = this%last_rank - this%first_rank + 1
    end subroutine

    !> set last process that does this part
    subroutine set_last_rank( this, rank )
      class(ph_part), intent(inout) :: this
      !> global rank of process
      integer, intent(in) :: rank
      
      this%last_rank = max( 0, rank )
      this%num_ranks = this%last_rank - this%first_rank + 1
    end subroutine

    !> check if this part is assigned to requested process
    function is_my_rank( this, rank ) result( l )
      class(ph_part)  :: this
      !> global rank of process
      integer, intent(in) :: rank
      logical :: l

      l = (rank >= this%first_rank .and. rank <= this%last_rank)
    end function

    !> first irrep member of this part assigned to requested process
    function first_d_of_rank( this, rank ) result( id )
      class(ph_part) :: this
      !> global rank of process
      integer, intent(in) :: rank
      integer :: id
      id = 0
      if( this%is_my_rank( rank ) ) id = this%first_d(rank)
    end function

    !> last irrep member of this part assigned to requested process
    function last_d_of_rank( this, rank ) result( id )
      class(ph_part) :: this
      !> global rank of process
      integer, intent(in) :: rank
      integer :: id
      id = -1
      if( this%is_my_rank( rank ) ) id = this%last_d(rank)
    end function

    !> first \({\bf k}\) point of this part assigned to requested process
    function first_k_of_rank( this, rank ) result( ik )
      class(ph_part) :: this
      !> global rank of process
      integer, intent(in) :: rank
      integer :: ik
      ik = 0
      if( this%is_my_rank( rank ) ) ik = this%first_k(rank)
    end function

    !> last \({\bf k}\) point of this part assigned to requested process
    function last_k_of_rank( this, rank ) result( ik )
      class(ph_part) :: this
      !> global rank of process
      integer, intent(in) :: rank
      integer :: ik
      ik = -1
      if( this%is_my_rank( rank ) ) ik = this%last_k(rank)
    end function

    !> get offset in an array with last dimensions \(N_{\bf k}({\bf q})\) and \(d_I\) 
    !> for a given pair of \({\bf k}\) point and irrep member
    function get_dk_offset( this, id, ik ) result( off )
      class(ph_part) :: this
      !> member of irrep
      integer, intent(in) :: id
      !> k-point index
      integer, intent(in) :: ik
      integer :: off
      off = (id - 1) * this%nkpt + ik - 1
    end function

    !> get computational load (coming from this part) assigned to requested process
    function load_of_rank( this, rank ) result( load )
      class(ph_part) :: this
      !> global rank of process
      integer, intent(in) :: rank
      integer :: load
      load = (this%last_d_of_rank( rank ) - this%first_d_of_rank( rank ) + 1) * &
             (this%last_k_of_rank( rank ) - this%first_k_of_rank( rank ) + 1)
    end function

    !> This subroutine distributes the computational load of this part as evenly as possible
    !> among all processes that work on this part.
    !> If multiple processes are working on this part, then this routine descides
    !> whether to first parallelize loops over irrep members or over \({\bf k}\) points.
    subroutine distribute_load( this )
      class(ph_part) :: this

      integer :: k_load, d_load
      integer :: i, j, k, l, m ,n
      integer :: dsize, ksize, nd, nk, id, ik, irank

      if( allocated( this%first_d ) ) deallocate( this%first_d )
      allocate( this%first_d(this%first_rank:this%last_rank), source=0 )
      if( allocated( this%last_d ) ) deallocate( this%last_d )
      allocate( this%last_d(this%first_rank:this%last_rank), source=-1 )
      if( allocated( this%first_k ) ) deallocate( this%first_k )
      allocate( this%first_k(this%first_rank:this%last_rank), source=0 )
      if( allocated( this%last_k ) ) deallocate( this%last_k )
      allocate( this%last_k(this%first_rank:this%last_rank), source=-1 )

      ! get load if k-loop is first parallelized
      i = ceiling( dble( this%nkpt ) / this%num_ranks )   ! max k-points per rank
      j = max( 1, this%num_ranks / this%nkpt )            ! max ranks for d-loop
      k = ceiling( dble( this%dirrep ) / j )              ! max members per rank
      k_load = k * i                                      ! max load per rank

      ! get load if d-loop is first parallelized
      l = ceiling( dble( this%dirrep ) / this%num_ranks ) ! max members per rank
      m = max( 1, this%num_ranks / this%dirrep )          ! max ranks for k-loop
      n = ceiling( dble( this%nkpt ) / m )                ! max k-points per rank
      d_load = n * l                                      ! max load per rank

      ! set number and size of blocks
      if( k_load <= d_load ) then
        dsize = k; ksize = i
      else
        dsize = l; ksize = n
      end if
      nd = ceiling( dble( this%dirrep ) / dsize )
      nk = ceiling( dble( this%nkpt ) / ksize )

      ! distribute blocks
      irank = this%first_rank
      id = 0; j = nd
      do while( id < this%dirrep )
        m = ceiling( dble( this%dirrep - id ) / j )
        ik = 0; k = nk
        do while( ik < this%nkpt )
          n = ceiling( dble( this%nkpt - ik ) / k )
          this%first_d(irank) = id + 1
          this%last_d(irank ) = id + m
          this%first_k(irank) = ik + 1
          this%last_k(irank)  = ik + n
          irank = irank + 1
          ik = ik + n
          k = k - 1
        end do
        id = id + m
        j = j - 1
      end do
    end subroutine

    !> set MPI communicator, rank and size
    subroutine set_mpi( mpi, comm )
      !> MPI object
      type(mpiinfo) :: mpi
      !> local MPI communicator
      integer, intent(in) :: comm

#ifdef MPI
      mpi%comm = comm
      if( comm /= MPI_COMM_NULL ) then
        call mpi_comm_size( mpi%comm, mpi%procs, mpi%ierr )
        call mpi_comm_rank( mpi%comm, mpi%rank, mpi%ierr )
      else
        mpi%procs = 0
        mpi%rank = -1
      end if
#endif
    end subroutine

end module phonons_parallelization
