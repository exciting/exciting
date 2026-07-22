!> exciting wrappers for mpi_gather and mpi_gatherv
module mod_mpi_gather
  use mod_mpi_env, only: mpiinfo
#ifdef MPI
  use mpi_f08, only: mpi_gather, mpi_gatherv, mpi_comm, MPI_Datatype, &
    MPI_DOUBLE_COMPLEX, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, MPI_REAL, MPI_INTEGER8, MPI_LOR
  use modmpi, only: terminate_mpi_env
#ifdef USE_MPI4_LARGE_COUNTS
  use mpi_f08, only: MPI_COUNT_KIND, MPI_ADDRESS_KIND
#endif

#endif
  use precision, only: dp, sp, i32, long_int

  implicit none 
  private
  
  public :: xmpi_gather, xmpi_gatherv, xmpi_get_gatherv_counts

  interface xmpi_gather
    module procedure :: mpi_gather_integer_i32_r1
    module procedure :: mpi_gather_real_dp_r2
    module procedure :: mpi_gather_complex_dp_r3
  end interface

  interface xmpi_gatherv
    module procedure :: mpi_gatherv_integer_i32_r1
    module procedure :: mpi_gatherv_real_dp_r2
    module procedure :: mpi_gatherv_complex_dp_r3
  end interface

contains

  !> Helper subroutine to gather array sizes from all ranks
  subroutine xmpi_get_gatherv_counts( mpi_env, local_send_count, receive_counts )
    type(mpiinfo), intent(in) :: mpi_env
    integer(long_int), intent(in) :: local_send_count
    integer(long_int), allocatable, intent(out) :: receive_counts(:)
#ifdef MPI
    integer(i32) :: ierr
    if ( mpi_env%is_root ) then
      allocate( receive_counts(mpi_env%procs) )
      call mpi_gather( local_send_count, 1, MPI_INTEGER8, receive_counts, 1, &
        MPI_INTEGER8, mpi_env%root, mpi_comm(mpi_env%comm), ierr )
    else
      allocate( receive_counts(1) ) 
      call mpi_gather( local_send_count, 1, MPI_INTEGER8, receive_counts, 1, &
        MPI_INTEGER8, mpi_env%root, mpi_comm(mpi_env%comm), ierr )
    end if
#else
    allocate( receive_counts(1), source=local_send_count )
#endif
  end subroutine

! === RANK 1 (1D Arrays) - e.g. Integer i32 ===
#define TYPE1 integer
#define PRECISION1 i32
#define MPIDATATYPE1 MPI_INTEGER
#define RANK_DECL (:)
#define FULL_EXT integer_i32_r1
#define ALLOC_GATHER(SND, RCV, PROCS) allocate( RCV(size(SND) * PROCS) )
#define ALLOC_GATHERV(SND, RCV, TOTAL) allocate( RCV(TOTAL) )
#include "mpi_gather_template.inc"

! === RANK 2 (2D Arrays) - e.g. Real dp ===
#define TYPE1 real
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_PRECISION
#define RANK_DECL (:,:)
#define FULL_EXT real_dp_r2
#define ALLOC_GATHER(SND, RCV, PROCS) allocate( RCV(size(SND, 1), size(SND, 2) * PROCS) )
#define ALLOC_GATHERV(SND, RCV, TOTAL) allocate( RCV(size(SND, 1), int(TOTAL / size(SND, 1), kind=i32)) )
#include "mpi_gather_template.inc"

! === RANK 3 (3D Arrays) - e.g. Complex dp ===
#define TYPE1 complex
#define PRECISION1 dp
#define MPIDATATYPE1 MPI_DOUBLE_COMPLEX
#define RANK_DECL (:,:,:)
#define FULL_EXT complex_dp_r3
#define ALLOC_GATHER(SND, RCV, PROCS) allocate( RCV(size(SND, 1), size(SND, 2), size(SND, 3) * PROCS) )
#define ALLOC_GATHERV(SND, RCV, TOTAL) allocate( RCV(size(SND, 1), size(SND, 2), int(TOTAL / (size(SND, 1) * size(SND, 2)), kind=i32)) )
#include "mpi_gather_template.inc"

end module