!> OpenMP utility wrapper module
module mod_omp_utils
#ifdef USEOMP
  use omp_lib
#endif

  implicit none
  private

  public :: omp_max_threads
  public :: omp_num_threads
  public :: omp_thread_num
  public :: omp_num_procs

contains

  !> Get OMP number of maximal available threads, i.e. the value set by `export OMP_NUM_THREADS`
  !> serial default value is 1
  integer function omp_max_threads()
    implicit none
#ifdef USEOMP
    omp_max_threads = omp_get_max_threads()
#else
    omp_max_threads = 1
#endif
  end function omp_max_threads


  !> Get OMP number of currently active threads. Only call from within a parallel region.
  !> serial default value is 1
  integer function omp_num_threads()
    implicit none
#ifdef USEOMP
    omp_num_threads = omp_get_num_threads()
#else
    omp_num_threads = 1
#endif
  end function omp_num_threads


  !> Get current OMP thread index, can be from 0 to N - 1 if omp_max_threads = N
  !> Master thread has always index zero
  !> serial default is zero
  integer function omp_thread_num()
    implicit none
#ifdef USEOMP
    omp_thread_num = omp_get_thread_num()
#else
    omp_thread_num = 0
#endif
  end function omp_thread_num


  !> Get currently available hardware-threads
  !> serial default value is 1
  integer function omp_num_procs()
    implicit none
#ifdef USEOMP
    omp_num_procs = omp_get_num_procs()
#else
    omp_num_procs = 1
#endif
  end function omp_num_procs

end module mod_omp_utils
