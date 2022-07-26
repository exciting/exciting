!> FORTRAN 95 interface for LAPACK diagonalization routines.
module diagonalization_interfaces
  use precision, only: dp

  interface

      ! Tridiagonal matrix
      subroutine dstedc(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info)
        import :: dp 
        implicit none
        character, intent(in)    :: compz
        integer,   intent(in)    :: n
        real(dp),  intent(inout) :: d(*)
        real(dp),  intent(inout) :: e(*)
        integer,   intent(in)    :: ldz
        real(dp),  intent(inout) :: z(ldz, *)
        real(dp),  intent(out)   :: work(*)
        integer,   intent(in)    :: lwork
        integer,   intent(out)   :: iwork(*)
        integer,   intent(in)    :: liwork
        integer,   intent(out)   :: info
      end subroutine dstedc

  end interface

end module diagonalization_interfaces