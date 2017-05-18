module m_invertzmat
  use modmpi 
  use modscl

  implicit none

  private

  public :: zinvert, dzinvert

  contains

    !BOP
    ! !ROUTINE: zinvert
    ! !INTERFACE:
    subroutine zinvert(zmat)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN/OUT:
    !   complex(8) :: zmat  ! On entry: Matrix.
    !                       ! On exit : Inverted matrix.
    !
    ! !DESCRIPTION:
    !   Inverts complex matrix.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      ! Arguments
      complex(8), intent(inout) :: zmat(:,:)

      integer(4) :: m, n
      integer(4) :: info, lwork
      integer(4), allocatable :: ipiv(:)
      complex(8), allocatable :: work(:)

      m = size(zmat,1)
      n = size(zmat,2)

      if(m /= n) then 
        write(*,*) "Error(zivert): m /= n"
        call terminate
      end if

      allocate(ipiv(m))

      ! LU factorization of general complex matrix A
      call zgetrf(m, n, zmat, m, ipiv, info)
      if(info /= 0) then
        write(*,'("Error(zinvert): ZGETRF has returned non-zero info.")')
        if(info < 0) then
          write(*,'("  ZGETRF: Incorrect input.")')
        else
          write(*,'("  ZGETRF: U is exactly singular.")')
        end if
        call terminate
      end if
      
      ! Compute Inverse with LU factorized A

      ! Workspace query 
      lwork = -1
      allocate(work(3))
      call zgetri(m, zmat, m, ipiv, work, lwork, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))

      ! Computation
      call zgetri(m, zmat, m, ipiv, work, lwork, info)
      if(info /= 0) then
        write(*,'("Error(zinvert): ZGETRI has returned non-zero info.")')
        if(info < 0) then
          write(*,'("  ZGETRF: Incorrect input.")')
        else
          write(*,'("  ZGETRF: Matrix singular.")')
        end if
        call terminate
      end if
      deallocate(work)
      deallocate(ipiv)

    end subroutine zinvert
    !EOC

    !BOP
    ! !ROUTINE: dzinvert
    ! !INTERFACE:
    subroutine dzinvert(zmat)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN/OUT:
    !   type(dzmat) :: zmat ! On entry: Distributed matrix to invert.
    !                       ! On exit : Inverted matrix.
    !
    ! !DESCRIPTION:
    !   Takes a distributed general complex matrix and inverts it.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      type(dzmat), intent(inout) :: zmat

      integer(4) :: info
      integer(4) :: ia, ja
      integer(4) :: lwork, liwork
      integer(4), allocatable :: ipiv(:)
      integer(4), allocatable :: iwork(:)
      complex(8), allocatable :: work(:)

      if(zmat%isdistributed .eqv. .false.) then 
        call zinvert(zmat%za)
        return
      end if

#ifdef SCAL
      allocate(ipiv(zmat%nrows_loc+zmat%mblck))

      ! LU factorization of general complex matrix A
      ia = 1
      ja = 1
      call pzgetrf(zmat%nrows, zmat%ncols, zmat%za, ia, ja, zmat%desc,&
        & ipiv, info)
      if(info /= 0) then
        write(*,'("Error(dzinvert): PZGETRF has returned non-zero info.")')
        if(info < 0) then
          write(*,'("  PZGETRF: Incorrect input.")')
        else
          write(*,'("  PZGETRF: U is exactly singular.")')
        end if
        call terminate
      end if
      
      ! Compute Inverse with LU factorized A

      ! Workspace query 
      liwork = -1
      allocate(iwork(3))
      allocate(work(3))
      call pzgetri(zmat%nrows, zmat%za, ia, ja, zmat%desc, ipiv,&
        & work, lwork, iwork, liwork, info)
      lwork = int(work(1))
      liwork = iwork(1)
      deallocate(work,iwork)
      allocate(work(lwork),iwork(liwork))

      ! Computation
      call pzgetri(zmat%nrows, zmat%za, ia, ja, zmat%desc, ipiv,&
        & work, lwork, iwork, liwork, info)
      if(info /= 0) then
        write(*,'("Error(dzinvert): PZGETRI has returned non-zero info.")')
        if(info < 0) then
          write(*,'("  PZGETRF: Incorrect input.")')
        else
          write(*,'("  PZGETRF: Matrix singular.")')
        end if
        call terminate
      end if
      deallocate(work,iwork)
      deallocate(ipiv)
#endif
    end subroutine dzinvert
    !EOC

end module m_invertzmat
