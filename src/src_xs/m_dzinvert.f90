module m_dzinvert
  use modmpi 
  use modscl

  implicit none

  contains

    !BOP
    ! !ROUTINE: dhesolver
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

      integer(4) :: info, lwork
      integer(4), allocatable :: ipiv(:)
      complex(8), allocatable :: work(:)
#ifdef SCAL
      integer(4) :: ia, ja, liwork
      integer(4), allocatable :: iwork(:)
#endif

#ifdef SCAL
      allocate(ipiv(zmat%nrows_loc+zmat%mblck))
#else
      allocate(ipiv(zmat%nrows))
#endif

      !! LU factorization of general complex matrix A
#ifdef SCAL
      ia = 1
      ja = 1
      call pzgetrf(zmat%nrows, zmat%ncols, zmat%za, ia, ja, zmat%desc,&
        & ipiv, info)
      if(info /= 0) then
        write(*,'("dzinvert (ERROR): PZGETRF has returned non-zero info.")')
        if(info < 0) then
          write(*,'("  PZGETRF: Incorrect input.")')
        else
          write(*,'("  PZGETRF: U is exactly singular.")')
        end if
        call terminate
      end if
#else
      call zgetrf(zmat%nrows, zmat%ncols, zmat%za, zmat%nrows, ipiv, info)
      if(info /= 0) then
        write(*,'("dzinvert (ERROR): ZGETRF has returned non-zero info.")')
        if(info < 0) then
          write(*,'("  ZGETRF: Incorrect input.")')
        else
          write(*,'("  ZGETRF: U is exactly singular.")')
        end if
        call terminate
      end if
#endif
      
      !! Compute Inverse with LU factorized A

      ! Workspace query 
      lwork = -1
      allocate(work(3))
#ifdef SCAL
      liwork = -1
      allocate(iwork(3))
      call pzgetri(zmat%nrows, zmat%za, ia, ja, zmat%desc, ipiv,&
        & work, lwork, iwork, liwork, info)
      lwork = int(work(1))
      liwork = iwork(1)
      deallocate(work,iwork)
      allocate(work(lwork),iwork(liwork))
#else
      call zgetri(zmat%nrows, zmat%za, zmat%nrows, ipiv, work, lwork, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
#endif

      ! Computation
#ifdef SCAL
      call pzgetri(zmat%nrows, zmat%za, ia, ja, zmat%desc, ipiv,&
        & work, lwork, iwork, liwork, info)
      if(info /= 0) then
        write(*,'("dzinvert (ERROR): PZGETRI has returned non-zero info.")')
        if(info < 0) then
          write(*,'("  PZGETRF: Incorrect input.")')
        else
          write(*,'("  PZGETRF: Matrix singular.")')
        end if
        call terminate
      end if
      deallocate(work,iwork)
#else
      call zgetri(zmat%nrows, zmat%za, zmat%nrows, ipiv, work, lwork, info)
      if(info /= 0) then
        write(*,'("dzinvert (ERROR): ZGETRI has returned non-zero info.")')
        if(info < 0) then
          write(*,'("  ZGETRF: Incorrect input.")')
        else
          write(*,'("  ZGETRF: Matrix singular.")')
        end if
        call terminate
      end if
      deallocate(work)
#endif
    end subroutine dzinvert
    !EOC
    
end module m_dzinvert
