module m_hesolver
  use modmpi

  implicit none

  contains

    !BOP
    ! !ROUTINE: hesolver
    ! !INTERFACE:
    subroutine hesolver(hemat, eval, evec, i1, i2, v1, v2, found)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   complex(8) :: hemat(:,:)       ! Upper triangular part of a hermitian matrix
    !   integer(4), optional :: i1, i2 ! Index range of eigen solutions 
    !   real(8), optional :: v1, v2    ! Range of eigenvalues to search for
    ! OUT:
    !   real(8) :: eval(:)                ! Real valued eigenvalues in ascending order 
    !   complex(8), optional :: evec(:,:) ! Corresponding eigenvectors as columns
    !   integer(4), optional :: found     ! Number of found eigen-solutions
    !
    ! !DESCRIPTION:
    !   Wrapper for Lapack's zheevr routine.
    !   Takes a upper triangular part of an hermitian matrix and finds
    !   eigenvalues and/or eigenvectors using the lapack routine {\tt zheevr}.
    ! 
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      complex(8), intent(in) :: hemat(:,:)
      real(8), intent(out) :: eval(:)
      ! Optional arguments
      complex(8), intent(out), optional :: evec(:,:)
      integer(4), intent(in), optional :: i1, i2
      real(8), intent(in), optional :: v1, v2
      integer(4), intent(out), optional :: found

      ! Local variables
      character(1) :: jobzchar, rangechar
      integer(4) :: solsize, hematsize
      real(8) :: vl, vu, abstol
      integer(4) :: il, iu, lwork, info, lrwork, liwork, nsol
      complex(8) :: evecdummy(3,3)

      ! Allocatable arrays
      complex(8), allocatable :: work(:)
      real(8), allocatable :: rwork(:)
      integer, allocatable :: iwork(:), isuppz(:)

      ! External functions
      real(8), external :: dlamch

      ! Tolerance parameter
      !abstol = 2.d0 * dlamch('s')
      abstol = dlamch('Safe minimum')

      hematsize = size(hemat,1)

      ! Input checking
      if(hematsize /= size(hemat,2)) then 
        write(*,*) "Error (hermitian_solver): Matrix needs to be square."
        call terminate
      end if
      if(size(eval) /= hematsize) then 
        write(*,*) "Error (hermitian_solver): Dimension mismatch between eval and mat."
        call terminate
      end if
      if(present(evec)) then
        if(size(evec,1) /= hematsize) then 
          write(*,*) "Error (hermitian_solver): Dimension mismatch between evec and mat."
          call terminate
        end if
      end if
      if((present(i1) .or. present(i2)) .and. (present(v1) .or. present(v2))) then
        write(*,*) "Error (hermitian_solver): I and V specified."
        call terminate
      end if
      if(present(v1) .and. .not. present(v2)) then 
        write(*,*) "Error (hermitian_solver): Specify whole interval."
        call terminate
      end if

      !! RANGES:
      ! Type 'I'
      if(present(i1) .or. present(i2)) then
        il = 1
        iu = hematsize
        rangechar = 'I'
        if(present(i1)) il = i1
        if(present(i2)) iu = i2
        ! Sanity check
        if(il < 1 .or. iu < 1) then 
          write(*,*) "Error (hermitian_solver): iu and il need to be positive."
          call terminate
        end if
        if(il > iu) then 
          write(*,*) "Error (hermitian_solver): il > iu."
          call terminate
        end if
        if(iu > hematsize) then 
          write(*,*) "Error (hermitian_solver): iu-il+1 > nrows."
          write(*,*) "range:", il, iu
          write(*,*) "hematsize:", hematsize
          call terminate
        end if
        allocate(isuppz(2*max(1,iu-il+1)))
      end if
      ! Type 'V'
      if(present(v1) .and. present(v2)) then
        rangechar = 'V'
        vl = v1
        vu = v2
        if(vl > vu) then
          write(*,*) "Error (hermitian_solver): vl > vu"
          write(*,*) "range:", vl, vu
          call terminate
        end if
        allocate(isuppz(2*max(1,hematsize)))
      end if
      ! Type 'A'
      if(.not. (present(i1) .or. present(i2)) .and. .not. present(v1)) then
        rangechar = 'A'
        allocate(isuppz(2*max(1,hematsize)))
      end if

      !! JOBZ
      if(present(evec)) then 
        jobzchar = 'V'
      else
        jobzchar = 'N'
      end if
      ! Get number of requested solutions
      if (rangechar=='I') then
        nsol=iu-il+1
      elseif (rangechar=='V') then
        nsol=vu-vl
      else
        nsol=hematsize
      end if
      ! if the full matrix has to be diagonalized, we use zheevd instead of zheevr
      ! zheevr is found to be instable for full-range calculation in the case of 
      ! spin-polarized RPA calculations
      if ((nsol==hematsize) .and. (jobzchar=='V')) then
        ! Get optimal work array lengths
        call workspacequery_zheevd(jobzchar, rangechar)
        allocate(work(lwork), rwork(lrwork), iwork(liwork))
        
        call zheevd(jobzchar,'u', hematsize, hemat, hematsize, eval, work, lwork, rwork,&
          & lrwork, iwork, liwork, info)
        evec(:,:)=hemat(:,:)

        if(info .ne. 0) then
          write(*,*)
          write(*, '("Error (hermitian_solver): zheevd returned non-zero info:", i6)') info
          call errorinspect(info)
          write(*,*)
          call terminate
        end if

        if(present(found)) then 
          found=hematsize
        end if
      else
        ! Get optimal work array lengths
        call workspacequery_zheevr(jobzchar, rangechar)
        allocate(work(lwork), rwork(lrwork), iwork(liwork))
      
        ! Diagonalize
        if(jobzchar == 'N') then 
          call zheevr(jobzchar, rangechar, 'u', hematsize, hemat, hematsize, vl, vu, il, iu, &
            & abstol, solsize, eval, evecdummy, hematsize, isuppz, work, lwork, rwork,&
            & lrwork, iwork, liwork, info)
        else
          call zheevr(jobzchar, rangechar, 'u', hematsize, hemat, hematsize, vl, vu, il, iu, &
            & abstol, solsize, eval, evec, hematsize, isuppz, work, lwork, rwork,&
            & lrwork, iwork, liwork, info)
        end if

        if(info .ne. 0) then
          write(*,*)
          write(*, '("Error (hermitian_solver): zheevr returned non-zero info:", i6)') info
          call errorinspect(info)
          write(*,*)
          call terminate
        end if

        if(present(found)) then 
          found = solsize
        end if
      end if

      deallocate(isuppz)
      deallocate(work, rwork, iwork)

      contains

        subroutine workspacequery_zheevr(jobztype, rangetype)

          character(1), intent(in) :: jobztype, rangetype

          lwork=-1
          lrwork=-1
          liwork=-1

          allocate(work(1), rwork(1), iwork(1))

          if(jobztype == 'N') then 
            call zheevr(jobztype, rangetype, 'u', hematsize, hemat, hematsize, vl, vu, il, iu, &
              & abstol, solsize, eval, evecdummy, hematsize, isuppz, work, lwork, rwork,&
              & lrwork, iwork, liwork, info)
          else
            call zheevr(jobztype, rangetype, 'u', hematsize, hemat, hematsize, vl, vu, il, iu, &
              & abstol, solsize, eval, evec, hematsize, isuppz, work, lwork, rwork,&
              & lrwork, iwork, liwork, info)
          end if

          ! Adjust workspace
          lwork=int(work(1))
          lrwork=int(rwork(1))
          liwork=int(iwork(1))

          deallocate(work, rwork, iwork)

        end subroutine workspacequery_zheevr

        subroutine workspacequery_zheevd(jobztype, rangetype)

          character(1), intent(in) :: jobztype, rangetype
          
          lwork=-1
          lrwork=-1
          liwork=-1
                      
          allocate(work(1), rwork(1), iwork(1))
          
          call zheevd(jobzchar,'u', hematsize, hemat, hematsize, eval, work, lwork, rwork,&
          & lrwork, iwork, liwork, info)
          
          ! Adjust workspace
          lwork=int(work(1))
          lrwork=int(rwork(1))
          liwork=int(iwork(1))

          deallocate(work, rwork, iwork)

        end subroutine workspacequery_zheevd

        subroutine errorinspect(ierror)
          integer(4) :: ierror

          if( ierror < 0) then
            write(*,'("dhesolver Error cause: Invalid input")')
          else 
            write(*,'("dhesolver Error cause: Internal error")')
          end if
        end subroutine errorinspect

    end subroutine hesolver 
    !EOC

end module m_hesolver
