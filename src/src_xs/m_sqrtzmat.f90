module m_sqrtzmat

  use modmpi
  use modscl
  use m_dhesolver
  use m_dzmatmult

  implicit none

  private

  complex(8), parameter :: zone = (1.0d0, 0.0d0)
  complex(8), parameter :: zzero = (0.0d0, 0.0d0)

  public :: sqrtzmat_hepd, sqrtdzmat_hepd

  contains

    subroutine sqrtzmat_hepd(hepdmat)

      complex(8), intent(inout) :: hepdmat(:,:)

      integer(4) :: m, n, i, j
      complex(8), allocatable :: evecs(:,:)
      real(8), allocatable :: evals(:)

      m = size(hepdmat,1)
      n = size(hepdmat,2)

      if(m /= n) then 
        write(*,*) "Error(sqrtmat_hepd): m /= n"
        call terminate
      end if

      allocate(evecs(m,m))
      allocate(evals(m))
      evals = 0.0d0

      ! Diagonalize hermitian matrix
      call hesolver(hepdmat, evals, evec=evecs)

      ! Take square root of eigenvalues
      if(any(evals < 0.0d0)) then 
        write(*,*) "Error(sqrtmat_hepd): Matrix is not positive definit"
        write(*,'(E10.3)') evals
        call terminate
      else
        ! Take D^{1/4} for: A^1/2 = Q D^1/2 Q^H = (Q D^1/4)*(Q D^1/4)^H
        evals = sqrt(sqrt(evals))
      end if

      ! Make Q D^1/4
      do j = 1, m
        do i = 1, m
          evecs(i,j) = evecs(i,j)*evals(j)
        end do
      end do
      
      ! Construct square root matrix
      call zgemm('N','C', m, m, m, zone, evecs, m,&
        & evecs, m, zzero, hepdmat, m)

      deallocate(evecs)
      deallocate(evals)

    end subroutine sqrtzmat_hepd

    subroutine sqrtdzmat_hepd(hepdmat, binfo, eecs)

      type(dzmat), intent(inout) :: hepdmat
      type(blacsinfo), intent(in) :: binfo
      integer(4), intent(in), optional :: eecs

      integer(4) :: m, n, i, j, jg, clustersize
      type(dzmat) :: evecs
      real(8), allocatable :: evals(:)

      character(*), parameter :: thisname = "sqrtdzmat_hepd"
      
      if(present(eecs)) then
        clustersize = eecs
      else
        clustersize = 3
      end if

      if(hepdmat%isdistributed .eqv. .false.) then
        call sqrtzmat_hepd(hepdmat%za)
        return
      end if

#ifdef SCAL
      if(hepdmat%context /= binfo%context) then
        if(binfo%isroot) then 
          write(*,*) "Error(sqrtdzmat_hepd): Wrong BLACS context."
          write(*,*) "  mat%context =", hepdmat%context
          write(*,*) "  binfo%context =", binfo%context
        end if
        call blacs_barrier(binfo%context, "A")
        call terminate
      end if

      m = hepdmat%nrows
      n = hepdmat%ncols

      if(m /= n) then 
        if(binfo%isroot) then 
          write(*,*) "Error(sqrtdzmat_hepd): m /= n"
        end if
        call blacs_barrier(binfo%context, "A")
        call terminate
      end if

      allocate(evals(m))
      evals = 0.0d0

      ! Diagonalize hermitian matrix
      call new_dzmat(evecs, m, n, binfo, hepdmat%mblck, hepdmat%nblck)
      call dhesolver(hepdmat, evals, binfo, evec=evecs, eecs=clustersize)

      !if(binfo%isroot) then 
      !  write(*,'("Info(",a,"): Passed diagonalizaion")') trim(thisname)
      !end if

      ! Take square root of eigenvalues
      if(any(evals < 0.0d0)) then 
        if(binfo%isroot) then 
          write(*,'("Error(",a,"):&
            & Matrix is not positive definit")') trim(thisname)
          write(*,'(E10.3)') evals
        end if
        call blacs_barrier(binfo%context, "A")
        call terminate
      else
        ! Take D^{1/4} for: A^1/2 = Q D^1/2 Q^H = (Q D^1/4)*(Q D^1/4)^H
        evals = sqrt(sqrt(evals))
      end if

      ! Make Q D^1/4
      do j = 1, evecs%ncols_loc
        jg = evecs%c2g(j)
        do i = 1, evecs%nrows_loc
          evecs%za(i,j) = evecs%za(i,j)*evals(jg)
        end do
      end do

      !if(binfo%isroot) then 
      !  write(*,'("Info(",a,"): Q*D^1/4 passed")') trim(thisname)
      !end if

      !write(*,'("Info(",a," at rank ", i3"): evecs descriptor = ", 9i5)')&
      !  & trim(thisname), mpiglobal%rank, evecs%desc
      !write(*,'("Info(",a," at rank ", i3"): hepdmat descriptor = ", 9i5)')&
      !  & trim(thisname), mpiglobal%rank, hepdmat%desc
      
      ! Construct square root matrix A^1/2 = (Q D^1/4)*( Q D^1/4)^H
      call dzmatmult(evecs, evecs, hepdmat, transb='C')

      !write(*,'("Info(",a," at rank ", i3"): Passed matrix mult.")')&
      !  & trim(thisname), mpiglobal%rank

      !if(binfo%isroot) then 
      !  write(*,'("Info(",a,"): Matrix mult. passed")') trim(thisname)
      !end if

      deallocate(evals)
      call del_dzmat(evecs)
#else
      write(*,'("Error(",a,"): -DSCAL was not specified but matrix is distributed.")')&
        & trim(thisname)
      call terminate
#endif

    end subroutine sqrtdzmat_hepd

end module m_sqrtzmat
